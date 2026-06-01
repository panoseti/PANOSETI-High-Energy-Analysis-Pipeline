#!/usr/bin/env python3

# pcapdecoder.py: A pure Python decoder for PANOSETI data stored in PCAP format.
# Does not rely on any external packages.

# Author: Stephen Fegan <sfegan@llr.in2p3.fr> (2026-05-10)
# Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

import struct
from collections import namedtuple
import os
import glob
import datetime
import bisect

# Define the structure for a PANOSETI Science Packet
# Using namedtuple for a "struct-like" behavior as requested.
SciencePacket = namedtuple('SciencePacket', [
    'filename', 'file_offset', 'file_captured_packet_index', 'file_science_packet_index',
    'pcap_time', 'pcap_sec', 'pcap_nsec', 'incl_len', 'sender_ip',
    'acq_mode', 'packet_ver', 'packet_num', 'board_loc',
    'telescope_id', 'quabo_id', 'tai', 'nanosec', 
    'event_time', 'event_time_sec', 'event_time_nsec', 'event_time_good',
    'flags', 'gti_index', 'gti_event_time', 'pix_data', 
])

def parse_time(t):
    """
    Attempts to parse a time value into a Unix timestamp.
    Accepts floats/ints (Unix epoch) or strings (various ISO-like formats).
    Assumes UTC if no timezone is specified.
    """
    if isinstance(t, (int, float)):
        return float(t)
    if isinstance(t, str):
        # Try a few common formats
        for fmt in ["%Y-%m-%dT%H:%M:%S.%f", "%Y-%m-%d %H:%M:%S.%f", 
                    "%Y-%m-%dT%H:%M:%S", "%Y-%m-%d %H:%M:%S",
                    "%Y-%m-%d"]:
            try:
                dt = datetime.datetime.strptime(t, fmt)
                if dt.tzinfo is None:
                    # Assume UTC as requested
                    dt = dt.replace(tzinfo=datetime.timezone.utc)
                return dt.timestamp()
            except ValueError:
                continue
        # Last resort: fromisoformat (available in Python 3.7+)
        try:
            dt = datetime.datetime.fromisoformat(t)
            if dt.tzinfo is None:
                dt = dt.replace(tzinfo=datetime.timezone.utc)
            return dt.timestamp()
        except (ValueError, AttributeError):
            pass
    raise ValueError(f"Could not parse time: {t}")

class GTIFilter:
    """
    Manages Good Time Intervals (GTI) for filtering packets.
    Normalizes, merges overlapping intervals, and provides efficient lookup.
    """
    def __init__(self, gti):
        if isinstance(gti, GTIFilter):
            self.intervals = gti.intervals
            self._starts = gti._starts
            self._current_idx = gti._current_idx
            return

        if gti is None:
            self.intervals = [(-float('inf'), float('inf'), True, 0)]
            self._starts = [-float('inf')]
            self._current_idx = 0
            return

        if isinstance(gti, dict):
            gti = [gti]
        
        # 1. Normalize and collect start/stop pairs
        raw_intervals = []
        for d in gti:
            start = parse_time(d.get('start', -float('inf')))
            # Support both 'stop' and 'end'
            stop = parse_time(d.get('stop') or d.get('end', float('inf')))
            if start < stop:
                raw_intervals.append((start, stop))
        
        # 2. Sort by start time
        raw_intervals.sort()
        
        # 3. Merge overlapping intervals
        merged = []
        if raw_intervals:
            curr_start, curr_stop = raw_intervals[0]
            for next_start, next_stop in raw_intervals[1:]:
                if next_start <= curr_stop:
                    curr_stop = max(curr_stop, next_stop)
                else:
                    merged.append((curr_start, curr_stop))
                    curr_start, curr_stop = next_start, next_stop
            merged.append((curr_start, curr_stop))
        
        # 4. Build a continuous set of intervals from -inf to inf
        self.intervals = []
        last_t = -float('inf')
        gti_idx = 0
        for start, stop in merged:
            if start > last_t:
                self.intervals.append((last_t, start, False, -1))
            self.intervals.append((start, stop, True, gti_idx))
            gti_idx += 1
            last_t = stop
        if last_t < float('inf'):
            self.intervals.append((last_t, float('inf'), False, -1))
        
        if not merged:
            self.intervals = [(-float('inf'), float('inf'), True, 0)]
            
        self._starts = [i[0] for i in self.intervals]
        self._current_idx = 0

    def get_info(self, t):
        """Returns (is_good, gti_index) for time t."""
        # Check cached interval first (optimizes for time-ordered data)
        start, stop, good, idx = self.intervals[self._current_idx]
        if start <= t < stop:
            return good, idx, start
        
        # Find correct interval using binary search
        idx_bin = bisect.bisect_right(self._starts, t) - 1
        self._current_idx = idx_bin
        start, _, good, idx = self.intervals[idx_bin]
        return good, idx, start

def quabo_timestamp_good(pcap_sec):
    return False # PANOSETI board times no good for now

def wr_to_unix(pkt_tai, pkt_nsec, pcap_sec, ignore_clock_desync=False):
    """
    Given a WR packet time (TAI) with only 10 bits of sec,
    and a Unix time that's within a few ms,
    return the complete WR time (in Unix time, not TAI).

    See: https://github.com/panoseti/panoseti/blob/master/control/utils/pff.py#L238
    """
    # 37 is the TAI-UTC offset. 
    # The packet TAI seconds is only 10 bits (0-1023).
    d = (pcap_sec - pkt_tai + 37) % 1024
    if d == 0:
        return True, pcap_sec, pkt_nsec, pcap_sec + pkt_nsec / 1e9
    elif d == 1:
        return True, pcap_sec - 1, pkt_nsec, pcap_sec - 1 + pkt_nsec / 1e9
    elif d == 1023:
        return True, pcap_sec + 1, pkt_nsec, pcap_sec + 1 + pkt_nsec / 1e9
    else:
        # The WR and DAQ clocks differ by > 1s => out of sync
        # Return 0 if ignore_clock_desync is False. Otherwise, return an approximation to the time.
        if ignore_clock_desync:
            approx_t = pcap_sec + pkt_nsec / 1e9
            return False, pcap_sec, pkt_nsec, approx_t
        else:
            raise Exception('WR and Unix times differ by > 1 sec: pkt_tai %d pcap_sec %d d %d' % (pkt_tai, pcap_sec, d))

class PcapDecoder:
    """
    A pure Python decoder for PANOSETI data stored in PCAP format.
    Does not rely on any external packages.

    See: Jamie's C++ code
    See: https://github.com/panoseti/panoseti/wiki/Quabo-packet-interface#science-packets
    """
    def __init__(self, filename, dest_port=60001):
        self.filename = filename
        self._file = open(filename, 'rb')
        self._is_pcapng = False
        self._endian = '<'
        self._ts_resol = 1000000 # Default to microseconds (10^-6)
        self.dest_port = dest_port
        self.file_captured_packet_index = 0
        self.file_science_packet_index = 0
        self._detect_format()

    def _detect_format(self):
        magic = self._file.read(4)
        if len(magic) < 4:
            raise EOFError("File too small to be a PCAP/PCAPNG file")
        
        # Standard PCAP
        if magic == b'\xa1\xb2\xc3\xd4':
            self._endian = '>'
            self._is_pcapng = False
        elif magic == b'\xd4\xc3\xb2\xa1':
            self._endian = '<'
            self._is_pcapng = False
        elif magic == b'\xa1\xb2\x3c\x4d': # PCAP nanosecond resolution
            self._endian = '>'
            self._is_pcapng = False
        elif magic == b'\x4d\x3c\xb2\xa1': # PCAP nanosecond resolution (swapped)
            self._endian = '<'
            self._is_pcapng = False
        # PCAPNG
        elif magic == b'\x0a\x0d\x0d\x0a':
            self._is_pcapng = True
            # PCAPNG magic is always this, but Byte Order Magic (BOM) follows
            self._file.read(4) # Skip Section Length
            bom = self._file.read(4)
            if bom == b'\x1a\x2b\x3c\x4d':
                self._endian = '>'
            elif bom == b'\x4d\x3c\x2b\x1a':
                self._endian = '<'
            else:
                raise ValueError(f"Invalid PCAPNG BOM: {bom.hex()}")
            # Seek back to start of file for easier block processing
            self._file.seek(0)
        else:
            raise ValueError(f"Unknown file magic: {magic.hex()}")

        if not self._is_pcapng:
            # Skip the rest of the 24-byte global header
            self._file.seek(24)

    def __iter__(self):
        if self._is_pcapng:
            return self._iter_pcapng()
        else:
            return self._iter_pcap()

    def _iter_pcap(self):
        while True:
            # Capture the start offset of the packet header
            offset = self._file.tell()
            
            # Read 16-byte packet header
            header_data = self._file.read(16)
            if not header_data:
                break
            if len(header_data) < 16:
                break # Truncated
            
            pcap_sec, pcap_usec, incl_len, orig_len = struct.unpack(f"{self._endian}IIII", header_data)
            packet_data = self._file.read(incl_len)
            
            packet = self._parse_packet(packet_data, pcap_sec, pcap_usec*1000, offset, self.file_captured_packet_index, self.file_science_packet_index)
            if packet:
                yield packet
                self.file_packet_index += 1

            self.file_captured_packet_index += 1

    def _iter_pcapng(self):
        while True:
            # Capture the start offset of the block
            offset = self._file.tell()
            
            block_header = self._file.read(8)
            if not block_header:
                break
            if len(block_header) < 8:
                break
            
            block_type, block_total_length = struct.unpack(f"{self._endian}II", block_header)
            block_data = self._file.read(block_total_length - 8)
            
            # Interface Description Block (IDB) is 0x00000001
            if block_type == 0x00000001:
                # Parse IDB for timestamp resolution (option 9)
                options = block_data[8:-4] # Skip link type, snaplen, and trailing length
                offset_opts = 0
                while offset_opts + 4 <= len(options):
                    ocode, olen = struct.unpack(f"{self._endian}HH", options[offset_opts:offset_opts+4])
                    if ocode == 0: break # End of options
                    if ocode == 9 and olen == 1:
                        resol = options[offset_opts+4]
                        self._ts_resol = 10**resol
                        break
                    offset_opts += 4 + (olen + 3)//4 * 4

            # Enhanced Packet Block (EPB) is 0x00000006
            if block_type == 0x00000006:
                # EPB structure:
                # Interface ID (4)
                # Timestamp (High) (4)
                # Timestamp (Low) (4)
                # Captured Packet Length (4)
                # Original Packet Length (4)
                # Packet Data (...)
                # Padding (to 4-byte boundary)
                # Block Total Length (4) - already read as part of block_data
                
                epb_header = struct.unpack(f"{self._endian}IIIII", block_data[:20])
                if_id, ts_high, ts_low, incl_len, orig_len = epb_header
                
                # Combine timestamps (PCAPNG usually uses microseconds or nanoseconds since epoch)
                timestamp = (ts_high << 32) | ts_low
                pcap_sec = timestamp // self._ts_resol
                # Calculate nsec based on resolution
                if(self._ts_resol == 1000000000): # ns resolution
                    pcap_nsec = timestamp % self._ts_resol
                elif self._ts_resol < 1000000000: # us probably
                    pcap_nsec = (timestamp % self._ts_resol) * (1000000000 // self._ts_resol)
                else: # better than ns reolution
                    pcap_nsec = (timestamp % self._ts_resol) // (self._ts_resol // 1000000000)

                packet_data = block_data[20:20+incl_len]
                packet = self._parse_packet(packet_data, pcap_sec, pcap_nsec, offset, self.file_captured_packet_index, self.file_science_packet_index)
                if packet:
                    yield packet
                    self.file_science_packet_index += 1
                self.file_captured_packet_index += 1

            # Skip other block types (SHB, IDB, etc. are already read into block_data)

    def _parse_packet(self, packet_data, pcap_sec, pcap_nsec, file_offset, file_captured_packet_index, file_science_packet_index):
        # Based on Quabo Packet Interface Wiki:
        # Science packets are sent to port 60001.
        # Total payload sizes: 
        #   - 528 bytes (16-bit modes)
        #   - 272 bytes (8-bit modes)
        # 42 bytes of Eth/IP/UDP headers are expected if captured on the wire.
        
        if len(packet_data) == 570: # 42 + 528
            payload = packet_data[42:]
        elif len(packet_data) == 314: # 42 + 272
            payload = packet_data[42:]
        else:
            # Not a recognized science packet size
            return None
        
        # Ensure it is an IPv4 UDP packet and destination port matches self.dest_port
        if len(packet_data) < 42:
            return None
        
        # Check EtherType is IPv4 (0x0800)
        eth_type = struct.unpack(">H", packet_data[12:14])[0]
        if eth_type != 0x0800:
            return None
        
        # Check IP version is 4, and extract header length
        ip_header_first_byte = packet_data[14]
        ip_version = ip_header_first_byte >> 4
        if ip_version != 4:
            return None
        
        ip_ihl = (ip_header_first_byte & 0x0F) * 4
        # Ensure protocol is UDP (17)
        # The protocol field is at offset 9 in IPv4 header (byte 23 of packet_data)
        ip_proto = packet_data[23]
        if ip_proto != 17:
            return None
        
        # UDP header starts after the IP header
        udp_start = 14 + ip_ihl
        if len(packet_data) < udp_start + 8:
            return None
        
        # Extract UDP destination port (bytes 2-3 of UDP header, big endian)
        dest_port = struct.unpack(">H", packet_data[udp_start+2 : udp_start+4])[0]
        
        if self.dest_port != 0 and dest_port != self.dest_port:
            return None
        
        # PANOSETI metadata (16 bytes)
        # acq_mode(1), packet_ver(1), packet_num(2), board_loc(2), TAI(4), nanosec(4), flags(2)
        # All multi-byte fields in Quabo packets are little-endian.
        meta_format = "<BBHHIIH"
        meta_size = struct.calcsize(meta_format)
        meta = struct.unpack(meta_format, payload[:meta_size])
        
        acq_mode = meta[0]
        
        # Pixel data starts at offset 16
        # Modes:
        # 0x01: Pulse Height (16-bit signed)
        # 0x02: Imaging (16-bit unsigned)
        # 0x06: Imaging (8-bit unsigned)
        
        if acq_mode == 0x01:
            pix_format = "<256h" # 16-bit signed
            pix_bytes = 512
        elif acq_mode == 0x02:
            pix_format = "<256H" # 16-bit unsigned
            pix_bytes = 512
        elif acq_mode == 0x06:
            pix_format = "<256B" # 8-bit unsigned
            pix_bytes = 256
        else:
            # Handle other modes (e.g., dual mode 0x03) if they follow 16-bit format
            # The C++ code used acq_mode 1 as signed, others as unsigned.
            pix_format = "<256H" if len(payload) >= 528 else "<256B"
            pix_bytes = 512 if len(payload) >= 528 else 256
        
        if len(payload) < meta_size + pix_bytes:
            return None

        pix_data = struct.unpack(pix_format, payload[meta_size:meta_size+pix_bytes])
        
        board_loc = meta[3]
        telescope_id = board_loc >> 2
        quabo_id = board_loc & 0x3
        
        tai = meta[4]
        nanosec = meta[5]

        pcap_time = pcap_sec + pcap_nsec / 1e9

        event_time_good = False
        if quabo_timestamp_good(pcap_sec):
            event_time_good, event_time_sec, event_time_nsec, event_time = wr_to_unix(tai, nanosec, pcap_sec, ignore_clock_desync=True)

        if not event_time_good:
            # Fall back to using PCAP time if WR time is not good
            event_time_sec = pcap_sec
            event_time_nsec = pcap_nsec
            event_time = pcap_time

        sender_ip = ".".join(str(b) for b in packet_data[26:30])

        return SciencePacket(
            filename=self.filename,
            file_offset=file_offset,
            file_captured_packet_index=file_captured_packet_index,
            file_science_packet_index=self.file_science_packet_index,
            pcap_time=pcap_time,
            pcap_sec=pcap_sec,
            pcap_nsec=pcap_nsec,
            incl_len=len(packet_data),
            acq_mode=acq_mode,
            packet_ver=meta[1],
            packet_num=meta[2],
            board_loc=board_loc,
            telescope_id=telescope_id,
            quabo_id=quabo_id,
            tai=tai,
            nanosec=nanosec,
            event_time_sec=event_time_sec,
            event_time_nsec=event_time_nsec,
            event_time=event_time,
            event_time_good=event_time_good,
            flags=meta[6],
            gti_index=0,
            gti_event_time=event_time,
            pix_data=pix_data,
            sender_ip=sender_ip
        )

    def close(self):
        self._file.close()

def get_panoseti_packets(filenames, gtis=None, verbose=False, dest_port=60001):
    """
    Convenience function to get packets from one or more PANOSETI PCAP/PCAPNG files.
    'filenames' can be a single filename (string), a glob pattern (string), 
    or a list/iterable of filenames.
    'gtis' can be a GTIFilter object, a dict with start/stop keys, or a list of such dicts.
    """
    gti_filter = gtis if isinstance(gtis, GTIFilter) else GTIFilter(gtis)
    if isinstance(filenames, str):
        if any(char in filenames for char in '*?['):
            files = sorted(glob.glob(filenames))
        else:
            files = [filenames]
    else:
        files = []
        for item in filenames:
            if isinstance(item, str) and any(char in item for char in '*?['):
                files.extend(sorted(glob.glob(item)))
            else:
                files.append(item)

    for filename in sorted(files):
        if(verbose):
            print(f"Processing file: {filename} (size: {os.path.getsize(filename)} bytes)")
        decoder = PcapDecoder(filename, dest_port=dest_port)
        try:
            for packet in decoder:
                is_good, gti_idx, gti_start_time = gti_filter.get_info(packet.event_time)
                if is_good:
                    gti_event_time = packet.event_time
                    if gti_start_time != -float('inf'):
                        gti_event_time -= gti_start_time
                    yield packet._replace(gti_index=gti_idx, gti_event_time=gti_event_time)
        finally:
            decoder.close()
