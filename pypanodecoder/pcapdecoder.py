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
_SciencePacketBase = namedtuple('SciencePacket', [
    'filename', 'file_offset', 'file_captured_packet_index', 'file_science_packet_index',
    'pcap_time', 'pcap_sec', 'pcap_nsec', 'pcap_len', 'pcap_src_ip',
    'acq_mode', 'packet_ver', 'packet_num', 'board_loc',
    'telescope_id', 'quabo_id', 'tai', 'nanosec',
    'event_time', 'event_time_sec', 'event_time_nsec', 'event_time_good',
    'flags', 'gti_index', 'gti_pcap_time', 'pix_data',
])

class SciencePacket(_SciencePacketBase):
    def __repr__(self):
        lines = [f"{type(self).__name__}("]
        for field in self._fields:
            value = getattr(self, field)
            if field == 'pix_data':
                lines.append(f"  {field}: <{len(value)} pixels>,")
            elif field in ['pcap_time', 'event_time', 'gti_pcap_time']:
                lines.append(f"  {field}: {value} (ns),")
            else:
                lines.append(f"  {field}: {value!r},")
        lines.append(")")
        return "\n".join(lines)

# Define the structure for a PANOSETI House Keeping Packet
_HKPacketBase = namedtuple('HKPacket', [
    'filename', 'file_offset', 'file_captured_packet_index', 'file_hk_packet_index',
    'pcap_time', 'pcap_sec', 'pcap_nsec', 'pcap_len', 'pcap_src_ip',
    'board_loc', 'telescope_id', 'quabo_id',
    'gti_index', 'gti_pcap_time', 'data'
])

class HKPacket(_HKPacketBase):
    def __repr__(self):
        lines = [f"{type(self).__name__}("]
        for field in self._fields:
            value = getattr(self, field)
            if field == 'data':
                lines.append(f"  {field}: <{len(value)} bytes>,")
            elif field in ['pcap_time', 'gti_pcap_time']:
                lines.append(f"  {field}: {value} (ns),")
            else:
                lines.append(f"  {field}: {value!r},")
        lines.append(")")
        return "\n".join(lines)

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

class QuaboTimeUnroller:
    """
    Handles unrolling and calibration of the 10-bit (1024s) Quabo TAI counter.
    Calculates a per-quabo offset to align hardware time with PCAP arrival time
    while maintaining hardware precision and monotonicity.
    """
    def __init__(self, rollover_sec=1024):
        self.rollover_sec = rollover_sec
        self.offset_ns = None
        self.last_tai = None
        self.last_pcap_ns = None

    def unroll(self, tai, nsec, pcap_ns):
        """
        Calculates monotonic nanoseconds since epoch for a packet.
        
        Args:
            tai (int): 10-bit TAI second counter from packet.
            nsec (int): Nanosecond counter from packet.
            pcap_ns (int): PCAP arrival time in nanoseconds.
        """
        # Hardware time within the current 1024s cycle
        hw_cycle_ns = int(tai) * 1000000000 + int(nsec)

        # 1. Check for gap or initial calibration
        if (self.offset_ns is None or 
            (self.last_pcap_ns is not None and abs(pcap_ns - self.last_pcap_ns) > self.rollover_sec * 1000000000)):
            # Calibrate: align this hardware cycle to the current PCAP second.
            # We align to the start of the PCAP second to keep it simple.
            pcap_sec_ns = (pcap_ns // 1000000000) * 1000000000
            self.offset_ns = pcap_sec_ns - (int(tai) * 1000000000)
            # Alternatively, to be closer to "middle of the second" as suggested:
            # self.offset_ns = (pcap_ns // 1000000000) * 1000000000 - (int(tai) * 1000000000)
            # This ensures event_time_sec is initially close to pcap_sec.

        # 2. Detect Rollover
        # If TAI dropped significantly (e.g. 1023 -> 0), increment offset
        if self.last_tai is not None and tai < self.last_tai - (self.rollover_sec // 2):
            self.offset_ns += self.rollover_sec * 1000000000
        # Handle case where packets might arrive slightly out of order near rollover
        elif self.last_tai is not None and tai > self.last_tai + (self.rollover_sec // 2):
            self.offset_ns -= self.rollover_sec * 1000000000

        self.last_tai = tai
        self.last_pcap_ns = pcap_ns
        
        return hw_cycle_ns + self.offset_ns

def quabo_timestamp_good(pcap_sec):
    return False # PANOSETI board times no good for now

class PcapDecoder:
    """
    A pure Python decoder for PANOSETI data stored in PCAP format.
    Does not rely on any external packages.

    See: Jamie's C++ code
    See: https://github.com/panoseti/panoseti/wiki/Quabo-packet-interface#science-packets
    """
    def __init__(self, filename, science_port=60001, hk_port=60000, decode_hk=False, acq_mode=0x01):
        self.filename = filename
        self._file = open(filename, 'rb')
        self._is_pcapng = False
        self._endian = '<'
        self._ts_resol = 1000000 # Default to microseconds (10^-6)
        self.science_port = science_port
        self.hk_port = hk_port
        self.decode_hk = decode_hk
        if acq_mode is None:
            self._acq_mode_filter = None
        elif isinstance(acq_mode, int):
            self._acq_mode_filter = {acq_mode}
        else:
            try:
                self._acq_mode_filter = set(acq_mode)
            except TypeError:
                self._acq_mode_filter = {acq_mode}
        self.file_captured_packet_index = 0
        self.file_science_packet_index = 0
        self.file_hk_packet_index = 0
        self._time_unrollers = {} # board_loc -> QuaboTimeUnroller
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

            packet = self._parse_packet(packet_data, pcap_sec, pcap_usec*1000, offset, self.file_captured_packet_index, self.file_science_packet_index, self.file_hk_packet_index)
            if packet:
                yield packet
                if isinstance(packet, SciencePacket):
                    self.file_science_packet_index += 1
                elif isinstance(packet, HKPacket):
                    self.file_hk_packet_index += 1

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
                packet = self._parse_packet(packet_data, pcap_sec, pcap_nsec, offset, self.file_captured_packet_index, self.file_science_packet_index, self.file_hk_packet_index)
                if packet:
                    yield packet
                    if isinstance(packet, SciencePacket):
                        self.file_science_packet_index += 1
                    elif isinstance(packet, HKPacket):
                        self.file_hk_packet_index += 1
                self.file_captured_packet_index += 1

            # Skip other block types (SHB, IDB, etc. are already read into block_data)

    def _parse_packet(self, packet_data, pcap_sec, pcap_nsec, file_offset, file_captured_packet_index, file_science_packet_index, file_hk_packet_index):
        # Based on Quabo Packet Interface Wiki:
        # Science packets are sent to port 60001 (by default).
        # House Keeping packets are sent to port 60000 (by default).

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
        ip_proto = packet_data[23]
        if ip_proto != 17:
            return None

        # UDP header starts after the IP header
        udp_start = 14 + ip_ihl
        if len(packet_data) < udp_start + 8:
            return None

        # Extract UDP destination port (bytes 2-3 of UDP header, big endian)
        dest_port = struct.unpack(">H", packet_data[udp_start+2 : udp_start+4])[0]

        sender_ip_bytes = packet_data[26:30]
        sender_ip = ".".join(str(b) for b in sender_ip_bytes)
        payload = packet_data[udp_start+8:]

        # Handle House Keeping packets
        if self.decode_hk and dest_port == self.hk_port:
            return self._decode_hk_packet(payload, pcap_sec, pcap_nsec, sender_ip_bytes, sender_ip, file_offset, file_captured_packet_index, file_hk_packet_index, len(packet_data))

        # Handle Science packets
        if dest_port == self.science_port:
            return self._decode_science_packet(payload, pcap_sec, pcap_nsec, sender_ip, file_offset, file_captured_packet_index, file_science_packet_index, len(packet_data))

        return None

    def _decode_hk_packet(self, payload, pcap_sec, pcap_nsec, sender_ip_bytes, sender_ip, file_offset, file_captured_packet_index, file_hk_packet_index, incl_len):
        # Decode board_loc, telescope_id, and quabo_id from source IP address
        # 2 LSBs are quabo_id and next 7 LSBs are telescope_id
        # All 9 of these together are board_loc
        octet3 = sender_ip_bytes[2]
        octet4 = sender_ip_bytes[3]
        board_loc = ((octet3 & 0x02) << 8) | octet4
        quabo_id = board_loc & 0x03
        telescope_id = (board_loc >> 2) & 0xff

        pcap_time_ns = int(pcap_sec) * 1000000000 + int(pcap_nsec)

        return HKPacket(
            filename=self.filename,
            file_offset=file_offset,
            file_captured_packet_index=file_captured_packet_index,
            file_hk_packet_index=file_hk_packet_index,
            pcap_time=pcap_time_ns,
            pcap_sec=pcap_sec,
            pcap_nsec=pcap_nsec,
            pcap_len=incl_len,
            pcap_src_ip=sender_ip,
            board_loc=board_loc,
            telescope_id=telescope_id,
            quabo_id=quabo_id,
            gti_index=0,
            gti_pcap_time=pcap_time_ns,
            data=payload
        )

    def _decode_science_packet(self, payload, pcap_sec, pcap_nsec, sender_ip, file_offset, file_captured_packet_index, file_science_packet_index, incl_len):
        # Science packets Total payload sizes:
        #   - 528 bytes (16-bit modes)
        #   - 272 bytes (8-bit modes)

        if len(payload) != 528 and len(payload) != 272:
            return None

        # PANOSETI metadata (16 bytes)
        meta_format = "<BBHHIIH"
        meta_size = struct.calcsize(meta_format)
        meta = struct.unpack(meta_format, payload[:meta_size])

        acq_mode = meta[0]

        if self._acq_mode_filter is not None and acq_mode not in self._acq_mode_filter:
            return None

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

        pcap_time_ns = int(pcap_sec) * 1000000000 + int(pcap_nsec)
        
        # Use stateful unroller to handle 1024s wraparound and calibration
        if board_loc not in self._time_unrollers:
            self._time_unrollers[board_loc] = QuaboTimeUnroller()
        
        event_time_ns = self._time_unrollers[board_loc].unroll(tai, nanosec, pcap_time_ns)

        return SciencePacket(
            filename=self.filename,
            file_offset=file_offset,
            file_captured_packet_index=file_captured_packet_index,
            file_science_packet_index=file_science_packet_index,
            pcap_time=pcap_time_ns,
            pcap_sec=pcap_sec,
            pcap_nsec=pcap_nsec,
            pcap_len=incl_len,
            pcap_src_ip=sender_ip,
            acq_mode=acq_mode,
            packet_ver=meta[1],
            packet_num=meta[2],
            board_loc=board_loc,
            telescope_id=telescope_id,
            quabo_id=quabo_id,
            tai=tai,
            nanosec=nanosec,
            event_time_sec=int(event_time_ns // 1000000000),
            event_time_nsec=int(event_time_ns % 1000000000),
            event_time=event_time_ns,
            event_time_good=True, # Calibrated offline
            flags=meta[6],
            gti_index=0,
            gti_pcap_time=pcap_time_ns,
            pix_data=pix_data
        )

    def close(self):
        self._file.close()

def get_panoseti_packets(filenames, gtis=None, verbose=False, science_port=60001, decode_hk=False, acq_mode=0x01):
    """
    Convenience function to get packets from one or more PANOSETI PCAP/PCAPNG files.
    'filenames' can be a single filename (string), a glob pattern (string),
    or a list/iterable of filenames.
    'gtis' can be a GTIFilter object, a dict with start/stop keys, or a list of such dicts.
    """
    gti_filter = gtis if isinstance(gtis, GTIFilter) else GTIFilter(gtis)
    if isinstance(filenames, str):
        if any(char in filenames for char in '*?['):
            files = sorted(glob.glob(filenames, recursive=True))
        else:
            files = [filenames]
    else:
        files = []
        for item in filenames:
            if isinstance(item, str) and any(char in item for char in '*?['):
                files.extend(sorted(glob.glob(item, recursive=True)))
            else:
                files.append(item)

    for filename in sorted(files):
        if(verbose):
            print(f"Processing file: {filename} (size: {os.path.getsize(filename)} bytes)")
        decoder = PcapDecoder(filename, science_port=science_port, decode_hk=decode_hk, acq_mode=acq_mode)
        try:
            for packet in decoder:
                # Use pcap_time for GTI filtering for both Science and HK packets
                is_good, gti_idx, gti_start_time = gti_filter.get_info(packet.pcap_time / 1e9)
                if is_good:
                    # Use pcap_time for gti_pcap_time since event_time is unreliable
                    # gti_pcap_time is now in nanoseconds
                    gti_pcap_time = packet.pcap_time
                    if gti_start_time != -float('inf'):
                        gti_pcap_time -= int(gti_start_time * 1e9)
                    yield packet._replace(gti_index=gti_idx, gti_pcap_time=gti_pcap_time)
        finally:
            decoder.close()
