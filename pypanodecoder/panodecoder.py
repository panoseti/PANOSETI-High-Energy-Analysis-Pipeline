#!/usr/bin/env python3

# panodecoder.py: A pure Python decoder for PANOSETI data stored in PCAP format.
# Does not rely on any external packages.

# Author: Steve Fegan (2026-05-10)

import struct
from collections import namedtuple
import os
import glob

# Define the structure for a PANOSETI Science Event
# Using namedtuple for a "struct-like" behavior as requested.
ScienceEvent = namedtuple('ScienceEvent', [
    'filename', 'file_offset', 'packet_idx',
    'pcap_sec', 'pcap_usec', 'incl_len',
    'acq_mode', 'packer_ver', 'packet_no', 'boardloc',
    'telescope_id', 'quabo_id',
    'tai', 'nanosec', 'event_time', 'event_time_sec', 'event_time_nsec', 'event_time_good',
    'dummy', 'pix_data'
])

def wr_to_unix(pkt_tai, pkt_nsec, tv_sec, ignore_clock_desync=False):
    """
    Given a WR packet time (TAI) with only 10 bits of sec,
    and a Unix time that's within a few ms,
    return the complete WR time (in Unix time, not TAI).
    """
    # 37 is the TAI-UTC offset. 
    # The packet TAI seconds is only 10 bits (0-1023).
    d = (tv_sec - pkt_tai + 37) % 1024
    if d == 0:
        return True, tv_sec, pkt_nsec, tv_sec + pkt_nsec / 1e9
    elif d == 1:
        return True, tv_sec - 1, pkt_nsec, tv_sec - 1 + pkt_nsec / 1e9
    elif d == 1023:
        return True, tv_sec + 1, pkt_nsec, tv_sec + 1 + pkt_nsec / 1e9
    else:
        # The WR and DAQ clocks differ by > 1s => out of sync
        # Return 0 if ignore_clock_desync is False. Otherwise, return an approximation to the time.
        if ignore_clock_desync:
            approx_t = tv_sec + pkt_nsec / 1e9
            return False, tv_sec, pkt_nsec, approx_t
        else:
            raise Exception('WR and Unix times differ by > 1 sec: pkt_tai %d tv_sec %d d %d' % (pkt_tai, tv_sec, d))

class PanosetiPcapDecoder:
    """
    A pure Python decoder for PANOSETI data stored in PCAP format.
    Does not rely on any external packages.

    See: Jamie's C++ code
    See: https://github.com/panoseti/panoseti/wiki/Quabo-packet-interface#science-packets
    """
    def __init__(self, filename):
        self.filename = filename
        self._file = open(filename, 'rb')
        self._is_pcapng = False
        self._endian = '<'
        self._ts_resol = 1000000 # Default to microseconds (10^-6)
        self.packet_idx = 0
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
            
            ts_sec, ts_usec, incl_len, orig_len = struct.unpack(f"{self._endian}IIII", header_data)
            packet_data = self._file.read(incl_len)
            
            event = self._parse_packet(packet_data, ts_sec, ts_usec, offset, self.packet_idx)
            if event:
                yield event
                self.packet_idx += 1

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
                ts_sec = timestamp // self._ts_resol
                # Calculate usec based on resolution
                if self._ts_resol == 1000000:
                    ts_usec = timestamp % 1000000
                elif self._ts_resol > 1000000:
                    ts_usec = (timestamp % self._ts_resol) // (self._ts_resol // 1000000)
                else:
                    ts_usec = (timestamp % self._ts_resol) * (1000000 // self._ts_resol)
                
                packet_data = block_data[20:20+incl_len]
                event = self._parse_packet(packet_data, ts_sec, ts_usec, offset, self.packet_idx)
                if event:
                    yield event
                    self.packet_idx += 1
            # Skip other block types (SHB, IDB, etc. are already read into block_data)

    def _parse_packet(self, packet_data, ts_sec, ts_usec, file_offset, packet_idx):
        # Based on Quabo Packet Interface Wiki:
        # Science packets are sent to port 60001.
        # Total payload sizes: 
        #   - 528 bytes (16-bit modes)
        #   - 272 bytes (8-bit modes)
        # 42 bytes of Eth/IP/UDP headers are expected if captured on the wire.
        
        # We handle cases where we might have just the payload or the full capture.
        if len(packet_data) == 570: # 42 + 528
            payload = packet_data[42:]
        elif len(packet_data) == 528:
            payload = packet_data
        elif len(packet_data) == 314: # 42 + 272
            payload = packet_data[42:]
        elif len(packet_data) == 272:
            payload = packet_data
        else:
            # Not a recognized science packet size
            return None
        
        # PANOSETI metadata (16 bytes)
        # acq_mode(1), packer_ver(1), packet_no(2), boardloc(2), TAI(4), nanosec(4), dummy(2)
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
        
        boardloc = meta[3]
        telescope_id = boardloc >> 2
        quabo_id = boardloc & 0x3
        
        tai = meta[4]
        nanosec = meta[5]
        event_time_good, event_time_sec, event_time_nsec, event_time = wr_to_unix(tai, nanosec, ts_sec, ignore_clock_desync=True)

        return ScienceEvent(
            filename=self.filename,
            file_offset=file_offset,
            packet_idx=packet_idx,
            pcap_sec=ts_sec,
            pcap_usec=ts_usec,
            incl_len=len(packet_data),
            acq_mode=acq_mode,
            packer_ver=meta[1],
            packet_no=meta[2],
            boardloc=boardloc,
            telescope_id=telescope_id,
            quabo_id=quabo_id,
            tai=tai,
            nanosec=nanosec,
            event_time_sec=event_time_sec,
            event_time_nsec=event_time_nsec,
            event_time=event_time,
            event_time_good=event_time_good,
            dummy=meta[6],
            pix_data=pix_data
        )

    def close(self):
        self._file.close()

def get_panoseti_events(filenames):
    """
    Convenience function to get events from one or more PANOSETI PCAP/PCAPNG files.
    'filenames' can be a single filename (string), a glob pattern (string), 
    or a list/iterable of filenames.
    """
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

    for filename in files:
        decoder = PanosetiPcapDecoder(filename)
        try:
            for event in decoder:
                yield event
        finally:
            decoder.close()

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <pcap_file_or_pattern> [...]")
        sys.exit(1)
    
    count = 0
    # Process all arguments passed on the command line
    for event in get_panoseti_events(sys.argv[1:]):
        if count % 1000 == 0:
            # Use _asdict() for pretty printing if needed, or just repr
            print(f"Event {count}: {event._replace(pix_data='...')}")
        count += 1
    print(f"Total events: {count}")
