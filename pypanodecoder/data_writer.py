#!/usr/bin/env python3

# data_writer.py: Write PANOSETI data in PCAPNG and PFF formats.

# Author: Stephen Fegan <sfegan@llr.in2p3.fr> (2026-06-14)
# Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

import os
import struct
import json
import datetime
import numpy as np
import platform

SITES = {
    'gattini': {
        'base_ip': '192.168.3.248',
        'daq_ip': '192.168.0.4',
        'scope': 'Gattini'
    },
    'winter': {
        'base_ip': '192.168.3.244',
        'daq_ip': '192.168.0.6',
        'scope': 'Winter'
    },
    'fern': {
        'base_ip': '192.168.3.240',
        'daq_ip': '192.168.0.9',
        'scope': 'Fern'
    },
    'pti': {
        'base_ip': '192.168.3.232',
        'daq_ip': '192.168.0.8',
        'scope': 'PTI'
    },
}

def get_site_info(module_id):
    """
    Returns site information for a given module ID.
    """
    if module_id is not None:
        for name, info in SITES.items():
            ip = info['base_ip'].split('.')
            mod = (int(ip[-2])*256 + int(ip[-1]))//4
            if mod == module_id:
                res = info.copy()
                res['module_id'] = mod
                res['name'] = name
                return res
    
    return {
        'scope': 'Unknown',
        'module_id': module_id if module_id is not None else 'unknown',
        'name': 'unknown',
        'daq_ip': '10.0.1.100',
        'base_ip': '10.0.1.0'
    }

def templated_filename(template, site_info, ts_utc, seq_id):
    """
    Builds a filename from a template string using site and time information.
    """
    dt = datetime.datetime.fromtimestamp(ts_utc, tz=datetime.timezone.utc)
    
    return template.format(
        scope=site_info.get('scope', 'Unknown'),
        name=site_info.get('name', 'unknown'),
        date=dt.strftime('%Y%m%d'),
        time=dt.strftime('%H%M%S'),
        isotime=dt.strftime('%Y-%m-%dT%H:%M:%SZ'),
        module=site_info.get('module_id', 'unknown'),
        seqid=seq_id
    )

class PcapWriter:
    """
    Writes PANOSETI Science packets or CameraEvents to a PCAPNG file.
    """
    def __init__(self, filename=None, module_id=None, scope=None, 
                 daq_ip=None, quabo_base_ip=None, 
                 filename_template='onsky_{scope}_{date}_{time}.pcapng',
                 app_name='pypanodecoder', comment=None):
        """
        Args:
            filename (str, optional): Target filename. If None, generated from template.
            module_id (int, optional): Module ID for metadata and filename.
            scope (str, optional): Override scope name for filename.
            daq_ip (str, optional): Override destination IP address for PCAP packets.
            quabo_base_ip (str, optional): Override base IP for quabos.
            filename_template (str): Template for automatic filename generation.
            app_name (str): Application name for PCAPNG header.
            comment (str, optional): Comment/Provenance for PCAPNG header.
        """
        self.filename = filename
        self.module_id = module_id
        self.app_name = app_name
        self.comment = comment
        
        # Get site info from module_id
        self.site_info = get_site_info(module_id)
        if scope: self.site_info['scope'] = scope
        if daq_ip: self.site_info['daq_ip'] = daq_ip
        if quabo_base_ip: self.site_info['base_ip'] = quabo_base_ip
        
        self.filename_template = filename_template
        self.file_handle = None
        self.seq_id = 0
        self.packet_count = 0

    def _open_file(self, timestamp):
        if self.file_handle is not None:
            return
            
        if self.filename is None:
            self.filename = templated_filename(self.filename_template, self.site_info, timestamp, self.seq_id)
        
        self.file_handle = open(self.filename, 'wb')
        self._write_shb()
        self._write_idb()

    def _write_shb(self):
        # Section Header Block (SHB)
        hw_str = platform.machine().encode('utf-8')
        hw_pad = (4 - (len(hw_str) % 4)) % 4
        os_str = f"Python {platform.python_version()} ({platform.system()})".encode('utf-8')
        os_pad = (4 - (len(os_str) % 4)) % 4
        appl_str = self.app_name.encode('utf-8')
        appl_pad = (4 - (len(appl_str) % 4)) % 4
        
        # Options: shb_hardware (2), shb_os (3), shb_userappl (4), opt_comment (1), opt_endofopt (0)
        options = (
            struct.pack('<HH', 2, len(hw_str)) + hw_str + b'\x00' * hw_pad +
            struct.pack('<HH', 3, len(os_str)) + os_str + b'\x00' * os_pad +
            struct.pack('<HH', 4, len(appl_str)) + appl_str + b'\x00' * appl_pad
        )
        
        if self.comment:
            comment_str = self.comment.encode('utf-8')
            comment_pad = (4 - (len(comment_str) % 4)) % 4
            options += struct.pack('<HH', 1, len(comment_str)) + comment_str + b'\x00' * comment_pad
            
        options += struct.pack('<HH', 0, 0)
        
        shb_body = struct.pack('<IHHq', 0x1A2B3C4D, 1, 0, -1)
        shb_len = 28 + len(options)
        shb_block = struct.pack('<II', 0x0A0D0D0A, shb_len) + shb_body + options + struct.pack('<I', shb_len)
        self.file_handle.write(shb_block)

    def _write_idb(self):
        # Interface Description Block (IDB)
        # LinkType 1 = Ethernet. Value 9 for if_tsresol means 10^-9 (nanoseconds).
        idb_body = struct.pack('<HHI', 1, 0, 65535)
        options = (
            struct.pack('<HHB', 9, 1, 9) + b'\x00' * 3 +
            struct.pack('<HH', 0, 0)
        )
        idb_len = 20 + len(options)
        idb_block = struct.pack('<II', 0x00000001, idb_len) + idb_body + options + struct.pack('<I', idb_len)
        self.file_handle.write(idb_block)

    def _get_quabo_ip(self, board_loc):
        # last 10 bits are board_loc
        base_parts = [int(x) for x in self.site_info['base_ip'].split('.')]
        base_int = (base_parts[0] << 24) | (base_parts[1] << 16) | (base_parts[2] << 8) | base_parts[3]
        ip_int = (base_int & ~0x3FF) | (board_loc & 0x3FF)
        return struct.pack('>I', ip_int)

    def _get_daq_ip_bytes(self):
        ip = self.site_info['daq_ip']
        parts = [int(x) for x in ip.split('.')]
        return struct.pack('>I', (parts[0] << 24) | (parts[1] << 16) | (parts[2] << 8) | parts[3])

    def write_science_packet(self, pcap_sec, pcap_nsec, board_loc, acq_mode, pkt_num, tai, nanosec, pix_data, flags=1):
        """
        Writes a single Science packet to the PCAPNG file.
        """
        if self.file_handle is None:
            self._open_file(pcap_sec)
            
        # 1. Ethernet Header (14 bytes)
        eth_hdr = b'\x00'*6 + b'\x00'*6 + b'\x08\x00'
        
        # 2. IP Header (20 bytes)
        src_ip = self._get_quabo_ip(board_loc)
        dest_ip = self._get_daq_ip_bytes()
        
        # Determine pixel bytes length
        if acq_mode == 0x01: # 16-bit signed
            pix_bytes = struct.pack(f'<{len(pix_data)}h', *pix_data)
        elif acq_mode == 0x02: # 16-bit unsigned
            pix_bytes = struct.pack(f'<{len(pix_data)}H', *pix_data)
        elif acq_mode == 0x06: # 8-bit unsigned
            pix_bytes = struct.pack(f'<{len(pix_data)}B', *pix_data)
        else:
            # Fallback
            if len(pix_data) == 256:
                pix_bytes = struct.pack(f'<256H', *pix_data)
            else:
                pix_bytes = bytes(pix_data)

        payload_len = 8 + 16 + len(pix_bytes)
        ip_total_len = 20 + payload_len
        # IP header: 0x45 (v4, IHL 5), TOS 0, Total Len, ID 0, Flags/Frag 0x4000 (DF), TTL 64, Proto 17 (UDP), Checksum 0, IPs
        ip_hdr = struct.pack('>BBHHHBBH4s4s', 
                             0x45, 0, ip_total_len, 0, 0x4000, 64, 17, 0, src_ip, dest_ip)
        
        # 3. UDP Header (8 bytes)
        udp_hdr = struct.pack('>HHHH', 60001, 60001, payload_len, 0)
        
        # 4. PANOSETI Science Header (16 bytes)
        sci_hdr = struct.pack('<BBHHIIH', acq_mode, 1, pkt_num % 65536, board_loc, int(tai), int(nanosec), flags)
        
        full_packet = eth_hdr + ip_hdr + udp_hdr + sci_hdr + pix_bytes
        
        # 5. PCAPNG Enhanced Packet Block (EPB)
        incl_len = len(full_packet)
        orig_len = incl_len
        ts_ns = int(pcap_sec) * 1_000_000_000 + int(pcap_nsec)
        ts_high = ts_ns >> 32
        ts_low = ts_ns & 0xFFFFFFFF
        
        padding = (4 - (incl_len % 4)) % 4
        block_tot_len = 28 + incl_len + padding + 4
        
        epb_header = struct.pack('<IIIIIII', 0x00000006, block_tot_len, 0, ts_high, ts_low, incl_len, orig_len)
        epb_footer = struct.pack('<I', block_tot_len)
        
        self.file_handle.write(epb_header + full_packet + b'\x00' * padding + epb_footer)
        self.packet_count += 1

    def write_science_packet_instance(self, pkt):
        """
        Writes a SciencePacket instance to the PCAPNG file.
        """
        self.write_science_packet(
            pcap_sec=pkt.pcap_sec,
            pcap_nsec=pkt.pcap_nsec,
            board_loc=pkt.board_loc,
            acq_mode=pkt.acq_mode,
            pkt_num=pkt.packet_num,
            tai=pkt.tai,
            nanosec=pkt.nanosec,
            pix_data=pkt.pix_data,
            flags=pkt.flags
        )

    def write_camera_event(self, event):
        """
        Writes all quabos of a CameraEvent to the PCAPNG file.
        Reconstructs packets from common fields.
        """
        module_id = getattr(event, 'module_id', None)
        if module_id is None:
            module_id = event.telescope_id
            
        for qid in range(4):
            if not np.isnan(event.quabo_pcap_sec[qid]):
                board_loc = (module_id << 2) | qid
                
                # Try to get pixel data for this quabo
                pix_data = None
                if qid in event.packets:
                    pix_data = getattr(event.packets[qid], 'pix_data', None)
                
                if pix_data is None:
                    # Reconstruct from the 32x32 image
                    image = event.get_image()
                    from .eventbuilder import CameraEvent
                    pix_data = CameraEvent.extract_quabo_pixels(image, qid)
                
                self.write_science_packet(
                    pcap_sec=event.quabo_pcap_sec[qid],
                    pcap_nsec=event.quabo_pcap_nsec[qid],
                    board_loc=board_loc,
                    acq_mode=event.acq_mode,
                    pkt_num=event.packet_num,
                    tai=event.quabo_pkt_sec[qid],
                    nanosec=event.quabo_pkt_nsec[qid],
                    pix_data=pix_data,
                    flags=1 # Reconstructed events are usually "good"
                )

    def close(self):
        if self.file_handle:
            self.file_handle.close()
            self.file_handle = None

class PffWriter:
    """
    Writes CameraEvents to a PFF file.
    """
    def __init__(self, filename=None, module_id=None, scope=None,
                 filename_template='start_{isotime}.dp_ph1024.bpp_2.module_{module}.seqno_{seqid}.pff'):
        """
        Args:
            filename (str, optional): Target filename. If None, generated from template.
            module_id (int, optional): Module ID for metadata and filename.
            scope (str, optional): Override scope name for filename.
            filename_template (str): Template for automatic filename generation.
        """
        self.filename = filename
        self.module_id = module_id
        
        # Get site info from module_id
        self.site_info = get_site_info(module_id)
        if scope: self.site_info['scope'] = scope
        
        self.filename_template = filename_template
        self.file_handle = None
        self.seq_id = 0

    def _open_file(self, timestamp):
        if self.file_handle is not None:
            return
            
        if self.filename is None:
            self.filename = templated_filename(self.filename_template, self.site_info, timestamp, self.seq_id)
        
        self.file_handle = open(self.filename, 'wb')

    def write_camera_event(self, event):
        """
        Writes a CameraEvent to the PFF file.
        """
        if self.file_handle is None:
            self._open_file(event.event_time)
            
        # PFF record size: 491 (JSON) + 1 (*) + 2048 (Image) = 2540
        header = {}
        for i in range(4):
            if not np.isnan(event.quabo_pcap_sec[i]):
                header[f"quabo_{i}"] = {
                    'tv_sec': int(event.quabo_pcap_sec[i]),
                    'tv_usec': int(event.quabo_pcap_nsec[i] / 1000),
                    'pkt_tai': int(event.quabo_pkt_sec[i]),
                    'pkt_nsec': int(event.quabo_pkt_nsec[i])
                }
        
        header_json = json.dumps(header)
        header_bytes = header_json.encode('utf-8').ljust(491, b' ')
        
        # Image data (32x32)
        image = event.get_image()
        # Rotate back to PFF orientation
        image_to_write = np.flip(image)
        image_bytes = image_to_write.astype('<i2').tobytes()
        
        self.file_handle.write(header_bytes + b'*' + image_bytes)

    def close(self):
        if self.file_handle:
            self.file_handle.close()
            self.file_handle = None
