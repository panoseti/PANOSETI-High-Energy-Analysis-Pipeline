#!/usr/bin/env python3

# eventbuilder.py: Build full camera events from PANOSETI Science packets.

# Author: Steve Fegan (2026-05-11)

import struct
import os
import sys
import numpy as np

# Ensure we can import from the current directory
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.append(current_dir)

try:
    from panodecoder import get_panoseti_events, ScienceEvent
except ImportError:
    # Fallback if imported from outside the package
    from .panodecoder import get_panoseti_events, ScienceEvent

class PanosetiCameraEvent:
    """
    Represents a full or partial camera event consisting of packets from up to 4 quabos.
    A camera event is defined as a group of packets from the same telescope (scope)
    that share the same hardware timestamp (TAI/nanosec).
    """
    def __init__(self, telescope_id, first_packet):
        self.telescope_id = telescope_id
        # Map of quabo_id (0-3) to ScienceEvent
        self.packets = {first_packet.quabo_id: first_packet}
        
        # Timing from the first packet that started this event
        self.pcap_sec = first_packet.pcap_sec
        self.pcap_usec = first_packet.pcap_usec
        self.start_pcap_time = self.pcap_sec + self.pcap_usec / 1e6
        
        # Hybridized event time (seconds since Unix epoch)
        self.event_time = first_packet.event_time
        
        # Metadata (taking from the first packet)
        self.acq_mode = first_packet.acq_mode
        self.tai = first_packet.tai
        self.nanosec = first_packet.nanosec
        self.packet_num = first_packet.packet_num
        self.board_loc = first_packet.board_loc

    def add_packet(self, packet):
        """Adds a packet to the current camera event."""
        self.packets[packet.quabo_id] = packet

    def get_image(self):
        """
        Reconstructs the 32x32 camera image as a NumPy array.
        Applies rotations and positions according to the logic in makeevents.h.
        
        Returns:
            np.ndarray: A 32x32 array of pixel values.
        """
        image = np.zeros((32, 32), dtype=np.int32)
        
        for quabo_id, packet in self.packets.items():
            # Initial quabo image before rotation
            pix_data = np.array(packet.pix_data).reshape((16, 16))
            
            # In makeevents.h:
            # xbin = 16 - (i // 16), ybin = 1 + (i % 16)
            # Translating to 0-indexed [x, y] coordinates: x = 15 - row, y = col
            q_img = np.flip(pix_data, axis=0)
            
            # Apply rotations based on quabo index (quadrant)
            # 0: Top Right (180 deg CW)
            # 1: Top Left (90 deg CW)
            # 2: Bottom Left (No rotation)
            # 3: Bottom Right (270 deg CW)
            
            if quabo_id == 0: # Top Right
                q_img = np.flip(q_img) # 180 deg
            elif quabo_id == 1: # Top Left
                q_img = q_img.T
                q_img = np.flip(q_img, axis=1) # 90 deg CW
            elif quabo_id == 2: # Bottom Left
                pass # No rotation
            elif quabo_id == 3: # Bottom Right
                q_img = q_img.T
                q_img = np.flip(q_img, axis=0) # 270 deg CW
            
            # Determine offsets in 32x32 image
            xoff, yoff = 0, 0
            if quabo_id == 0: # Top Right
                xoff, yoff = 16, 16
            elif quabo_id == 1: # Top Left
                xoff, yoff = 0, 16
            elif quabo_id == 2: # Bottom Left
                xoff, yoff = 0, 0
            elif quabo_id == 3: # Bottom Right
                xoff, yoff = 16, 0
            
            image[xoff:xoff+16, yoff:yoff+16] = q_img
            
        return image

    def __repr__(self):
        quabos = sorted(list(self.packets.keys()))
        return f"<PanosetiCameraEvent tel={self.telescope_id} quabos={quabos} time={self.event_time:.9f} pcap_time={self.start_pcap_time:.6f}>"

class PanosetiEventBuilder:
    """
    Builds camera events from a stream of PANOSETI Science packets.
    Uses hybridized event_time for grouping and PCAP arrival time for flushing.
    """
    def __init__(self, max_pcap_tdiff=1.0, max_event_tdiff=1e-6):
        """
        Args:
            max_pcap_tdiff (float): Maximum PCAP arrival time difference (seconds) 
                                     to wait for remaining quabos of an event.
                                     Default 1.0s handles significant network buffering.
            max_event_tdiff (float): Maximum hybridized event time difference (seconds)
                                    to consider quabos as part of the same event.
                                    Default is 1e-6s (1us).
        """
        self.max_pcap_tdiff = max_pcap_tdiff
        self.max_event_tdiff = max_event_tdiff
        # active_events maps telescope_id -> { event_time -> PanosetiCameraEvent }
        self.active_events = {}

    def process_packet(self, packet):
        """
        Process a single ScienceEvent and yield any completed PanosetiCameraEvents.
        """
        pcap_time = packet.pcap_sec + packet.pcap_usec / 1e6
        event_time = packet.event_time
        telescope_id = packet.telescope_id
        quabo_id = packet.quabo_id
        
        # 1. Flush any events that have timed out (based on PCAP arrival time)
        for tid in list(self.active_events.keys()):
            for et in list(self.active_events[tid].keys()):
                event = self.active_events[tid][et]
                if pcap_time - event.start_pcap_time > self.max_pcap_tdiff:
                    yield self.active_events[tid].pop(et)
            if not self.active_events[tid]:
                del self.active_events[tid]
        
        # 2. Match the current packet to an active event
        if telescope_id not in self.active_events:
            self.active_events[telescope_id] = {}
            
        # Search for an event with "close" event_time
        matched_event_et = None
        for et in self.active_events[telescope_id]:
            if abs(event_time - et) <= self.max_event_tdiff:
                # Found a match. Check if we already have this quabo.
                if quabo_id not in self.active_events[telescope_id][et].packets:
                    matched_event_et = et
                    break
                else:
                    # Duplicate quabo for this hybridized time. 
                    # This implies the previous event is finished or this is a separate trigger.
                    # We flush the old one.
                    yield self.active_events[telescope_id].pop(et)
                    break

        if matched_event_et is not None:
            self.active_events[telescope_id][matched_event_et].add_packet(packet)
            if len(self.active_events[telescope_id][matched_event_et].packets) == 4:
                yield self.active_events[telescope_id].pop(matched_event_et)
        else:
            # No matching event, start a new one
            self.active_events[telescope_id][event_time] = PanosetiCameraEvent(telescope_id, packet)

    def flush(self):
        """
        Yield all remaining active events. Should be called at the end of the packet stream.
        """
        for tid in list(self.active_events.keys()):
            for et in list(self.active_events[tid].keys()):
                yield self.active_events[tid].pop(et)
        self.active_events.clear()

def get_camera_events(filenames, max_pcap_tdiff=1.0, max_event_tdiff=1e-6, gti=None):
    """
    Convenience generator to get full camera events from one or more PANOSETI pcap files.
    
    Args:
        filenames (str or list): One or more paths to .pcap or .pcapng files, or a glob pattern.
        max_pcap_tdiff (float): Maximum PCAP arrival time difference.
        max_event_tdiff (float): Maximum hybridized event time difference (s).
        gti (dict or list): Good Time Intervals for filtering events.
        
    Yields:
        PanosetiCameraEvent: Reconstructed camera events.
    """
    builder = PanosetiEventBuilder(max_pcap_tdiff=max_pcap_tdiff, max_event_tdiff=max_event_tdiff)
    for packet in get_panoseti_events(filenames, gti=gti):
        for event in builder.process_packet(packet):
            yield event
    for event in builder.flush():
        yield event

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <pcap_file_or_pattern> [...]")
        sys.exit(1)
    
    count = 0
    for event in get_camera_events(sys.argv[1:]):
        if count % 100 == 0:
            print(event)
        count += 1
    print(f"Total camera events found: {count}")
