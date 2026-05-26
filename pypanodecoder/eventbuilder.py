#!/usr/bin/env python3

# eventbuilder.py: Build full camera events from PANOSETI Science packets.

# Author: Stephen Fegan <sfegan@llr.in2p3.fr> (2026-05-11)
# Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

import struct
import os
import sys
import numpy as np

from .pcapdecoder import get_panoseti_packets, SciencePacket, GTIFilter

class CameraEvent:
    """
    Represents a full or partial camera event consisting of packets from up to 4 quabos.
    A camera event is defined as a group of packets from the same telescope (scope)
    that share the same hardware timestamp (TAI/nanosec).
    """
    def __init__(self, telescope_id, first_packet):
        self.telescope_id = telescope_id
        # Map of quabo_id (0-3) to SciencePacket
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
        self.gti_index = first_packet.gti_index
        self.gti_event_time = first_packet.gti_event_time

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
        return f"<CameraEvent tel={self.telescope_id} quabos={quabos} time={self.event_time:.9f} gti={self.gti_index} pcap_time={self.start_pcap_time:.6f}>"

class FilterInfo(dict):
    """
    Dict-like object that records applied filters and is callable to apply new filters.
    """
    def __init__(self, parent, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._parent = parent

    def __call__(self, min_quabos=None, min_delta_t=None):
        return self._parent._apply_filter(min_quabos=min_quabos, min_delta_t=min_delta_t)

class CameraImages:
    """
    Container for PANOSETI camera images and metadata.
    Allows access via attributes (e.g., .images) or keys (e.g., ['images']).
    """
    def __init__(self, images, event_times, gti_indexes, 
                 gti_event_times, quabo_masks, gtis=None, events=None, dtype=None, filter=None):
        self.images = np.array(images, dtype=dtype) if dtype else images
        self.event_times = event_times
        self.gti_indexes = gti_indexes
        self.gti_event_times = gti_event_times
        self.quabo_masks = quabo_masks
        self.gtis = gtis if gtis is not None else {}
        self.events = events if events is not None else []
        self.filter = FilterInfo(self, filter or {})

    @classmethod
    def concatenate(cls, objects):
        """Concatenates multiple CameraImages objects into one."""
        if not objects:
            raise ValueError("concatenate() requires a non-empty list of CameraImages objects")
        
        merged_gtis = {}
        for o in objects:
            merged_gtis.update(o.gtis)
            
        merged_filter = {}
        for o in objects:
            f = getattr(o, 'filter', {})
            for k, v in f.items():
                if v is not None:
                    if k in merged_filter:
                        merged_filter[k] = min(merged_filter[k], v)
                    else:
                        merged_filter[k] = v

        return cls(
            images=np.concatenate([o.images for o in objects], axis=2),
            event_times=np.concatenate([o.event_times for o in objects]),
            gti_indexes=np.concatenate([o.gti_indexes for o in objects]),
            gti_event_times=np.concatenate([o.gti_event_times for o in objects]),
            quabo_masks=np.concatenate([o.quabo_masks for o in objects]),
            gtis=merged_gtis,
            events=[ev for o in objects for ev in o.events],
            filter=merged_filter
        )

    @property
    def unique_gti_indexes(self):
        """Returns the sorted unique GTI indexes present in this object."""
        return np.unique(self.gti_indexes)

    def filter_gtis(self, gti_index):
        """Returns a new CameraImages object filtered for a specific GTI index."""
        mask = (self.gti_indexes == gti_index)
        
        filtered_events = []
        if self.events:
            filtered_events = [self.events[i] for i, m in enumerate(mask) if m]
            
        return CameraImages(
            images=self.images[:, :, mask],
            event_times=self.event_times[mask],
            gti_indexes=self.gti_indexes[mask],
            gti_event_times=self.gti_event_times[mask],
            quabo_masks=self.quabo_masks[mask],
            gtis={gti_index: self.gtis[gti_index]} if gti_index in self.gtis else {},
            events=filtered_events,
            filter=dict(self.filter)
        )

    def _apply_filter(self, min_quabos=None, min_delta_t=None):
        """
        Filter camera images based on specific cuts.
        
        Args:
            min_quabos (int, optional): Minimum number of Quabos required in the image.
            min_delta_t (float, optional): Minimum time (seconds) since the last event.
                                           Used to filter out high-frequency spikes.
                                           
        Returns:
            CameraImages: A new CameraImages object with filtered images and metadata.
        """
        keep = np.ones(len(self.event_times), dtype=bool)
        
        if min_quabos is not None:
            masks = np.asarray(self.quabo_masks, dtype=int)
            num_quabos = (
                (masks & 1) + 
                ((masks >> 1) & 1) + 
                ((masks >> 2) & 1) + 
                ((masks >> 3) & 1)
            )
            keep &= (num_quabos >= min_quabos)
            
        if min_delta_t is not None and len(self.event_times) > 0:
            delta_ts = np.concatenate(([np.inf], np.diff(self.event_times)))
            keep &= (delta_ts >= min_delta_t)
            
        filtered_events = []
        if self.events:
            filtered_events = [self.events[i] for i, m in enumerate(keep) if m]
            
        new_filter_dict = dict(self.filter)
        if min_quabos is not None:
            existing_mq = new_filter_dict.get('min_quabos')
            new_filter_dict['min_quabos'] = max(existing_mq, min_quabos) if existing_mq is not None else min_quabos
        if min_delta_t is not None:
            existing_mdt = new_filter_dict.get('min_delta_t')
            new_filter_dict['min_delta_t'] = max(existing_mdt, min_delta_t) if existing_mdt is not None else min_delta_t

        return CameraImages(
            images=self.images[:, :, keep],
            event_times=self.event_times[keep],
            gti_indexes=self.gti_indexes[keep],
            gti_event_times=self.gti_event_times[keep],
            quabo_masks=self.quabo_masks[keep],
            gtis=self.gtis,
            events=filtered_events,
            filter=new_filter_dict
        )

    def apply_pedestal_corrections(self, pcorr):
        """
        Applies a time-dependent polynomial correction to each pixel's image values.
        Correction for pixel (i, j) at event k is -polyval(pcorr[i, j, :], gti_event_times[k]).
        
        Args:
            pcorr (np.ndarray): 3D array of shape (32, 32, N_coeffs) containing 
                                polynomial coefficients (highest degree first).
                                
        Returns:
            CameraImages: A new object with adjusted images (dtype float).
        """
        num_coeffs = pcorr.shape[2]
        
        # Horner's method for efficient polynomial evaluation across the stack
        # res = p[0]*x^n + p[1]*x^(n-1) + ... + p[n]
        res = np.full(self.images.shape, pcorr[:, :, 0][..., np.newaxis], dtype=float)
        for k in range(1, num_coeffs):
            res = res * self.gti_event_times + pcorr[:, :, k][..., np.newaxis]
            
        return CameraImages(
            images=self.images - res,
            event_times=self.event_times,
            gti_indexes=self.gti_indexes,
            gti_event_times=self.gti_event_times,
            quabo_masks=self.quabo_masks,
            gtis=self.gtis,
            events=self.events,
            filter=dict(self.filter)
        )

    def __getitem__(self, key):
        return getattr(self, key)

    def __repr__(self):
        return f"<CameraImages events={len(self.event_times)} gtis={list(self.gtis.keys())}>"

class CameraEventBuilder:
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
        # active_events maps telescope_id -> { event_time -> CameraEvent }
        self.active_events = {}

    def process_packet(self, packet):
        """
        Process a single SciencePacket and yield any completed CameraEvents.
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
            self.active_events[telescope_id][event_time] = CameraEvent(telescope_id, packet)

    def flush(self):
        """
        Yield all remaining active events. Should be called at the end of the packet stream.
        """
        for tid in list(self.active_events.keys()):
            for et in list(self.active_events[tid].keys()):
                yield self.active_events[tid].pop(et)
        self.active_events.clear()

def get_camera_events(filenames, max_pcap_tdiff=1.0, max_event_tdiff=1e-6, gtis=None, verbose=False):
    """
    Convenience generator to get full camera events from one or more PANOSETI pcap files.
    
    Args:
        filenames (str or list): One or more paths to .pcap or .pcapng files, or a glob pattern.
        max_pcap_tdiff (float): Maximum PCAP arrival time difference.
        max_event_tdiff (float): Maximum hybridized event time difference (s).
        gtis (GTIFilter, dict or list): Good Time Intervals for filtering events.
        
    Yields:
        CameraEvent: Reconstructed camera events.
    """
    builder = CameraEventBuilder(max_pcap_tdiff=max_pcap_tdiff, max_event_tdiff=max_event_tdiff)
    for packet in get_panoseti_packets(filenames, gtis=gtis, verbose=verbose):
        for event in builder.process_packet(packet):
            yield event
    for event in builder.flush():
        yield event

def load_camera_images(filenames, gtis=None, min_quabos=None, store_camera_events=False, **kwargs):
    """
    Load all camera images from one or more PANOSETI pcap files.
    
    Args:
        filenames (str or list): One or more paths to .pcap or .pcapng files, or a glob pattern.
        gtis (GTIFilter, dict or list, optional): Good Time Intervals.
        min_quabos (int, optional): Minimum number of quabo packets required to form an image.
        store_camera_events (bool): If True, store the CameraEvent instances.

    Returns:
        CameraImages: Object containing the 32x32xNevent images and metadata.
    """
    images = []
    event_times = []
    gti_indexes = []
    gti_event_times = []
    quabo_masks = []
    events_list = []

    # Get GTI info
    gti_filter = gtis if isinstance(gtis, GTIFilter) else GTIFilter(gtis)
    gtis_dict = {}
    for start, stop, good, idx in gti_filter.intervals:
        if good:
            gtis_dict[idx] = {'start': start, 'stop': stop}

    for event in get_camera_events(filenames, gtis=gti_filter, **kwargs):
        if min_quabos is None or len(event.packets) >= min_quabos:
            images.append(event.get_image())
            event_times.append(event.event_time)
            gti_indexes.append(event.gti_index)
            gti_event_times.append(event.gti_event_time)
            
            # Calculate quabo bitmask
            mask = 0
            for qid in event.packets.keys():
                mask |= (1 << qid)
            quabo_masks.append(mask)
            
            if store_camera_events:
                events_list.append(event)

    if images:
        sort_idx = np.argsort(event_times)
        images = np.stack(images, axis=2)[:, :, sort_idx]
        event_times = np.array(event_times)[sort_idx]
        gti_indexes = np.array(gti_indexes)[sort_idx]
        gti_event_times = np.array(gti_event_times)[sort_idx]
        quabo_masks = np.array(quabo_masks)[sort_idx]
        if store_camera_events:
            events_list = [events_list[i] for i in sort_idx]
        else:
            events_list = None
    else:
        images = np.zeros((32, 32, 0))
        event_times = np.array([])
        gti_indexes = np.array([])
        gti_event_times = np.array([])
        quabo_masks = np.array([])
        events_list = None

    filter_dict = {}
    if min_quabos is not None:
        filter_dict['min_quabos'] = min_quabos

    return CameraImages(
        images=images,
        event_times=event_times,
        gti_indexes=gti_indexes,
        gti_event_times=gti_event_times,
        quabo_masks=quabo_masks,
        gtis=gtis_dict,
        events=events_list,
        filter=filter_dict
    )
