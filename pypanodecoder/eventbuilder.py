#!/usr/bin/env python3

# eventbuilder.py: Build full camera events from PANOSETI Science packets.

# Author: Stephen Fegan <sfegan@llr.in2p3.fr> (2026-05-11)
# Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

import struct
import os
import sys
import json
import glob
import numpy as np

from .pcapdecoder import get_panoseti_packets, SciencePacket, GTIFilter, wr_to_unix

class CameraEvent:
    """
    Represents a full or partial camera event consisting of packets from up to 4 quabos.
    A camera event is defined as a group of packets from the same telescope (scope)
    that share the same hardware timestamp (TAI/nanosec).
    """
    def __init__(self, telescope_id, first_packet, module_id=None):
        self.telescope_id = telescope_id
        self.module_id = module_id
        # Map of quabo_id (0-3) to SciencePacket
        self.packets = {first_packet.quabo_id: first_packet}
        
        # Timing from the first packet that started this event (nanoseconds)
        self.pcap_sec = first_packet.pcap_sec
        self.pcap_nsec = first_packet.pcap_nsec
        self.pcap_time = first_packet.pcap_time
        
        # Hybridized event time (nanoseconds since Unix epoch)
        self.event_time = first_packet.event_time
        
        # Metadata (taking from the first packet)
        self.acq_mode = first_packet.acq_mode
        self.tai = first_packet.tai
        self.nanosec = first_packet.nanosec
        self.packet_num = first_packet.packet_num
        self.board_loc = first_packet.board_loc
        self.gti_index = first_packet.gti_index
        self.gti_pcap_time = first_packet.gti_pcap_time

        # Per-quabo detailed timing arrays
        self.quabo_pcap_sec = np.full(4, np.nan)
        self.quabo_pcap_nsec = np.full(4, np.nan)
        self.quabo_pkt_sec = np.full(4, np.nan)
        self.quabo_pkt_nsec = np.full(4, np.nan)
        self.quabo_pcap_time = np.zeros(4, dtype=np.int64)
        self.quabo_event_time = np.zeros(4, dtype=np.int64)
        
        self._set_quabo_timing(first_packet)

    def _set_quabo_timing(self, packet):
        qid = packet.quabo_id
        self.quabo_pcap_sec[qid] = packet.pcap_sec
        self.quabo_pcap_nsec[qid] = packet.pcap_nsec
        self.quabo_pkt_sec[qid] = packet.tai
        self.quabo_pkt_nsec[qid] = packet.nanosec
        self.quabo_pcap_time[qid] = packet.pcap_time
        self.quabo_event_time[qid] = packet.event_time

    @staticmethod
    def extract_quabo_pixels(image, qid):
        """
        Extracts 16x16 quadrant from 32x32 image and reverses tiling.
        """
        if qid == 0:
            pix = np.flip(image[16:32, 16:32].T, axis=1)
        elif qid == 1:
            pix = image[16:32, 0:16]
        elif qid == 2:
            pix = np.flip(image[0:16, 0:16].T, axis=0)
        elif qid == 3:
            pix = np.flip(image[0:16, 16:32])
        else:
            return np.zeros(256)
        return pix.flatten()

    def add_packet(self, packet):
        """Adds a packet to the current camera event."""
        self.packets[packet.quabo_id] = packet
        self._set_quabo_timing(packet)

    @staticmethod
    def make_image(quabo0_pix=None, quabo1_pix=None, quabo2_pix=None, quabo3_pix=None,
                   transpose=False, dtype=None, missing_value=0):
        """
        Combines 4 quabo pixel arrays (256 elements each) into a 32x32 camera image.
        
        Args:
            quabo0_pix to quabo3_pix (array-like, optional): 256-pixel arrays.
            transpose (bool): Whether to transpose the resulting quadrants.
            dtype (numpy.dtype, optional): Data type of the output image. 
                                           If None, uses type of the first present quabo.
            missing_value: Value to use for quabos not present. Default 0.
            
        Returns:
            np.ndarray: A 32x32 array of pixel values.
        """
        # Determine dtype if not provided
        if dtype is None:
            for pix in [quabo0_pix, quabo1_pix, quabo2_pix, quabo3_pix]:
                if pix is not None:
                    dtype = np.asanyarray(pix).dtype
                    break
            if dtype is None:
                dtype = np.float64

        image = np.full((32, 32), missing_value, dtype=dtype)
        if transpose:
            if quabo0_pix is not None:
                pix = np.array(quabo0_pix).reshape((16, 16))
                image[16:32, 16:32] = np.flip(pix, axis=1)
            if quabo1_pix is not None:
                pix = np.array(quabo1_pix).reshape((16, 16))
                image[0:16, 16:32] = pix.T
            if quabo2_pix is not None:
                pix = np.array(quabo2_pix).reshape((16, 16))
                image[0:16, 0:16] = np.flip(pix, axis=0)
            if quabo3_pix is not None:
                pix = np.array(quabo3_pix).reshape((16, 16))
                image[16:32, 0:16] = np.flip(pix).T
        else:
            if quabo0_pix is not None:
                pix = np.array(quabo0_pix).reshape((16, 16))
                image[16:32, 16:32] = np.flip(pix, axis=1).T
            if quabo1_pix is not None:
                pix = np.array(quabo1_pix).reshape((16, 16))
                image[16:32, 0:16] = pix
            if quabo2_pix is not None:
                pix = np.array(quabo2_pix).reshape((16, 16))
                image[0:16, 0:16] = np.flip(pix, axis=0).T
            if quabo3_pix is not None:
                pix = np.array(quabo3_pix).reshape((16, 16))
                image[0:16, 16:32] = np.flip(pix)

        return image

    def get_image(self, transpose=False, dtype=None, missing_value=0):
        """
        Reconstructs the 32x32 camera image as a NumPy array.
        
        Returns:
            np.ndarray: A 32x32 array of pixel values.
        """
        return self.make_image(
            quabo0_pix=self.packets[0].pix_data if 0 in self.packets else None,
            quabo1_pix=self.packets[1].pix_data if 1 in self.packets else None,
            quabo2_pix=self.packets[2].pix_data if 2 in self.packets else None,
            quabo3_pix=self.packets[3].pix_data if 3 in self.packets else None,
            transpose=transpose,
            dtype=dtype,
            missing_value=missing_value
        )

    def __repr__(self):
        quabos = sorted(list(self.packets.keys()))
        return f"<CameraEvent tel={self.telescope_id} quabos={quabos} time={self.event_time}ns gti={self.gti_index} pcap_time={self.pcap_time}ns>"

class CameraImages:
    """
    Container for PANOSETI camera images and metadata.
    Allows access via attributes (e.g., .images) or keys (e.g., ['images']).
    """
    def __init__(self, images, event_times, gti_indexes, 
                 gti_pcap_times, quabo_masks, gtis=None, events=None, 
                 dtype=None, filter=None, source=None,
                 quabo_pcap_sec=None, quabo_pcap_nsec=None,
                 quabo_pkt_sec=None, quabo_pkt_nsec=None,
                 pcap_times=None, quabo_pcap_time=None, quabo_event_time=None):
        self.images = np.array(images, dtype=dtype) if dtype else images
        self.event_times = event_times
        self.pcap_times = pcap_times
        self.gti_indexes = gti_indexes
        self.gti_pcap_times = gti_pcap_times
        self.quabo_masks = quabo_masks
        self.gtis = gtis if gtis is not None else {}
        self.events = events if events is not None else []
        self.filter = dict(filter) if filter is not None else {}
        # source can be a string (applied to all GTIs) or a dict {gti_index: source_name}
        if isinstance(source, str):
            self.source = {idx: source for idx in self.unique_gti_indexes}
        else:
            self.source = dict(source) if source is not None else {}

        # Per-quabo timing (Nx4 arrays)
        self.quabo_pcap_sec = quabo_pcap_sec
        self.quabo_pcap_nsec = quabo_pcap_nsec
        self.quabo_pkt_sec = quabo_pkt_sec
        self.quabo_pkt_nsec = quabo_pkt_nsec
        self.quabo_pcap_time = quabo_pcap_time
        self.quabo_event_time = quabo_event_time

    @classmethod
    def concatenate(cls, objects):
        """Concatenates multiple CameraImages objects into one."""
        if not objects:
            raise ValueError("concatenate() requires a non-empty list of CameraImages objects")
        
        merged_gtis = {}
        merged_source = {}
        for o in objects:
            merged_gtis.update(o.gtis)
            merged_source.update(o.source)
            
        merged_filter = {}
        for o in objects:
            f = getattr(o, 'filter', {})
            for k, v in f.items():
                if v is not None:
                    if k in merged_filter:
                        merged_filter[k] = min(merged_filter[k], v)
                    else:
                        merged_filter[k] = v

        def _concat_timing(attr):
            parts = [getattr(o, attr) for o in objects if getattr(o, attr) is not None]
            return np.concatenate(parts) if parts else None

        return cls(
            images=np.concatenate([o.images for o in objects], axis=2),
            event_times=np.concatenate([o.event_times for o in objects]),
            pcap_times=_concat_timing('pcap_times'),
            gti_indexes=np.concatenate([o.gti_indexes for o in objects]),
            gti_pcap_times=np.concatenate([o.gti_pcap_times for o in objects]),
            quabo_masks=np.concatenate([o.quabo_masks for o in objects]),
            gtis=merged_gtis,
            events=[ev for o in objects for ev in o.events],
            filter=merged_filter,
            source=merged_source,
            quabo_pcap_sec=_concat_timing('quabo_pcap_sec'),
            quabo_pcap_nsec=_concat_timing('quabo_pcap_nsec'),
            quabo_pkt_sec=_concat_timing('quabo_pkt_sec'),
            quabo_pkt_nsec=_concat_timing('quabo_pkt_nsec'),
            quabo_pcap_time=_concat_timing('quabo_pcap_time'),
            quabo_event_time=_concat_timing('quabo_event_time')
        )

    @property
    def unique_gti_indexes(self):
        """Returns the sorted unique GTI indexes present in this object."""
        return np.unique(self.gti_indexes)

    @property
    def data_source(self):
        """Returns 'PCAP', 'PFF', or 'MIXED' based on the sources of all GTIs."""
        sources = set(self.source.values())
        if len(sources) > 1:
            return 'MIXED'
        elif len(sources) == 1:
            return list(sources)[0]
        else:
            return None

    def filter_gtis(self, gti_index):
        """Returns a new CameraImages object filtered for a specific GTI index."""
        mask = (self.gti_indexes == gti_index)
        
        filtered_events = []
        if self.events:
            filtered_events = [self.events[i] for i, m in enumerate(mask) if m]
            
        return CameraImages(
            images=self.images[:, :, mask],
            event_times=self.event_times[mask],
            pcap_times=self.pcap_times[mask] if self.pcap_times is not None else None,
            gti_indexes=self.gti_indexes[mask],
            gti_pcap_times=self.gti_pcap_times[mask],
            quabo_masks=self.quabo_masks[mask],
            gtis={gti_index: self.gtis[gti_index]} if gti_index in self.gtis else {},
            events=filtered_events,
            filter=dict(self.filter),
            source={gti_index: self.source.get(gti_index)},
            quabo_pcap_sec=self.quabo_pcap_sec[mask] if self.quabo_pcap_sec is not None else None,
            quabo_pcap_nsec=self.quabo_pcap_nsec[mask] if self.quabo_pcap_nsec is not None else None,
            quabo_pkt_sec=self.quabo_pkt_sec[mask] if self.quabo_pkt_sec is not None else None,
            quabo_pkt_nsec=self.quabo_pkt_nsec[mask] if self.quabo_pkt_nsec is not None else None,
            quabo_pcap_time=self.quabo_pcap_time[mask] if self.quabo_pcap_time is not None else None,
            quabo_event_time=self.quabo_event_time[mask] if self.quabo_event_time is not None else None
        )

    def filter_events(self, min_quabos=None, min_delta_t=None):
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
            
        if min_delta_t is not None and self.pcap_times is not None and len(self.pcap_times) > 0:
            delta_ts = np.concatenate(([np.inf], np.diff(self.pcap_times / 1e9)))
            keep &= (delta_ts >= min_delta_t)
        elif min_delta_t is not None and len(self.event_times) > 0:
            delta_ts = np.concatenate(([np.inf], np.diff(self.event_times / 1e9)))
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

        # We keep only sources for GTIs that still have events
        new_gti_indexes = np.unique(self.gti_indexes[keep])
        new_source = {idx: self.source[idx] for idx in new_gti_indexes if idx in self.source}

        return CameraImages(
            images=self.images[:, :, keep],
            event_times=self.event_times[keep],
            pcap_times=self.pcap_times[keep] if self.pcap_times is not None else None,
            gti_indexes=self.gti_indexes[keep],
            gti_pcap_times=self.gti_pcap_times[keep],
            quabo_masks=self.quabo_masks[keep],
            gtis=self.gtis,
            events=filtered_events,
            filter=new_filter_dict,
            source=new_source,
            quabo_pcap_sec=self.quabo_pcap_sec[keep] if self.quabo_pcap_sec is not None else None,
            quabo_pcap_nsec=self.quabo_pcap_nsec[keep] if self.quabo_pcap_nsec is not None else None,
            quabo_pkt_sec=self.quabo_pkt_sec[keep] if self.quabo_pkt_sec is not None else None,
            quabo_pkt_nsec=self.quabo_pkt_nsec[keep] if self.quabo_pkt_nsec is not None else None
        )

    def apply_pedestal_corrections(self, pcorr):
        """
        Applies a time-dependent polynomial correction to each pixel's image values.
        Correction for pixel (i, j) at event k is -polyval(pcorr[i, j, :], gti_pcap_times[k]).
        
        Args:
            pcorr (np.ndarray): 3D array of shape (32, 32, N_coeffs) containing 
                                polynomial coefficients (highest degree first).
                                
        Returns:
            CameraImages: A new object with adjusted images (dtype float).
        """
        num_coeffs = pcorr.shape[2]
        
        # Convert gti_pcap_times to seconds for numerical stability and compatibility with fit
        t_sec = self.gti_pcap_times / 1e9
        res = np.full(self.images.shape, pcorr[:, :, 0][..., np.newaxis], dtype=float)
        for k in range(1, num_coeffs):
            res = res * t_sec + pcorr[:, :, k][..., np.newaxis]

            
        return CameraImages(
            images=self.images - res,
            event_times=self.event_times,
            pcap_times=self.pcap_times,
            gti_indexes=self.gti_indexes,
            gti_pcap_times=self.gti_pcap_times,
            quabo_masks=self.quabo_masks,
            gtis=self.gtis,
            events=self.events,
            filter=dict(self.filter),
            source=dict(self.source),
            quabo_pcap_sec=self.quabo_pcap_sec,
            quabo_pcap_nsec=self.quabo_pcap_nsec,
            quabo_pkt_sec=self.quabo_pkt_sec,
            quabo_pkt_nsec=self.quabo_pkt_nsec
        )

    def map_gtis(self, functor, *args, **kwargs):
        """
        Applies a functor (callable) to each GTI-filtered CameraImages subset,
        and concatenates the results.
        
        Args:
            functor (callable): A function or callable object that takes a 
                                CameraImages instance as its first argument and 
                                returns a modified/new CameraImages instance.
            *args: Additional positional arguments passed to the functor.
            **kwargs: Additional keyword arguments passed to the functor.
            
        Returns:
            CameraImages: The concatenated result of applying the functor 
                           to each GTI.
        """
        results = []
        for gti_idx in self.unique_gti_indexes:
            gti_images = self.filter_gtis(gti_idx)
            res = functor(gti_images, *args, **kwargs)
            results.append(res)
        return CameraImages.concatenate(results)

    def __getitem__(self, key):
        return getattr(self, key)

    def __repr__(self):
        return f"<CameraImages events={len(self.event_times)} gtis={list(self.gtis.keys())}>"

class CameraEventBuilder:
    """
    Builds camera events from a stream of PANOSETI Science packets.
    Uses hybridized event_time for grouping and PCAP arrival time for flushing.
    """
    def __init__(self, max_pcap_tdiff=1.0, max_event_tdiff=1e-6, module_id=None, use_packet_num=False):
        """
        Args:
            max_pcap_tdiff (float): Maximum PCAP arrival time difference (seconds) 
                                     to wait for remaining quabos of an event.
                                     Default 1.0s handles significant network buffering.
            max_event_tdiff (float): Maximum hybridized event time difference (seconds)
                                    to consider quabos as part of the same event.
                                    Default is 1e-6s (1us).
            module_id (int, optional): Module ID to assign to camera events.
            use_packet_num (bool): If True, merge quabos into events based on packet_num.
        """
        self.max_pcap_tdiff_ns = int(max_pcap_tdiff * 1e9)
        self.max_event_tdiff_ns = int(max_event_tdiff * 1e9)
        self.module_id = module_id
        self.use_packet_num = use_packet_num
        # active_events maps telescope_id -> { key -> CameraEvent }
        # key is packet_num (if use_packet_num) or event_time
        self.active_events = {}

    def process_packet(self, packet):
        """
        Process a single SciencePacket and yield any completed CameraEvents.
        """
        pcap_time = packet.pcap_time
        event_time = packet.event_time
        telescope_id = packet.telescope_id
        quabo_id = packet.quabo_id
        packet_num = packet.packet_num
        
        # 1. Flush any events that have timed out (based on PCAP arrival time)
        if telescope_id in self.active_events:
            for key in list(self.active_events[telescope_id].keys()):
                event = self.active_events[telescope_id][key]
                if pcap_time - event.pcap_time > self.max_pcap_tdiff_ns:
                    yield self.active_events[telescope_id].pop(key)
        
        # 2. Match the current packet to an active event
        if telescope_id not in self.active_events:
            self.active_events[telescope_id] = {}
            
        # Search for an event with matching criteria
        matched_key = None
        if self.use_packet_num:
            if packet_num in self.active_events[telescope_id]:
                event = self.active_events[telescope_id][packet_num]
                if quabo_id not in event.packets:
                    matched_key = packet_num
                else:
                    # Duplicate quabo for this packet_num. Flush the old one.
                    yield self.active_events[telescope_id].pop(packet_num)
        else:
            # Search for an event with "close" event_time
            for et in self.active_events[telescope_id]:
                if abs(event_time - et) <= self.max_event_tdiff_ns:
                    event = self.active_events[telescope_id][et]
                    if quabo_id not in event.packets:
                        matched_key = et
                        break
                    else:
                        # Duplicate quabo for this hybridized time. Flush the old one.
                        yield self.active_events[telescope_id].pop(et)
                        break

        if matched_key is not None:
            self.active_events[telescope_id][matched_key].add_packet(packet)
            if len(self.active_events[telescope_id][matched_key].packets) == 4:
                yield self.active_events[telescope_id].pop(matched_key)
        else:
            # No matching event, start a new one
            key = packet_num if self.use_packet_num else event_time
            self.active_events[telescope_id][key] = CameraEvent(telescope_id, packet, module_id=self.module_id)

    def flush(self):
        """
        Yield all remaining active events. Should be called at the end of the packet stream.
        """
        for tid in list(self.active_events.keys()):
            # Dict keys iteration is insertion order in Python 3.7+
            for key in list(self.active_events[tid].keys()):
                yield self.active_events[tid].pop(key)
        self.active_events.clear()

def get_camera_events(filenames, max_pcap_tdiff=1.0, max_event_tdiff=1e-6, gtis=None, verbose=False, module_id=None, use_packet_num=False):
    """
    Convenience generator to get full camera events from one or more PANOSETI pcap files.
    
    Args:
        filenames (str or list): One or more paths to .pcap or .pcapng files, or a glob pattern.
        max_pcap_tdiff (float): Maximum PCAP arrival time difference.
        max_event_tdiff (float): Maximum hybridized event time difference (s).
        gtis (GTIFilter, dict or list): Good Time Intervals for filtering events.
        module_id (int, optional): Module ID to assign to camera events.
        use_packet_num (bool): If True, merge quabos into events based on packet_num.
        
    Yields:
        CameraEvent: Reconstructed camera events.
    """
    builder = CameraEventBuilder(max_pcap_tdiff=max_pcap_tdiff, max_event_tdiff=max_event_tdiff, module_id=module_id, use_packet_num=use_packet_num)
    for packet in get_panoseti_packets(filenames, gtis=gtis, verbose=verbose):
        for event in builder.process_packet(packet):
            yield event
    for event in builder.flush():
        yield event

def load_pcap_camera_images(filenames, gtis=None, min_quabos=None, store_camera_events=False, module_id=None, use_packet_num=False, no_sort=False, **kwargs):
    """
    Load all camera images from one or more PANOSETI pcap files.
    
    Args:
        filenames (str or list): One or more paths to .pcap or .pcapng files, or a glob pattern.
        gtis (GTIFilter, dict or list, optional): Good Time Intervals.
        min_quabos (int, optional): Minimum number of quabo packets required to form an image.
        store_camera_events (bool): If True, store the CameraEvent instances.
        module_id (int, optional): Module ID to assign to camera events.
        use_packet_num (bool): If True, merge quabos into events based on packet_num.
        no_sort (bool): If True, do not sort images by event_time.

    Returns:
        CameraImages: Object containing the 32x32xNevent images and metadata.
    """
    images = []
    event_times = []
    pcap_times = []
    gti_indexes = []
    gti_pcap_times = []
    quabo_masks = []
    events_list = []
    q_pcap_sec = []
    q_pcap_nsec = []
    q_pkt_sec = []
    q_pkt_nsec = []
    q_pcap_time = []
    q_event_time = []

    # Get GTI info
    gti_filter = gtis if isinstance(gtis, GTIFilter) else GTIFilter(gtis)
    gtis_dict = {}
    for start, stop, good, idx in gti_filter.intervals:
        if good:
            gtis_dict[idx] = {'start': start, 'stop': stop}

    for event in get_camera_events(filenames, gtis=gti_filter, module_id=module_id, use_packet_num=use_packet_num, **kwargs):
        if min_quabos is None or len(event.packets) >= min_quabos:
            images.append(event.get_image())
            event_times.append(event.event_time)
            pcap_times.append(event.pcap_time)
            gti_indexes.append(event.gti_index)
            gti_pcap_times.append(event.gti_pcap_time)

            # Calculate quabo bitmask
            mask = 0
            for qid in event.packets.keys():
                mask |= (1 << qid)
            quabo_masks.append(mask)

            # Detailed timing
            q_pcap_sec.append(event.quabo_pcap_sec)
            q_pcap_nsec.append(event.quabo_pcap_nsec)
            q_pkt_sec.append(event.quabo_pkt_sec)
            q_pkt_nsec.append(event.quabo_pkt_nsec)
            q_pcap_time.append(event.quabo_pcap_time)
            q_event_time.append(event.quabo_event_time)

            if store_camera_events:
                events_list.append(event)

    if images:
        if no_sort:
            sort_idx = np.arange(len(event_times))
        else:
            sort_idx = np.argsort(event_times)

        images = np.stack(images, axis=2)[:, :, sort_idx]
        event_times = np.array(event_times)[sort_idx]
        pcap_times = np.array(pcap_times)[sort_idx]
        gti_indexes = np.array(gti_indexes)[sort_idx]
        gti_pcap_times = np.array(gti_pcap_times)[sort_idx]
        quabo_masks = np.array(quabo_masks)[sort_idx]
        q_pcap_sec = np.array(q_pcap_sec)[sort_idx]
        q_pcap_nsec = np.array(q_pcap_nsec)[sort_idx]
        q_pkt_sec = np.array(q_pkt_sec)[sort_idx]
        q_pkt_nsec = np.array(q_pkt_nsec)[sort_idx]
        q_pcap_time = np.array(q_pcap_time)[sort_idx]
        q_event_time = np.array(q_event_time)[sort_idx]
        if store_camera_events:
            events_list = [events_list[i] for i in sort_idx]
        else:
            events_list = None
    else:
        images = np.zeros((32, 32, 0))
        event_times = np.array([])
        pcap_times = np.array([])
        gti_indexes = np.array([])
        gti_pcap_times = np.array([])
        quabo_masks = np.array([])
        q_pcap_sec = np.zeros((0, 4))
        q_pcap_nsec = np.zeros((0, 4))
        q_pkt_sec = np.zeros((0, 4))
        q_pkt_nsec = np.zeros((0, 4))
        q_pcap_time = np.zeros((0, 4), dtype=np.int64)
        q_event_time = np.zeros((0, 4), dtype=np.int64)
        events_list = None

    filter_dict = {}
    if min_quabos is not None:
        filter_dict['min_quabos'] = min_quabos

    return CameraImages(
        images=images,
        event_times=event_times,
        pcap_times=pcap_times,
        gti_indexes=gti_indexes,
        gti_pcap_times=gti_pcap_times,
        quabo_masks=quabo_masks,
        gtis=gtis_dict,
        events=events_list,
        filter=filter_dict,
        source='PCAP',
        quabo_pcap_sec=q_pcap_sec,
        quabo_pcap_nsec=q_pcap_nsec,
        quabo_pkt_sec=q_pkt_sec,
        quabo_pkt_nsec=q_pkt_nsec,
        quabo_pcap_time=q_pcap_time,
        quabo_event_time=q_event_time
    )
def load_pff_camera_images(filenames, gtis=None, min_quabos=None, store_camera_events=False, verbose=False, module_id=None, no_sort=False):
    """
    Load all camera images from one or more PANOSETI PFF files.
    
    Args:
        filenames (str or list): One or more paths to .pff files, or a glob pattern.
        gtis (GTIFilter, dict or list, optional): Good Time Intervals.
        min_quabos (int, optional): Minimum number of quabo packets required to form an image.
        store_camera_events (bool): If True, store the CameraEvent instances (mocked).
        verbose (bool): If True, print progress.
        module_id (int, optional): Module ID to assign to camera events.
        no_sort (bool): If True, do not sort images by event_time.

    Returns:
        CameraImages: Object containing the 32x32xNevent images and metadata.
    """
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
    files = sorted(files)

    images = []
    event_times = []
    pcap_times = []
    gti_indexes = []
    gti_pcap_times = []
    quabo_masks = []
    events_list = []
    q_pcap_sec = []
    q_pcap_nsec = []
    q_pkt_sec = []
    q_pkt_nsec = []
    q_pcap_time = []
    q_event_time = []

    gti_filter = gtis if isinstance(gtis, GTIFilter) else GTIFilter(gtis)
    gtis_dict = {}
    for start, stop, good, idx in gti_filter.intervals:
        if good:
            gtis_dict[idx] = {'start': start, 'stop': stop}

    RECORD_SIZE = 2540
    JSON_LEN = 491
    
    # State for unrolling TAI (board_loc -> QuaboTimeUnroller)
    time_unrollers = {}

    for filename in files:
        if verbose:
            print(f"Reading {filename}...")
        with open(filename, 'rb') as f:
            record_idx = 0
            while True:
                data = f.read(RECORD_SIZE)
                if len(data) < RECORD_SIZE:
                    break
                
                header_data = data[:JSON_LEN].decode('utf-8').strip()
                try:
                    header = json.loads(header_data)
                except json.JSONDecodeError:
                    record_idx += 1
                    continue
                
                present_quabos = []
                quabo_mask = 0
                pcap_time_ns = 0
                event_time_ns = 0
                # Initialize timing for this event (4 quabos)
                evt_pcap_sec = np.full(4, np.nan)
                evt_pcap_nsec = np.full(4, np.nan)
                evt_pkt_sec = np.full(4, np.nan)
                evt_pkt_nsec = np.full(4, np.nan)
                evt_pcap_time = np.zeros(4, dtype=np.int64)
                evt_event_time = np.zeros(4, dtype=np.int64)

                for i in range(4):
                    q_key = f"quabo_{i}"
                    q_info = header.get(q_key, {})
                    if q_info.get('tv_sec', 0) != 0:
                        present_quabos.append(i)
                        quabo_mask |= (1 << i)
                        
                        # board location for unroller
                        module_id_val = module_id if module_id is not None else 0
                        board_loc = (module_id_val << 2) | i

                        # capture time
                        this_q_pcap_sec = q_info['tv_sec']
                        this_q_pcap_usec = q_info['tv_usec']
                        this_q_pcap_time_ns = int(this_q_pcap_sec) * 1000000000 + int(this_q_pcap_usec) * 1000
                        
                        # hardware time
                        this_q_pkt_tai = q_info.get('pkt_tai', 0)
                        this_q_pkt_nsec = q_info.get('pkt_nsec', 0)
                        
                        # Determine event time using stateful unroller
                        if board_loc not in time_unrollers:
                            time_unrollers[board_loc] = QuaboTimeUnroller()
                        
                        this_q_event_time_ns = time_unrollers[board_loc].unroll(this_q_pkt_tai, this_q_pkt_nsec, this_q_pcap_time_ns)
                        
                        # Use the first available quabo to set event-wide times
                        if pcap_time_ns == 0:
                            pcap_time_ns = this_q_pcap_time_ns
                            event_time_ns = this_q_event_time_ns

                        evt_pcap_sec[i] = this_q_pcap_sec
                        evt_pcap_nsec[i] = this_q_pcap_usec * 1000
                        evt_pkt_sec[i] = this_q_pkt_tai
                        evt_pkt_nsec[i] = this_q_pkt_nsec
                        evt_pcap_time[i] = this_q_pcap_time_ns
                        evt_event_time[i] = this_q_event_time_ns
                
                if not present_quabos:
                    record_idx += 1
                    continue
                
                if min_quabos is not None and len(present_quabos) < min_quabos:
                    record_idx += 1
                    continue
                
                is_good, gti_idx, gti_start = gti_filter.get_info(pcap_time_ns / 1e9)
                if not is_good:
                    record_idx += 1
                    continue
                
                # Image data starts after JSON header and '*'
                image_raw = np.frombuffer(data[JSON_LEN+1:], dtype='<i2')
                image = image_raw.reshape((32, 32)).copy()
                
                # Rotate 180 degrees to match CameraImages standard
                image = np.flip(image)
                
                if store_camera_events:
                    class MockPacket:
                        def __init__(self, qid, q_info, module_id, gti_idx, gti_start, event_time_ns, pcap_time_ns, pix_data, packet_num):
                            self.quabo_id = qid
                            self.pcap_sec = q_info['tv_sec']
                            self.pcap_nsec = q_info['tv_usec'] * 1000
                            self.pcap_time = pcap_time_ns
                            self.event_time = event_time_ns
                            self.acq_mode = 0x01
                            self.tai = q_info.get('pkt_tai', 0)
                            self.nanosec = q_info.get('pkt_nsec', 0)
                            self.packet_num = packet_num
                            self.board_loc = ((module_id or 0) << 2) | qid
                            self.gti_index = gti_idx
                            self.gti_pcap_time = pcap_time_ns if gti_start == -float('inf') else pcap_time_ns - int(gti_start * 1e9)
                            self.pix_data = pix_data
                    event = None
                    for i in present_quabos:
                        q_info = header.get(f"quabo_{i}")
                        pix_data = CameraEvent.extract_quabo_pixels(image, i)
                        p = MockPacket(i, q_info, module_id, gti_idx, gti_start, event_time_ns, pcap_time_ns, pix_data, record_idx)
                        if event is None:
                            event = CameraEvent(module_id if module_id is not None else 0, p, module_id=module_id)
                        else:
                            event.add_packet(p)
                    # Force the event-wide pcap_time
                    event.pcap_time = pcap_time_ns
                    event.pcap_sec = int(pcap_time_ns // 1000000000)
                    event.pcap_nsec = int(pcap_time_ns % 1000000000)
                    events_list.append(event)

                images.append(image)
                event_times.append(event_time_ns)
                pcap_times.append(pcap_time_ns)
                gti_indexes.append(gti_idx)
                gti_pcap_times.append(pcap_time_ns if gti_start == -float('inf') else pcap_time_ns - int(gti_start * 1e9))
                quabo_masks.append(quabo_mask)
                q_pcap_sec.append(evt_pcap_sec)
                q_pcap_nsec.append(evt_pcap_nsec)
                q_pkt_sec.append(evt_pkt_sec)
                q_pkt_nsec.append(evt_pkt_nsec)
                q_pcap_time.append(evt_pcap_time)
                q_event_time.append(evt_event_time)
                record_idx += 1

    if images:
        if no_sort:
            sort_idx = np.arange(len(event_times))
        else:
            sort_idx = np.argsort(event_times)
            
        images = np.stack(images, axis=2)[:, :, sort_idx]
        event_times = np.array(event_times)[sort_idx]
        pcap_times = np.array(pcap_times)[sort_idx]
        gti_indexes = np.array(gti_indexes)[sort_idx]
        gti_pcap_times = np.array(gti_pcap_times)[sort_idx]
        quabo_masks = np.array(quabo_masks)[sort_idx]
        q_pcap_sec = np.array(q_pcap_sec)[sort_idx]
        q_pcap_nsec = np.array(q_pcap_nsec)[sort_idx]
        q_pkt_sec = np.array(q_pkt_sec)[sort_idx]
        q_pkt_nsec = np.array(q_pkt_nsec)[sort_idx]
        q_pcap_time = np.array(q_pcap_time)[sort_idx]
        q_event_time = np.array(q_event_time)[sort_idx]
        if store_camera_events:
            events_list = [events_list[i] for i in sort_idx]
        else:
            events_list = None
    else:
        images = np.zeros((32, 32, 0))
        event_times = np.array([])
        pcap_times = np.array([])
        gti_indexes = np.array([])
        gti_pcap_times = np.array([])
        quabo_masks = np.array([])
        q_pcap_sec = np.zeros((0, 4))
        q_pcap_nsec = np.zeros((0, 4))
        q_pkt_sec = np.zeros((0, 4))
        q_pkt_nsec = np.zeros((0, 4))
        q_pcap_time = np.zeros((0, 4), dtype=np.int64)
        q_event_time = np.zeros((0, 4), dtype=np.int64)
        events_list = None

    filter_dict = {}
    if min_quabos is not None:
        filter_dict['min_quabos'] = min_quabos

    return CameraImages(
        images=images,
        event_times=event_times,
        pcap_times=pcap_times,
        gti_indexes=gti_indexes,
        gti_pcap_times=gti_pcap_times,
        quabo_masks=quabo_masks,
        gtis=gtis_dict,
        events=events_list,
        filter=filter_dict,
        source='PFF',
        quabo_pcap_sec=q_pcap_sec,
        quabo_pcap_nsec=q_pcap_nsec,
        quabo_pkt_sec=q_pkt_sec,
        quabo_pkt_nsec=q_pkt_nsec,
        quabo_pcap_time=q_pcap_time,
        quabo_event_time=q_event_time
    )

def load_camera_images(filenames, gtis=None, priority="PCAP", **kwargs):
    """
    Unified loader for PANOSETI data (PCAP or PFF).
    Globs filenames, separates them by format, loads them, and merges them.
    If a GTI is present in both formats, 'priority' determines which is kept.

    Args:
        filenames (str or list): Glob pattern(s) or list of files.
        gtis (GTIFilter, dict or list, optional): Good Time Intervals.
        priority (str): "PCAP" or "PFF". Which format to prefer if GTIs overlap.
        **kwargs: Arguments passed to load_pcap_camera_images or load_pff_camera_images.

    Returns:
        CameraImages: The merged camera images.
    """
    if isinstance(filenames, str):
        if any(char in filenames for char in '*?['):
            all_files = sorted(glob.glob(filenames, recursive=True))
        else:
            all_files = [filenames]
    else:
        all_files = []
        for item in filenames:
            if isinstance(item, str) and any(char in item for char in '*?['):
                all_files.extend(sorted(glob.glob(item, recursive=True)))
            else:
                all_files.append(item)
    
    pcap_files = [f for f in all_files if f.lower().endswith(('.pcap', '.pcapng'))]
    pff_files = [f for f in all_files if f.lower().endswith('.pff')]

    images_pcap = None
    if pcap_files:
        images_pcap = load_pcap_camera_images(pcap_files, gtis=gtis, **kwargs)

    images_pff = None
    if pff_files:
        # load_pff_camera_images uses 'verbose' instead of '**kwargs' currently, 
        # but let's pass kwargs that it can handle (min_quabos, store_camera_events, verbose, module_id, no_sort)
        pff_kwargs = {k: v for k, v in kwargs.items() if k in ['min_quabos', 'store_camera_events', 'verbose', 'module_id', 'no_sort']}
        images_pff = load_pff_camera_images(pff_files, gtis=gtis, **pff_kwargs)

    if images_pcap and images_pff:
        # Check for overlapping GTIs
        pcap_gtis = set(images_pcap.unique_gti_indexes)
        pff_gtis = set(images_pff.unique_gti_indexes)
        overlap = pcap_gtis.intersection(pff_gtis)

        to_merge = []
        if priority.upper() == "PCAP":
            to_merge.append(images_pcap)
            for gidx in pff_gtis:
                if gidx not in overlap:
                    to_merge.append(images_pff.filter_gtis(gidx))
        else: # Priority PFF
            to_merge.append(images_pff)
            for gidx in pcap_gtis:
                if gidx not in overlap:
                    to_merge.append(images_pcap.filter_gtis(gidx))
        
        return CameraImages.concatenate(to_merge)
    
    return images_pcap or images_pff or CameraImages(np.zeros((32, 32, 0)), np.array([]), np.array([]), np.array([]), np.array([]), source='NONE', pcap_times=np.array([]))


