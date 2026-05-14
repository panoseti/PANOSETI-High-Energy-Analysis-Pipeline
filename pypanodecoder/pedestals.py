#!/usr/bin/env python3

# pedestals.py: Calculate pedestal values for PANOSETI Science packets.

# Author: Steve Fegan (2026-05-13)

import os
import sys
import numpy as np

# Ensure we can import from the current directory
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.append(current_dir)

try:
    from eventbuilder import get_camera_events, PanosetiCameraEvent, PanosetiCameraImages
except ImportError:
    # Fallback if imported from outside the package
    from .eventbuilder import get_camera_events, PanosetiCameraEvent, PanosetiCameraImages

class PanosetiChargeSpectra:
    """
    Container for PANOSETI charge spectra (pedestal histograms).
    Allows access via attributes (e.g., .pix) or keys (e.g., ['pix']).
    """
    def __init__(self, num_events, pix, sipm, quabo, camera):
        self.num_events = num_events
        self.pix = pix
        self.sipm = sipm
        self.quabo = quabo
        self.camera = camera

    def __getitem__(self, key):
        return getattr(self, key)

    def __repr__(self):
        return f"<PanosetiChargeSpectra events={self.num_events}>"

def _build_spectra(images, qmin_pix=-512, qmax_pix=4096, downsample=False):
    """Internal helper to build spectra from a 3D numpy array of images (32, 32, N)."""
    
    # Bin widths
    w_pix = 1
    w_sipm = 8 if downsample else 1
    w_quabo = 16 if downsample else 1
    w_camera = 32 if downsample else 1

    # Histogram for charge each pixel (32x32)
    qrange_pix = np.arange(qmin_pix, qmax_pix, w_pix)
    qhist_pix = np.zeros((32, 32, len(qrange_pix)), dtype=int)

    # Histograms for (summed) charge in each SiPM 
    qrange_sipm = np.arange(qmin_pix*64, qmax_pix*64, w_sipm)
    qhist_sipm = np.zeros((4, 4, len(qrange_sipm)), dtype=int)

    # Histograms for (summed) charge in each Quabo
    qrange_quabo = np.arange(qmin_pix*256, qmax_pix*256, w_quabo)
    qhist_quabo = np.zeros((2, 2, len(qrange_quabo)), dtype=int)

    # Histograms for (summed) charge in the full camera
    qrange_camera = np.arange(qmin_pix*1024, qmax_pix*1024, w_camera)
    qhist_camera = np.zeros(len(qrange_camera), dtype=int)

    num_events = images.shape[2]

    ii_pix, jj_pix = np.meshgrid(range(32),range(32))
    ii_sipm, jj_sipm = np.meshgrid(range(4),range(4))
    ii_quabo, jj_quabo = np.meshgrid(range(2),range(2))

    for n in range(num_events):
        image = np.maximum(np.minimum(images[:,:,n], qmax_pix-1), qmin_pix)
        qhist_pix[jj_pix, ii_pix, (image - qmin_pix) // w_pix] += 1

        image_sum = np.zeros((33,33),dtype=int)
        image_sum[1:,1:] = np.cumsum(np.cumsum(image, axis=0), axis=1)

        image_sipm = np.diff(np.diff(image_sum[::8,::8],axis=0),axis=1)
        # Ensure image_sipm is within range for binned indexing
        image_sipm = np.maximum(np.minimum(image_sipm, qmax_pix*64 - 1), qmin_pix*64)
        qhist_sipm[jj_sipm, ii_sipm, (image_sipm - qmin_pix*64) // w_sipm] += 1

        image_quabo = np.diff(np.diff(image_sum[::16,::16],axis=0),axis=1)
        image_quabo = np.maximum(np.minimum(image_quabo, qmax_pix*256 - 1), qmin_pix*256)
        qhist_quabo[jj_quabo, ii_quabo, (image_quabo - qmin_pix*256) // w_quabo] += 1

        image_camera = image_sum[32,32]
        image_camera = np.maximum(np.minimum(image_camera, qmax_pix*1024 - 1), qmin_pix*1024)
        qhist_camera[(image_camera - qmin_pix*1024) // w_camera] += 1

    return PanosetiChargeSpectra(
        num_events=num_events,
        pix=(qrange_pix, qhist_pix),
        sipm=(qrange_sipm, qhist_sipm),
        quabo=(qrange_quabo, qhist_quabo),
        camera=(qrange_camera, qhist_camera)
    )

def calculate_charge_spectra(camera_images, gti_indexes=None, combine_gtis=True, 
                             qmin_pix=-512, qmax_pix=4096, downsample=False):
    """
    Calculate the charge spectrum for each pixel from a PanosetiCameraImages object.
    
    Args:
        camera_images (PanosetiCameraImages): The loaded camera images and metadata.
        gti_indexes (list, optional): List of GTI indexes to include. If None, all are included.
        combine_gtis (bool): If True, combine all selected GTIs into one spectra object.
                            If False, return a dict of spectra objects keyed by GTI index.
        qmin_pix (int): Minimum charge for pixel histograms.
        qmax_pix (int): Maximum charge for pixel histograms.
        downsample (bool): If True, use larger bins for SiPM, Quabo, and camera histograms.

    Returns:
        PanosetiChargeSpectra or dict: The calculated spectra.
    """
    if gti_indexes is None:
        mask = np.ones(len(camera_images.gti_indexes), dtype=bool)
    else:
        mask = np.isin(camera_images.gti_indexes, gti_indexes)

    if combine_gtis:
        return _build_spectra(camera_images.images[:, :, mask], 
                              qmin_pix=qmin_pix, qmax_pix=qmax_pix, downsample=downsample)
    else:
        results = {}
        unique_gtis = np.unique(camera_images.gti_indexes[mask])
        for gti_idx in unique_gtis:
            gti_mask = (camera_images.gti_indexes == gti_idx)
            results[gti_idx] = _build_spectra(camera_images.images[:, :, gti_mask],
                                              qmin_pix=qmin_pix, qmax_pix=qmax_pix, downsample=downsample)
        return results

