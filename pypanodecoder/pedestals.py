#!/usr/bin/env python3

# pedestals.py: Calculate charge spectra and pedestal values

# Author: Stephen Fegan <sfegan@llr.in2p3.fr> (2026-05-13)
# Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

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

class PanosetiChargeHistogram:
    """
    Represents a multi-dimensional histogram of charge values.
    """
    def __init__(self, qcenter, qhist, bin_width):
        self.qcenter = qcenter
        self.qhist = qhist
        self.bin_width = bin_width

    @property
    def shape(self):
        """Returns the shape of the histogram grid (e.g., (32, 32))."""
        return self.qhist.shape[:-1]

    def cdf(self):
        """
        Returns the normalized cumulative distribution function (CDF) for all elements.
        
        Returns:
            tuple: (edges, cdf)
                edges: np.ndarray of bin edges (length N+1)
                cdf: np.ndarray of normalized cumulative counts (shape: grid_shape + (N+1))
        """
        # Cumulative sum along the range axis
        cumsum = np.cumsum(self.qhist, axis=-1)
        total_counts = cumsum[..., -1]
        
        # If qcenter are centers, left edges are qcenter - 0.5 * bin_width
        left_edges = self.qcenter - 0.5 * self.bin_width
        edges = np.append(left_edges, left_edges[-1] + self.bin_width)
        
        # Result CDF array
        cdf_shape = self.qhist.shape[:-1] + (self.qhist.shape[-1] + 1,)
        cdf = np.zeros(cdf_shape)
        
        good_mask = total_counts > 0
        if good_mask.any():
            cdf[good_mask, 1:] = cumsum[good_mask] / total_counts[good_mask, np.newaxis]
            
        return edges, cdf

    def quantiles(self, q):
        """
        Computes the quantiles for all elements in the histogram grid.
        
        Args:
            q (float or iterable): The quantile(s) to compute (e.g., 0.5 for median).
            
        Returns:
            list of np.ndarray: A list containing an array of quantile values for each q.
        """
        # Ensure q is an array
        q_vals = np.atleast_1d(q)
        
        edges, cdf = self.cdf()
        
        # Flatten for interpolation
        flat_cdf = cdf.reshape(-1, cdf.shape[-1])
        grid_shape = self.shape
        num_elements = flat_cdf.shape[0]
        
        # Mask for histograms with data
        # Last element of CDF is 1.0 for those with data
        good_mask = flat_cdf[:, -1] > 0
        
        results = []
        for qv in q_vals:
            # Output array for this quantile
            res = np.full(num_elements, np.nan)
            
            # For each histogram, interpolate qv
            for i in range(num_elements):
                if good_mask[i]:
                    res[i] = np.interp(qv, flat_cdf[i], edges)
            
            results.append(res.reshape(grid_shape))
            
        return results[0] if np.isscalar(q) else results

    def mean(self):
        """Returns the weighted mean for all elements in the histogram grid."""
        total = np.sum(self.qhist, axis=-1)
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.sum(self.qhist * self.qcenter, axis=-1) / total

    def var(self):
        """Returns the weighted variance for all elements in the histogram grid."""
        total = np.sum(self.qhist, axis=-1)
        m = self.mean()
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.sum(self.qhist * (self.qcenter - m[..., np.newaxis])**2, axis=-1) / total

    def std(self):
        """Returns the weighted standard deviation for all elements in the histogram grid."""
        return np.sqrt(self.var())

    def median(self):
        """Returns the median for all elements in the histogram grid."""
        return self.quantiles(0.5)

    def iqr(self):
        """Returns the Interquartile Range (IQR) for all elements in the histogram grid."""
        q1, q3 = self.quantiles([0.25, 0.75])
        return q3 - q1

    def winsorized_mean(self, limits=(0.05, 0.05)):
        """
        Computes the Winsorized mean for all elements in the histogram grid.
        
        Args:
            limits (tuple): The fraction to cut from the (low, high) ends.
        """
        v_low, v_high = self.quantiles([limits[0], 1.0 - limits[1]])
        total = np.sum(self.qhist, axis=-1)
        
        # Clip bin centers to the Winsorized limits for each grid point
        centers_clipped = np.clip(self.qcenter, v_low[..., np.newaxis], v_high[..., np.newaxis])
        
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.sum(self.qhist * centers_clipped, axis=-1) / total

    def winsorized_var(self, limits=(0.05, 0.05)):
        """
        Computes the Winsorized variance for all elements in the histogram grid.
        
        Args:
            limits (tuple): The fraction to cut from the (low, high) ends.
        """
        v_low, v_high = self.quantiles([limits[0], 1.0 - limits[1]])
        total = np.sum(self.qhist, axis=-1)
        
        wm = self.winsorized_mean(limits=limits)
        centers_clipped = np.clip(self.qcenter, v_low[..., np.newaxis], v_high[..., np.newaxis])
        
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.sum(self.qhist * (centers_clipped - wm[..., np.newaxis])**2, axis=-1) / total

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
    # Centers are the expected integer values
    qcenter_pix = np.arange(qmin_pix, qmax_pix, w_pix)
    qhist_pix = np.zeros((32, 32, len(qcenter_pix)), dtype=int)

    # Histograms for (summed) charge in each SiPM 
    qcenter_sipm = np.arange(qmin_pix*64, qmax_pix*64, w_sipm)
    qhist_sipm = np.zeros((4, 4, len(qcenter_sipm)), dtype=int)

    # Histograms for (summed) charge in each Quabo
    qcenter_quabo = np.arange(qmin_pix*256, qmax_pix*256, w_quabo)
    qhist_quabo = np.zeros((2, 2, len(qcenter_quabo)), dtype=int)

    # Histograms for (summed) charge in the full camera
    qcenter_camera = np.arange(qmin_pix*1024, qmax_pix*1024, w_camera)
    qhist_camera = np.zeros(len(qcenter_camera), dtype=int)

    num_events = images.shape[2]

    ii_pix, jj_pix = np.meshgrid(range(32),range(32))
    ii_sipm, jj_sipm = np.meshgrid(range(4),range(4))
    ii_quabo, jj_quabo = np.meshgrid(range(2),range(2))

    for n in range(num_events):
        image = images[:,:,n]
        # Index: (val - center_0) / w + 0.5. Since center_0 is qmin, it's (val - qmin) / w + 0.5
        idx_pix = np.floor((image - qmin_pix) / w_pix + 0.5).astype(int)
        idx_pix = np.clip(idx_pix, 0, len(qcenter_pix) - 1)
        qhist_pix[jj_pix, ii_pix, idx_pix] += 1

        image_sum = np.zeros((33,33),dtype=float)
        image_sum[1:,1:] = np.cumsum(np.cumsum(image, axis=0), axis=1)

        image_sipm = np.diff(np.diff(image_sum[::8,::8],axis=0),axis=1)
        idx_sipm = np.floor((image_sipm - qmin_pix*64) / w_sipm + 0.5).astype(int)
        idx_sipm = np.clip(idx_sipm, 0, len(qcenter_sipm) - 1)
        qhist_sipm[jj_sipm, ii_sipm, idx_sipm] += 1

        image_quabo = np.diff(np.diff(image_sum[::16,::16],axis=0),axis=1)
        idx_quabo = np.floor((image_quabo - qmin_pix*256) / w_quabo + 0.5).astype(int)
        idx_quabo = np.clip(idx_quabo, 0, len(qcenter_quabo) - 1)
        qhist_quabo[jj_quabo, ii_quabo, idx_quabo] += 1

        image_camera = image_sum[32,32]
        idx_camera = int(np.floor((image_camera - qmin_pix*1024) / w_camera + 0.5))
        idx_camera = np.clip(idx_camera, 0, len(qcenter_camera) - 1)
        qhist_camera[idx_camera] += 1

    return PanosetiChargeSpectra(
        num_events=num_events,
        pix=PanosetiChargeHistogram(qcenter_pix, qhist_pix, w_pix),
        sipm=PanosetiChargeHistogram(qcenter_sipm, qhist_sipm, w_sipm),
        quabo=PanosetiChargeHistogram(qcenter_quabo, qhist_quabo, w_quabo),
        camera=PanosetiChargeHistogram(qcenter_camera, qhist_camera, w_camera)
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

def apply_polynomial_pedestal_correction(camera_images, norder=0, quantiles=(0.15, 0.85)):
    """
    Fits a polynomial pedestal model to each pixel and subtracts it.
    The fit is performed independently for each GTI.
    Events outside the specified quantiles are excluded from the fit for robustness.
    
    Args:
        camera_images (PanosetiCameraImages): The image container.
        norder (int): The order of the polynomial to fit (0=constant, 1=linear, etc.).
        quantiles (tuple): The (low, high) quantile limits for including data in the fit.
        
    Returns:
        PanosetiCameraImages: A new container with pedestal-subtracted images.
    """
    results = []
    for gti_idx in camera_images.unique_gti_indexes:
        gti_images = camera_images.filter_gti(gti_idx)
        
        # Calculate bounds for fitting (robust estimation of the quiet pedestal)
        charge_spectra = calculate_charge_spectra(gti_images, combine_gtis=True)
        q_low, q_high = charge_spectra.pix.quantiles(quantiles)
        
        t = gti_images.gti_event_times
        pcorr = np.zeros((32, 32, norder + 1))
        
        for i in range(32):
            for j in range(32):
                y = gti_images.images[i, j, :]
                mask = (y >= q_low[i, j]) & (y <= q_high[i, j])
                
                if np.sum(mask) > norder:
                    pcorr[i, j, :] = np.polyfit(t[mask], y[mask], norder)
                else:
                    # Fallback to simple mean if not enough data for polynomial fit
                    pcorr[i, j, -1] = np.mean(y)
                    
        # Use the class method to apply the subtraction
        results.append(gti_images.apply_pedestal_corrections(pcorr))
        
    return PanosetiCameraImages.concatenate(results)

