#!/usr/bin/env python3

# pedestals.py: Calculate charge spectra and pedestal values

# Author: Stephen Fegan <sfegan@llr.in2p3.fr> (2026-05-13)
# Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

import os
import sys
import numpy as np

from .eventbuilder import CameraImages

class ChargeHistogram:
    """
    Represents a multi-dimensional histogram of charge values.
    """
    def __init__(self, qcenter, qhist, bin_width):
        self.qcenter = qcenter
        self.qhist = qhist
        self.bin_width = bin_width

    @classmethod
    def gaussian(cls, sigma, qcenter, mean=None):
        """
        Generates a Gaussian ChargeHistogram.
        
        Args:
            sigma (array-like): Standard deviation(s) of the Gaussian(s).
            qcenter (array-like): Bin centers for the histogram.
            mean (array-like, optional): Mean(s) of the Gaussian(s). 
                                       Defaults to None (all zeros).
            
        Returns:
            ChargeHistogram: A new histogram containing Gaussian densities.
        """
        from scipy.stats import norm
        sigma = np.asanyarray(sigma)
        if mean is None:
            mean = np.zeros_like(sigma)
        else:
            mean = np.asanyarray(mean)

        if len(qcenter) < 2:
            raise ValueError("qcenter must have at least two values to determine bin_width")
        
        bin_width = qcenter[1] - qcenter[0]
        edges = np.append(qcenter - 0.5 * bin_width, qcenter[-1] + 0.5 * bin_width)
        
        # Broadcast mean and sigma to (grid..., 1) and edges to (N+1,)
        # scipy.stats.norm.cdf handles the broadcasting
        cdf_edges = norm.cdf(edges, loc=mean[..., np.newaxis], scale=sigma[..., np.newaxis])
        
        # Probabilities in each bin
        probabilities = np.diff(cdf_edges, axis=-1)
        
        # Normalize as density: divide by bin_width
        qhist = probabilities / bin_width
        
        return cls(qcenter, qhist, bin_width)

    @property
    def shape(self):
        """Returns the shape of the histogram grid (e.g., (32, 32))."""
        return self.qhist.shape[:-1]

    def extract(self, *index):
        """
        Extracts a single element (or sub-grid) from the histogram grid.
        
        Args:
            index: The index or slice to extract.
            
        Returns:
            ChargeHistogram: A new histogram containing only the selected elements.
        """
        # qhist is (grid_shape..., num_bins)
        # We slice the grid dimensions; the last dimension (bins) is preserved.
        new_qhist = self.qhist[index]
        return ChargeHistogram(self.qcenter, new_qhist, self.bin_width)

    def __add__(self, other):
        """
        Sums two ChargeHistograms together.
        The histograms must have the same grid shape, bin centers, and bin width.
        """
        if not isinstance(other, ChargeHistogram):
            return NotImplemented
        
        if self.shape != other.shape:
            raise ValueError(f"Incompatible shapes: {self.shape} vs {other.shape}")
        
        if self.bin_width != other.bin_width:
            raise ValueError(f"Incompatible bin widths: {self.bin_width} vs {other.bin_width}")

        if not np.array_equal(self.qcenter, other.qcenter):
            # Fallback to check if they are "close enough" if needed, 
            # but usually they should be identical.
            raise ValueError("Bin centers do not match")
            
        return ChargeHistogram(self.qcenter, self.qhist + other.qhist, self.bin_width)

    def __mul__(self, other):
        """
        Performs a numerical convolution of two ChargeHistograms across the grid.
        This is mapped to the multiplication operator (*).
        
        The result represents the distribution of the sum of two independent 
        random variables, each distributed according to one of the histograms.
        
        Note: This assumes the binning is uniform and identical. 
        The resulting histogram will be centered around the sum of the means,
        but since we force it back into the same qcenter range, 
        values falling outside the range will be lost (truncated).
        """
        if not isinstance(other, ChargeHistogram):
            return NotImplemented
        
        if self.shape != other.shape:
            raise ValueError(f"Incompatible shapes: {self.shape} vs {other.shape}")

        if self.bin_width != other.bin_width:
            raise ValueError(f"Incompatible bin widths: {self.bin_width} vs {other.bin_width}")

        if not np.array_equal(self.qcenter, other.qcenter):
            raise ValueError("Bin centers do not match")

        # Numerical convolution along the last axis (the histogram bins)
        from scipy.signal import fftconvolve
        
        # qhist is (grid..., bins)
        new_qhist = fftconvolve(self.qhist, other.qhist, mode='full', axes=-1)
        
        # Clip small negative values from FFT floating-point artifacts
        new_qhist = np.maximum(new_qhist, 0)
        
        n_bins = len(self.qcenter)
        q0 = self.qcenter[0]
        
        # Alignment logic:
        # The first bin of 'full' convolution corresponds to value 2*q0.
        # We want to extract bins starting at value q0.
        # The index of value q0 in the 'full' array is (q0 - 2*q0) / bin_width = -q0 / bin_width.
        # This ensures that a delta function at value 0 (if present) acts as the identity.
        offset = int(round(-q0 / self.bin_width))
        
        if offset < 0:
            crop_start = 0
            fill_start = -offset
        else:
            crop_start = offset
            fill_start = 0
            
        final_qhist = np.zeros_like(self.qhist, dtype=float)
        
        crop_end = min(crop_start + n_bins - fill_start, new_qhist.shape[-1])
        fill_end = fill_start + (crop_end - crop_start)
        
        if fill_start < n_bins:
            final_qhist[..., fill_start:fill_end] = new_qhist[..., crop_start:crop_end]
            
        # Normalization to keep the y-scale (counts) roughly consistent.
        # Dividing by the sum of 'other' ensures that if 'other' is a delta-like
        # kernel, the total counts of 'self' are preserved.
        total_other = np.sum(other.qhist, axis=-1, keepdims=True)
        with np.errstate(divide='ignore', invalid='ignore'):
            final_qhist = np.where(total_other > 0, final_qhist / total_other, final_qhist)

        return ChargeHistogram(self.qcenter, final_qhist, self.bin_width)

    def normalize(self, density=True):
        """
        Normalizes the histogram.
        
        Args:
            density (bool): If True, the integral of the histogram will be 1 
                            (divided by sum and bin width).
                            If False, the sum of the bins will be 1.
                            
        Returns:
            ChargeHistogram: A new normalized histogram.
        """
        total = np.sum(self.qhist, axis=-1, keepdims=True)
        with np.errstate(divide='ignore', invalid='ignore'):
            new_qhist = np.where(total > 0, self.qhist / total, self.qhist).astype(float)
            if density:
                new_qhist /= self.bin_width
        return ChargeHistogram(self.qcenter, new_qhist, self.bin_width)

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

    def huber_location(self, loss='huber', scale=None):
        """
        Estimates the location (center) of the distribution using a robust loss function.
        
        Args:
            loss (str): The loss function to use ('huber', 'soft_l1', 'cauchy', 'arctan').
                        Defaults to 'huber'.
            scale (array-like, optional): The scale parameter for the loss function.
                                          Defaults to the square root of the Winsorized variance.
                                          
        Returns:
            np.ndarray: The estimated location for each point in the grid.
        """
        from scipy.optimize import minimize_scalar
        
        grid_shape = self.shape
        num_elements = int(np.prod(grid_shape)) if grid_shape else 1
        
        # Flatten for iteration
        flat_qhist = self.qhist.reshape(num_elements, -1)
        
        # Get initial guesses and scales
        x0_all = np.atleast_1d(self.median()).flatten()
        if scale is None:
            # Use a combination of IQR and Winsorized variance for scale
            # IQR / 1.349 is a robust estimate of sigma for Gaussian
            iqr_all = np.atleast_1d(self.iqr()).flatten() / 1.349
            wvar_all = np.sqrt(np.atleast_1d(self.winsorized_var())).flatten()
            
            # Take the smaller of the two, but ensure it's not zero
            scale_all = np.where((iqr_all > 0) & (iqr_all < wvar_all), iqr_all, wvar_all)
        else:
            scale_all = np.asanyarray(scale).flatten()
            if scale_all.size == 1:
                scale_all = np.full(num_elements, scale_all[0])

        # Loss functions rho(z) where z = (x - mu) / scale
        if loss == 'huber':
            rho = lambda z: np.where(np.abs(z) <= 1.0, 0.5 * z**2, np.abs(z) - 0.5)
        elif loss == 'soft_l1':
            rho = lambda z: 2 * (np.sqrt(1 + 0.5 * z**2) - 1)
        elif loss == 'cauchy':
            rho = lambda z: np.log(1 + 0.5 * z**2)
        elif loss == 'arctan':
            rho = lambda z: np.arctan(0.5 * z**2)
        else: # linear / squared
            rho = lambda z: 0.5 * z**2

        results = np.full(num_elements, np.nan)
        
        for i in range(num_elements):
            s_val = scale_all[i]
            if np.isnan(s_val) or s_val <= 0:
                s_val = self.bin_width
                
            x0 = x0_all[i]
            if np.isnan(x0):
                continue
            
            mask = flat_qhist[i] > 0
            if not np.any(mask):
                continue
                
            w = flat_qhist[i, mask]
            c = self.qcenter[mask]
            
            def objective(mu):
                return np.sum(w * rho((c - mu) / s_val))
            
            try:
                # Bracket around the median
                res = minimize_scalar(objective, bracket=(x0 - 5*s_val, x0, x0 + 5*s_val))
                results[i] = res.x
            except Exception:
                results[i] = x0
            
        return results.reshape(grid_shape) if grid_shape else results[0]

class ChargeSpectra:
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
        return f"<ChargeSpectra events={self.num_events}>"

    def normalize(self, density=True):
        """
        Normalizes all histograms in the container.
        
        Args:
            density (bool): If True, normalizes to a density (integral is 1).
                            If False, normalizes so the sum is 1.
        """
        return ChargeSpectra(
            self.num_events,
            self.pix.normalize(density=density),
            self.sipm.normalize(density=density),
            self.quabo.normalize(density=density),
            self.camera.normalize(density=density)
        )

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

    return ChargeSpectra(
        num_events=num_events,
        pix=ChargeHistogram(qcenter_pix, qhist_pix, w_pix),
        sipm=ChargeHistogram(qcenter_sipm, qhist_sipm, w_sipm),
        quabo=ChargeHistogram(qcenter_quabo, qhist_quabo, w_quabo),
        camera=ChargeHistogram(qcenter_camera, qhist_camera, w_camera)
    )

def calculate_charge_spectra(camera_images, gti_indexes=None, combine_gtis=True, 
                             qmin_pix=-512, qmax_pix=4096, downsample=False,
                             density=False):
    """
    Calculate the charge spectrum for each pixel from a CameraImages object.
    
    Args:
        camera_images (CameraImages): The loaded camera images and metadata.
        gti_indexes (list, optional): List of GTI indexes to include. If None, all are included.
        combine_gtis (bool): If True, combine all selected GTIs into one spectra object.
                            If False, return a dict of spectra objects keyed by GTI index.
        qmin_pix (int): Minimum charge for pixel histograms.
        qmax_pix (int): Maximum charge for pixel histograms.
        downsample (bool): If True, use larger bins for SiPM, Quabo, and camera histograms.
        density (bool): If True, normalize histograms to a density (integral is 1).

    Returns:
        ChargeSpectra or dict: The calculated spectra.
    """
    if gti_indexes is None:
        mask = np.ones(len(camera_images.gti_indexes), dtype=bool)
    else:
        mask = np.isin(camera_images.gti_indexes, gti_indexes)

    if combine_gtis:
        res = _build_spectra(camera_images.images[:, :, mask], 
                             qmin_pix=qmin_pix, qmax_pix=qmax_pix, downsample=downsample)
        if density:
            res = res.normalize(density=True)
        return res
    else:
        results = {}
        unique_gtis = np.unique(camera_images.gti_indexes[mask])
        for gti_idx in unique_gtis:
            gti_mask = (camera_images.gti_indexes == gti_idx)
            res = _build_spectra(camera_images.images[:, :, gti_mask],
                                 qmin_pix=qmin_pix, qmax_pix=qmax_pix, downsample=downsample)
            if density:
                res = res.normalize(density=True)
            results[gti_idx] = res
        return results

def apply_polynomial_pedestal_correction(camera_images, norder=0, quantiles=(0.15, 0.85)):
    """
    Fits a polynomial pedestal model to each pixel and subtracts it.
    The fit is performed independently for each GTI.
    Events outside the specified quantiles are excluded from the fit for robustness.
    
    Args:
        camera_images (CameraImages): The image container.
        norder (int): The order of the polynomial to fit (0=constant, 1=linear, etc.).
        quantiles (tuple): The (low, high) quantile limits for including data in the fit.
        
    Returns:
        CameraImages: A new container with pedestal-subtracted images.
    """
    def _correct_gti_images(gti_images):
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
        return gti_images.apply_pedestal_corrections(pcorr)
        
    return camera_images.map_gtis(_correct_gti_images)

def apply_constant_pedestal_correction(images, pedestal_calculator=None, **kwargs):
    """
    Calculates a constant pedestal correction value for each pixel and subtracts it.
    The calculation is performed independently for each GTI.
    
    Args:
        images (CameraImages): The image container.
        pedestal_calculator (str or callable, optional): 
            If a string, one of:
                - "quantile": Uses the specified quantile (default 0.5/median).
                - "huber", "cauchy", "soft_l1", "arctan": Uses the corresponding 
                  robust location estimator from ChargeHistogram.
            If None, defaults to "quantile".
            If a callable, it should take a CameraImages object and **kwargs 
            and return a (32, 32) array of pedestal values.
        **kwargs: Arguments passed to the calculator (e.g., quantile=0.5 or scale=1.0).
                  
    Returns:
        CameraImages: A new container with pedestal-subtracted images.
    """
    if pedestal_calculator is None or pedestal_calculator == "quantile":
        def internal_pedestal_calculator(gti_images, quantile=0.5, **unused_kwargs):
            charge_spectra = calculate_charge_spectra(gti_images, combine_gtis=True)
            return charge_spectra.pix.quantiles(quantile)
        pedestal_calculator = internal_pedestal_calculator
    elif isinstance(pedestal_calculator, str) and pedestal_calculator in ["huber", "cauchy", "soft_l1", "arctan"]:
        loss_name = pedestal_calculator
        def robust_pedestal_calculator(gti_images, **kwargs):
            charge_spectra = calculate_charge_spectra(gti_images, combine_gtis=True)
            return charge_spectra.pix.huber_location(loss=loss_name, **kwargs)
        pedestal_calculator = robust_pedestal_calculator

    def _correct_gti_images(gti_images):
        pedestal_val = pedestal_calculator(gti_images, **kwargs)
        pcorr = pedestal_val[:, :, np.newaxis]
        return gti_images.apply_pedestal_corrections(pcorr)
        
    return images.map_gtis(_correct_gti_images)
