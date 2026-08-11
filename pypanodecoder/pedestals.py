#!/usr/bin/env python3

# pedestals.py: Calculate charge spectra and pedestal values

# Author: Stephen Fegan <sfegan@llr.in2p3.fr> (2026-05-13)
# Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

import os
import json
import glob
import bisect
import time
import numpy as np

from .pcapdecoder import parse_time
from .eventbuilder import CameraImages, CameraEvent

class BaselineDatabase:
    """
    A database of Quabo pedestal baseline measurements.
    Loads measurements from quabo_ph_baseline.json files and maps UIDs to
    board locations using quabo_uids.json files.
    """
    def __init__(self, paths=None):
        """
        Args:
            paths (str, list, optional): Glob pattern(s) or list of files/directories
                                         to search for baseline files.
        """
        # board_loc -> [ (timestamp, baseline_array), ... ] sorted by timestamp
        self.baselines = {}
        if paths:
            self.load(paths)

    def load(self, paths):
        """
        Searches for and loads baseline measurements.

        Args:
            paths (str, list): Glob pattern(s) or list of files/directories.
        """
        if isinstance(paths, str):
            if any(c in paths for c in '*?['):
                files = glob.glob(paths, recursive=True)
            else:
                files = [paths]
        else:
            files = []
            for p in paths:
                if any(c in p for c in '*?['):
                    files.extend(glob.glob(p, recursive=True))
                else:
                    files.append(p)

        # Collect all quabo_ph_baseline.json files
        baseline_files = []
        for f in files:
            if os.path.isdir(f):
                baseline_files.extend(glob.glob(os.path.join(f, "**/quabo_ph_baseline.json"), recursive=True))
            elif f.endswith("quabo_ph_baseline.json"):
                baseline_files.append(f)

        for bfile in sorted(baseline_files):
            # Look for quabo_uids.json in the same directory
            bdir = os.path.dirname(bfile)
            ufile = os.path.join(bdir, "quabo_uids.json")

            if not os.path.exists(ufile):
                # Try parent directory as a fallback
                ufile = os.path.join(os.path.dirname(bdir), "quabo_uids.json")
                if not os.path.exists(ufile):
                    continue

            uid_map = self._load_uid_map(ufile)
            self._load_baseline_file(bfile, uid_map)

    def _load_uid_map(self, ufile):
        """Loads quabo_uids.json and returns a map of UID -> board_loc."""
        with open(ufile, 'r') as f:
            data = json.load(f)

        uid_to_loc = {}
        for dome in data.get('domes', []):
            for module in dome.get('modules', []):
                ip = module.get('ip_addr', '')
                try:
                    octets = [int(o) for o in ip.split('.')]
                    if len(octets) == 4:
                        # board_loc logic from pcapdecoder.py
                        board_loc_base = ((octets[2] & 0x03) << 8) | octets[3]
                        # Each module has 4 quabos, usually sequentially indexed                                                                                                                                           [352/1091]
                        for i, q in enumerate(module.get('quabos', [])):
                            uid = q.get('uid')
                            if uid:
                                uid_to_loc[uid] = board_loc_base + i
                except (ValueError, IndexError):
                    continue
        return uid_to_loc

    def _load_baseline_file(self, bfile, uid_map):
        """Loads a baseline file and adds its data to the database."""
        with open(bfile, 'r') as f:
            data = json.load(f)

        try:
            ts = parse_time(data.get('date'))
        except (ValueError, TypeError):
            # Use file modification time if date is missing or invalid
            ts = os.path.getmtime(bfile)

        for q in data.get('quabos', []):
            uid = q.get('uid')
            coefs = q.get('coefs')
            if uid in uid_map and coefs:
                loc = uid_map[uid]
                if loc not in self.baselines:
                    self.baselines[loc] = []

                # Check if this timestamp already exists for this board
                exists = False
                for existing_ts, _ in self.baselines[loc]:
                    if abs(existing_ts - ts) < 1e-3:
                        exists = True
                        break

                if not exists:
                    # Store as numpy array for easier manipulation
                    self.baselines[loc].append((ts, np.array(coefs)))
                    self.baselines[loc].sort(key=lambda x: x[0])

    def get_baseline(self, time, board_loc=None, telescope_id=None):
        """
        Retrieves the nearest baseline measurement before the given time.

        Args:
            time (float, str): UTC time (timestamp or ISO string).
            board_loc (int, optional): Specific board location.
            telescope_id (int, optional): Retrieve combined baseline image for a telescope.

        Returns:
            np.ndarray: 256-pixel array (if board_loc given) or
                        32x32 camera image (if telescope_id given).
                        Returns None if no baseline found.
        """
        ts_req = parse_time(time)

        if board_loc is not None:
            return self._get_board_baseline(ts_req, board_loc)

        if telescope_id is not None:
            # Reconstruct 32x32 image from the 4 quabos (0-3) of this telescope
            q_baselines = []
            for i in range(4):
                loc = (telescope_id << 2) | i
                q_baselines.append(self._get_board_baseline(ts_req, loc))

            # Check if we have at least one quabo
            if all(b is None for b in q_baselines):
                return None

            return CameraEvent.make_image(
                quabo0_pix=q_baselines[0],
                quabo1_pix=q_baselines[1],
                quabo2_pix=q_baselines[2],
                quabo3_pix=q_baselines[3]
            )

        return None

    def _get_board_baseline(self, ts, loc):
        if loc not in self.baselines or not self.baselines[loc]:
            return None

        # Find the index of the first baseline with timestamp > ts
        times = [b[0] for b in self.baselines[loc]]
        idx = bisect.bisect_right(times, ts)

        # We want the one BEFORE or AT ts, which is idx - 1
        if idx == 0:
            # No baseline before this time, return the earliest one anyway?
            # User says "return the nearest measured baseline that is before any given time"
            # If nothing is before, maybe we shouldn't return anything.
            return None

        return self.baselines[loc][idx - 1][1]

    def get_adjacent_times(self, time, board_loc=None, telescope_id=None):
        """
        Returns the times of the previous and next measurements.

        Returns:
            tuple: (prev_time, next_time). None if not available.
        """
        ts_req = parse_time(time)
        times = self.get_all_times(board_loc=board_loc, telescope_id=telescope_id)
        if not times:
            return (None, None)

        idx = bisect.bisect_left(times, ts_req)

        if idx < len(times) and abs(times[idx] - ts_req) < 1e-4:
            prev_idx = idx - 1
            next_idx = idx + 1
        else:
            prev_idx = idx - 1
            next_idx = idx

        prev_time = times[prev_idx] if prev_idx >= 0 else None
        next_time = times[next_idx] if next_idx < len(times) else None

        return (prev_time, next_time)

    def get_all_times(self, board_loc=None, telescope_id=None):
        """
        Returns a sorted list of all unique baseline measurement times.

        Args:
            board_loc (int, optional): Filter for a specific board.
            telescope_id (int, optional): Filter for a specific telescope (any of its 4 quabos).

        Returns:
            list: Sorted list of Unix timestamps (floats).
        """
        times = set()
        if board_loc is not None:
            if board_loc in self.baselines:
                times.update(b[0] for b in self.baselines[board_loc])
        elif telescope_id is not None:
            for i in range(4):
                loc = (telescope_id << 2) | i
                if loc in self.baselines:
                    times.update(b[0] for b in self.baselines[loc])
        else:
            for bl_list in self.baselines.values():
                times.update(b[0] for b in bl_list)

        return sorted(list(times))


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

    def downsample(self, factor):
        """
        Downsamples the histogram by combining adjacent bins.

        Args:
            factor (int): The number of bins to combine.

        Returns:
            ChargeHistogram: A new downsampled histogram.
        """
        if factor <= 1:
            return self

        n_bins = self.qhist.shape[-1]
        n_new_bins = n_bins // factor

        if n_new_bins == 0:
            # Return a single bin covering everything if factor > n_bins
            new_qcenter = np.array([self.qcenter.mean()])
            new_qhist = self.qhist.sum(axis=-1, keepdims=True)
            return ChargeHistogram(new_qcenter, new_qhist, self.bin_width * n_bins)

        truncated_n = n_new_bins * factor

        # Reshape and sum the histogram
        new_qhist = self.qhist[..., :truncated_n].reshape(self.shape + (n_new_bins, factor)).sum(axis=-1)

        # Reshape and average the bin centers
        new_qcenter = self.qcenter[:truncated_n].reshape(n_new_bins, factor).mean(axis=-1)

        return ChargeHistogram(new_qcenter, new_qhist, self.bin_width * factor)

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

    def to_log10(self, bins_per_decade=10, density=False):
        """
        Transforms the histogram into log10 space.

        Args:
            bins_per_decade (int): Number of bins per decade in log10 space.
            density (bool): If True, the resulting histogram will be a density
                            (divided by the log bin width).

        Returns:
            ChargeHistogram: A new histogram with log10(q) as bin centers.
        """
        # Filter positive qcenters
        mask = self.qcenter > 0
        if not np.any(mask):
            # Return an empty histogram if no positive values
            return ChargeHistogram(np.array([0.0]),
                                   np.zeros(self.shape + (1,)),
                                   1.0 / bins_per_decade)

        q_pos = self.qcenter[mask]
        log_q = np.log10(q_pos)

        log_bin_width = 1.0 / bins_per_decade

        # Define log10 edges covering the range of positive qcenters
        min_log = np.floor(np.min(log_q) * bins_per_decade) / bins_per_decade
        max_log = np.ceil(np.max(log_q) * bins_per_decade) / bins_per_decade
        log_edges = np.arange(min_log, max_log + 0.5 * log_bin_width, log_bin_width)

        if len(log_edges) < 2:
             log_edges = np.array([min_log, min_log + log_bin_width])

        new_qcenter = 0.5 * (log_edges[:-1] + log_edges[1:])

        # Original edges and CDF
        edges_orig, cdf_orig = self.cdf()

        # Interpolate CDF at target edges (in linear space)
        target_edges = 10**log_edges

        grid_shape = self.shape
        num_elements = int(np.prod(grid_shape)) if grid_shape else 1

        # Flatten for interpolation
        flat_cdf = cdf_orig.reshape(num_elements, -1)
        new_flat_cdf = np.zeros((num_elements, len(log_edges)))

        for i in range(num_elements):
            new_flat_cdf[i] = np.interp(target_edges, edges_orig, flat_cdf[i])

        # Probability in each log bin
        new_qhist_flat = np.diff(new_flat_cdf, axis=-1)

        # Multiply by original total sum to preserve scale (e.g. counts)
        total_sum = np.sum(self.qhist, axis=-1).reshape(num_elements, 1)
        new_qhist_flat *= total_sum

        if density:
            new_qhist_flat /= log_bin_width

        new_qhist = new_qhist_flat.reshape(grid_shape + (len(new_qcenter),))

        return ChargeHistogram(new_qcenter, new_qhist, log_bin_width)

    def cdf(self, reverse=False):
        """
        Returns the normalized cumulative distribution function (CDF) for all elements.

        Args:
            reverse (bool): If True, returns the reversed CDF (1 to 0),
                            representing the probability P(X >= x).

        Returns:
            tuple: (edges, cdf)
                edges: np.ndarray of bin edges (length N+1)
                cdf: np.ndarray of normalized cumulative counts (shape: grid_shape + (N+1))
        """
        # If qcenter are centers, left edges are qcenter - 0.5 * bin_width
        left_edges = self.qcenter - 0.5 * self.bin_width
        edges = np.append(left_edges, left_edges[-1] + self.bin_width)

        total_counts = np.sum(self.qhist, axis=-1)

        # Result CDF array
        cdf_shape = self.qhist.shape[:-1] + (self.qhist.shape[-1] + 1,)
        cdf = np.zeros(cdf_shape)

        good_mask = total_counts > 0
        if not good_mask.any():
            return edges, cdf

        if reverse:
            # P(X >= x)
            # cumsum from the right
            cumsum = np.cumsum(self.qhist[..., ::-1], axis=-1)[..., ::-1]
            cdf[good_mask, :-1] = cumsum[good_mask] / total_counts[good_mask, np.newaxis]
            # cdf[..., -1] is already 0
        else:
            # P(X <= x)
            # Cumulative sum along the range axis
            cumsum = np.cumsum(self.qhist, axis=-1)
            cdf[good_mask, 1:] = cumsum[good_mask] / total_counts[good_mask, np.newaxis]
            # cdf[..., 0] is already 0

        return edges, cdf

    def sf(self):
        """
        Returns the survival function (SF) for all elements in the histogram grid.
        SF is defined as 1 - CDF, representing P(X >= x).

        Returns:
            tuple: (edges, sf)
                edges: np.ndarray of bin edges (length N+1)
                sf: np.ndarray of normalized survival function values (shape: grid_shape + (N+1))
        """
        return self.cdf(reverse=True)

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

    def huber_location(self, loss='huber', scale=None, loss_parameter=1.0):
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
        loss_functions = {
            'huber':   lambda z: np.where(np.abs(z) <= loss_parameter, 0.5 * z**2, loss_parameter * np.abs(z) - 0.5 * loss_parameter**2),
            'soft_l1': lambda z: 2 * (np.sqrt(1 + 0.5 * z**2) - 1),
            'cauchy':  lambda z: np.log(1 + 0.5 * z**2),
            'arctan':  lambda z: np.arctan(0.5 * z**2),
            'linear':  lambda z: 0.5 * z**2
        }

        if not loss:
            loss = 'linear'

        if loss not in loss_functions:
            raise ValueError(f"Unknown loss function: {loss}. Supported: {list(loss_functions.keys())}")

        rho = loss_functions[loss]

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

    def huber_scale(self, loss='huber', location=None, loss_parameter=1.0):
        """
        Estimates the scale (standard deviation) of the distribution using a robust loss function.
        The result is calibrated to be consistent with the standard deviation for a Gaussian.

        Args:
            loss (str): The loss function to use ('huber', 'soft_l1', 'cauchy').
                        Defaults to 'huber'.
            location (array-like, optional): The location (center) parameter for the loss function.
                                             Defaults to the median.

        Returns:
            np.ndarray: The estimated scale for each point in the grid.
        """
        from scipy.optimize import minimize_scalar

        grid_shape = self.shape
        num_elements = int(np.prod(grid_shape)) if grid_shape else 1

        # Flatten for iteration
        flat_qhist = self.qhist.reshape(num_elements, -1)

        # Get initial locations and guesses for scale
        if location is None:
            mu_all = np.atleast_1d(self.median()).flatten()
        else:
            mu_all = np.asanyarray(location).flatten()
            if mu_all.size == 1:
                mu_all = np.full(num_elements, mu_all[0])

        # Robust initial scale guess
        iqr_all = np.atleast_1d(self.iqr()).flatten() / 1.349
        wvar_all = np.sqrt(np.atleast_1d(self.winsorized_var())).flatten()
        s0_all = np.where((iqr_all > 0) & (iqr_all < wvar_all), iqr_all, wvar_all)

        loss_functions = {
            'huber':   lambda z: np.where(np.abs(z) <= loss_parameter, 0.5 * z**2, loss_parameter * np.abs(z) - 0.5 * loss_parameter**2),
            'soft_l1': lambda z: 2 * (np.sqrt(1 + 0.5 * z**2) - 1),
            'cauchy':  lambda z: np.log(1 + 0.5 * z**2),
            'linear':  lambda z: 0.5 * z**2
        }

        from scipy.special import erf
        # Calibration constants beta = E[psi(Z)*Z] where Z ~ N(0, 1)
        # These constants ensure the estimator is consistent with Gaussian sigma.
        beta_map = {
            'huber':   erf(loss_parameter / np.sqrt(2.0)), # norm.cdf(k) - norm.cdf(-k)
            'soft_l1': 0.6808803445892556,
            'cauchy':  0.4842562477382487,
            'linear':  1.0
        }

        if not loss:
            loss = 'linear'

        if loss not in loss_functions:
            raise ValueError(f"Unknown loss function for scale estimation: {loss}. "
                             f"Supported: {list(loss_functions.keys())}")

        rho = loss_functions[loss]
        beta = beta_map[loss]

        results = np.full(num_elements, np.nan)

        for i in range(num_elements):
            s0 = s0_all[i]
            if np.isnan(s0) or s0 <= 0:
                s0 = self.bin_width

            mu = mu_all[i]
            if np.isnan(mu):
                continue

            mask = flat_qhist[i] > 0
            if not np.any(mask):
                continue

            w = flat_qhist[i, mask]
            c = self.qcenter[mask]
            w_sum = np.sum(w)

            # The objective function minimizes J(s) = sum(w * rho((c - mu) / s)) + beta * sum(w) * log(s)
            def objective(s):
                if s <= 0:
                    return np.inf
                return np.sum(w * rho((c - mu) / s)) + beta * w_sum * np.log(s)

            try:
                # Bounded optimization is more reliable for scale
                res = minimize_scalar(objective, bounds=(1e-6 * s0, 100 * s0), method='bounded')
                results[i] = res.x
            except Exception:
                results[i] = s0

        return results.reshape(grid_shape) if grid_shape else results[0]

    def huber_location_and_scale(self, loss='huber', loss_parameter=1.0):
        """
        Estimates the location (center) and scale (standard deviation) of the distribution
        simultaneously using a robust loss function.
        The scale result is calibrated to be consistent with the standard deviation for a Gaussian.

        Args:
            loss (str): The loss function to use ('huber', 'soft_l1', 'cauchy').
                        Defaults to 'huber'.

        Returns:
            tuple: (locations, scales)
                locations (np.ndarray): The estimated location for each point in the grid.
                scales (np.ndarray): The estimated scale for each point in the grid.
        """
        from scipy.optimize import minimize

        grid_shape = self.shape
        num_elements = int(np.prod(grid_shape)) if grid_shape else 1

        # Flatten for iteration
        flat_qhist = self.qhist.reshape(num_elements, -1)

        # Get initial locations and guesses for scale
        mu0_all = np.atleast_1d(self.median()).flatten()

        # Robust initial scale guess
        iqr_all = np.atleast_1d(self.iqr()).flatten() / 1.349
        wvar_all = np.sqrt(np.atleast_1d(self.winsorized_var())).flatten()
        s0_all = np.where((iqr_all > 0) & (iqr_all < wvar_all), iqr_all, wvar_all)

        loss_functions = {
            'huber':   lambda z: np.where(np.abs(z) <= loss_parameter, 0.5 * z**2, loss_parameter * np.abs(z) - 0.5 * loss_parameter**2),
            'soft_l1': lambda z: 2 * (np.sqrt(1 + 0.5 * z**2) - 1),
            'cauchy':  lambda z: np.log(1 + 0.5 * z**2),
            'linear':  lambda z: 0.5 * z**2
        }

        from scipy.special import erf
        # Calibration constants beta = E[psi(Z)*Z] where Z ~ N(0, 1)
        # These constants ensure the estimator is consistent with Gaussian sigma.
        beta_map = {
            'huber':   erf(loss_parameter / np.sqrt(2.0)), # norm.cdf(k) - norm.cdf(-k)
            'soft_l1': 0.6808803445892556,
            'cauchy':  0.4842562477382487,
            'linear':  1.0
        }

        if not loss:
            loss = 'linear'

        if loss not in loss_functions:
            raise ValueError(f"Unknown loss function for joint estimation: {loss}. "
                             f"Supported: {list(loss_functions.keys())}")

        rho = loss_functions[loss]
        beta = beta_map[loss]

        results_mu = np.full(num_elements, np.nan)
        results_s = np.full(num_elements, np.nan)

        for i in range(num_elements):
            s0 = s0_all[i]
            if np.isnan(s0) or s0 <= 0:
                s0 = self.bin_width
            
            mu0 = mu0_all[i]
            if np.isnan(mu0):
                continue

            mask = flat_qhist[i] > 0
            if not np.any(mask):
                continue

            w = flat_qhist[i, mask]
            c = self.qcenter[mask]
            w_sum = np.sum(w)

            # The objective function minimizes J(mu, s) = sum(w * rho((c - mu) / s)) + beta * sum(w) * log(s)
            def objective(params):
                mu, s = params
                if s <= 0:
                    return np.inf
                return np.sum(w * rho((c - mu) / s)) + beta * w_sum * np.log(s)

            try:
                res = minimize(objective, x0=[mu0, s0], bounds=[(None, None), (1e-6 * s0, 100 * s0)])
                results_mu[i] = res.x[0]
                results_s[i] = res.x[1]
            except Exception:
                results_mu[i] = mu0
                results_s[i] = s0

        res_mu = results_mu.reshape(grid_shape) if grid_shape else results_mu[0]
        res_s = results_s.reshape(grid_shape) if grid_shape else results_s[0]
        return res_mu, res_s

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

    def to_log10(self, bins_per_decade=10, density=False):
        """
        Transforms all histograms in the container into log10 space.

        Args:
            bins_per_decade (int): Number of bins per decade in log10 space.
            density (bool): If True, the resulting histograms will be densities.

        Returns:
            ChargeSpectra: A new container with log-transformed histograms.
        """
        return ChargeSpectra(
            self.num_events,
            self.pix.to_log10(bins_per_decade=bins_per_decade, density=density),
            self.sipm.to_log10(bins_per_decade=bins_per_decade, density=density),
            self.quabo.to_log10(bins_per_decade=bins_per_decade, density=density),
            self.camera.to_log10(bins_per_decade=bins_per_decade, density=density)
        )

    def downsample(self, factor):
        """
        Downsamples all histograms in the container.

        Args:
            factor (int): The number of bins to combine.

        Returns:
            ChargeSpectra: A new container with downsampled histograms.
        """
        return ChargeSpectra(
            self.num_events,
            self.pix.downsample(factor),
            self.sipm.downsample(factor),
            self.quabo.downsample(factor),
            self.camera.downsample(factor)
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

        # Convert gti_pcap_times from ns to seconds for numerical stability in polyfit
        t = gti_images.gti_pcap_times / 1e9
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

def calculate_spline_location_and_scale(camera_images, dtknot=600, nknot=None, loss='huber',
                                        max_iter=20, tol=1e-6, ridge=1e-8, loud=False, profile=False, loss_parameter=1.0):
    """
    Estimates the time-dependent location (center) and scale (standard deviation) of the
    distribution simultaneously using a robust loss function assuming a spline for location.
    The scale result is calibrated to be consistent with the standard deviation for a Gaussian.

    Performance notes vs. the scipy.optimize.minimize version:
      1. The cubic-spline basis matrix B (mu_interp = B @ knots) depends only on tknot_s and
         times_s, which are identical for every pixel, so it is built once instead of being
         reconstructed (as a CubicSpline object) on every objective evaluation of every pixel.
      2. The per-pixel optimization is solved via Iteratively Reweighted Least Squares (IRLS):
         each iteration is a weighted linear solve for the knot values (closed form) plus a 1D
         root-find for the scale, instead of a generic bounded quasi-Newton search with
         finite-difference gradients over (nknot + 1) parameters.
      3. All pixels are solved *simultaneously* rather than in a per-pixel Python loop:
         - The knot normal-equations (A, b) are built with two batched GEMMs using a
           precomputed pairwise-basis tensor G = einsum('ti,tj->tij', B, B), and solved with
           a single batched np.linalg.solve call instead of ~num_elements small (nknot,nknot)
           solves.
         - The scale root-find (previously a per-pixel scipy.optimize.brentq call) is replaced
           with a vectorized bisection over all active pixels at once. Bisection needs more
           iterations than Brent's method for the same tolerance, but collapsing ~1e6 tiny
           scalar calls (brentq's internal per-pixel evaluations) into ~60 calls on full-size
           arrays is a large net win at these array sizes. If the scale solve is later shown
           (via profiling) to dominate wall-clock, replacing bisection with a vectorized
           Illinois-modified regula falsi would cut ~60 iterations down to ~15-20.
         - Pixels are dropped out of the active set as soon as they converge, so later
           iterations only pay for pixels that still need work.

    Caveat: for convex losses (huber, soft_l1, linear) IRLS is a majorize-minimize algorithm
    and is guaranteed to decrease the same objective monotonically to a stationary point, so
    results should match the old optimizer to numerical precision. 'cauchy' is a non-convex
    (redescending) loss, so IRLS and BFGS can in principle converge to different local minima
    depending on the starting point -- worth spot-checking a subset of pixels against the old
    implementation before fully switching over, especially for 'cauchy'.

    Args:
        camera_images: images (CameraImages): The image container.
        dtknot: maximum time between knots
        nknot: force number of knots, overriding tknot if not None
        loss (str): The loss function to use ('huber', 'soft_l1', 'cauchy', 'linear').
                    Defaults to 'huber'.
        max_iter (int): maximum number of IRLS iterations per pixel.
        tol (float): relative convergence tolerance on knot values and scale between
                     IRLS iterations.
        ridge (float): small Tikhonov regularization added to the normal equations for
                       numerical stability (e.g. knots with little/no data nearby).

    Returns:
        tuple: (tknots, locations, scales)
            tknot (np.ndarray): Time of the knots. Shape: Nknot
            locations (np.ndarray): The estimated location for each point in the grid. Shape: 32x32xNknot
            scales (np.ndarray): The estimated scale for each point in the grid. Shape: 32x32
    """
    from scipy.interpolate import CubicSpline

    grid_shape = camera_images.images.shape[:-1]
    num_elements = int(np.prod(grid_shape)) if grid_shape else 1

    # Use seconds for knot spacing calculations to make tknot argument units clear
    times_s = camera_images.gti_pcap_times / 1e9
    t_start_s = np.min(times_s)
    t_end_s = np.max(times_s)
    if nknot is None:
        nknot = np.maximum(2, int(np.ceil((t_end_s - t_start_s) / dtknot)))
    tknot_s = np.linspace(t_start_s, t_end_s, nknot)

    if loud:
        import datetime
        dtknot_actual = int(np.round((t_end_s - t_start_s) / max(1, nknot - 1)))

        if camera_images.pcap_times is not None and len(camera_images.pcap_times) > 0:
            display_t_start_s = np.min(camera_images.pcap_times) / 1e9
        else:
            display_t_start_s = t_start_s

        start_time_str = datetime.datetime.fromtimestamp(display_t_start_s).strftime('%Y-%m-%d %H:%M:%S')
        prefix = f"{start_time_str} {nknot:3d}*{dtknot_actual:3d}s : "

    # Flatten for iteration
    flat_images = camera_images.images.reshape(num_elements, -1)
    n_time = flat_images.shape[1]

    # Get initial locations and robust guesses for scale
    mu0_all = np.median(flat_images, axis=-1)
    s0_all = (np.quantile(flat_images, 0.75, axis=-1) - np.quantile(flat_images, 0.25, axis=-1)) / 1.349

    # Build per-knot initial guesses by taking the median of samples within
    # a time window around each knot (window = +/- dtknot/2 seconds).
    mu0_knots = np.full((num_elements, nknot), np.nan)
    window = float(dtknot) / 2
    if window <= 0:
        total_dur = max(1e-12, (t_end_s - t_start_s))
        window = total_dur / max(2.0 * nknot, 1.0)

    for k in range(nknot):
        t_k = tknot_s[k]
        mask_t = (times_s >= (t_k - window)) & (times_s <= (t_k + window))
        if not np.any(mask_t):
            continue
        vals = flat_images[:, mask_t]
        mu0_knots[:, k] = np.median(vals, axis=-1)

    # Fill any knot windows that had no samples with the pixel's global median.
    nan_mask = np.isnan(mu0_knots)
    mu0_knots = np.where(nan_mask, mu0_all[:, None], mu0_knots)

    # --- Precompute the spline basis matrix once ---
    # CubicSpline(tknot_s, knots)(times_s) is linear in `knots` for fixed tknot_s/times_s,
    # so mu_interp = B @ knots. Build B by evaluating the spline for each unit knot vector,
    # once total, instead of reconstructing a CubicSpline per pixel per objective call.
    B = np.empty((n_time, nknot))
    for k in range(nknot):
        e_k = np.zeros(nknot)
        e_k[k] = 1.0
        B[:, k] = CubicSpline(tknot_s, e_k)(times_s)

    # Pairwise basis products, precomputed once. For a batch of pixels with IRLS weights
    # W (shape P x T), the weighted normal-equation matrices for every pixel in the batch
    # are obtained with a single GEMM: (W @ G).reshape(P, K, K) == B.T @ diag(w_p) @ B per row.
    G = np.einsum('ti,tj->tij', B, B).reshape(n_time, nknot * nknot)

    # IRLS weight functions w(z) = psi(z)/z (finite at z=0 for all four losses below).
    weight_functions = {
        'huber':   lambda z: np.divide(loss_parameter, np.abs(z), out=np.ones_like(z, dtype=float), where=np.abs(z) > loss_parameter),
        'soft_l1': lambda z: 1.0 / np.sqrt(1 + 0.5 * z ** 2),
        'cauchy':  lambda z: 1.0 / (1 + 0.5 * z ** 2),
        'linear':  lambda z: np.ones_like(z),
    }
    
    from scipy.special import erf
    # Calibration constants beta = E[psi(Z)*Z] where Z ~ N(0, 1); same values as the
    # original implementation, since the scale stationarity condition is
    # sum(psi(z)*z) = beta*qlen, i.e. sum(w(z)*z**2) = beta*qlen.
    beta_map = {
        'huber':   erf(loss_parameter / np.sqrt(2.0)),
        'soft_l1': 0.6808803445892556,
        'cauchy':  0.4842562477382487,
        'linear':  1.0,
    }

    if not loss:
        loss = 'linear'
    if loss not in weight_functions:
        raise ValueError(f"Unknown loss function for joint estimation: {loss}. "
                          f"Supported: {list(weight_functions.keys())}")

    wfun = weight_functions[loss]
    beta = beta_map[loss]
    ridge_matrix = ridge * np.eye(nknot)
    qlen = float(n_time)  # every pixel shares the same number of time samples

    def solve_scale_batched(R, s_guess, n_bisect=32, xtol_rel=1e-6):
        """
        Vectorized version of the original per-pixel brentq root-find:
            sum(w(r/s) * (r/s)**2) == beta * qlen  for s
        solved simultaneously for every row of R via bracket expansion + bisection.
        R: (P, T) residuals. s_guess: (P,) previous/initial scale estimate.
        Returns s: (P,) new scale estimate.

        This is the single most expensive call in the whole solve (it runs once, over the
        full active batch, on the first IRLS iteration), since each iteration evaluates
        wfun(z)*z**2 over the entire (P, T) residual array. The bracket spans 1e-3*s_guess
        to 1e3*s_guess, so matching brentq's old xtol=1e-6*s_guess needs log2(1e6/1e-3) ~ 30
        halvings -- not 60. An early-exit on bracket width additionally stops the loop as
        soon as every row is converged, rather than always spending the full budget.
        """
        def g(s):
            z = R / s[:, None]
            return np.sum(wfun(z) * z ** 2, axis=1) - beta * qlen

        lo = s_guess * 1e-3
        hi = s_guess * 1e3
        glo = g(lo)
        ghi = g(hi)

        for _ in range(12):
            bad = glo * ghi > 0
            if not np.any(bad):
                break
            lo = np.where(bad, lo * 0.1, lo)
            hi = np.where(bad, hi * 10.0, hi)
            glo = np.where(bad, g(lo), glo)
            ghi = np.where(bad, g(hi), ghi)

        # Pixels where a bracket still couldn't be found (shouldn't normally happen):
        # keep the previous scale estimate, same fallback as the original implementation.
        unbracketed = glo * ghi > 0

        a, b = lo.copy(), hi.copy()
        for _ in range(n_bisect):
            if np.all((b - a) < xtol_rel * s_guess):
                break
            mid = 0.5 * (a + b)
            gmid = g(mid)
            same_sign = (np.sign(gmid) == np.sign(glo))
            a = np.where(same_sign, mid, a)
            glo = np.where(same_sign, gmid, glo)
            b = np.where(same_sign, b, mid)

        s = 0.5 * (a + b)
        return np.where(unbracketed, s_guess, s)

    def solve_all_pixels(Q, mu0, s0, profile=False):
        """
        Batched IRLS over all valid pixels at once.
        Q: (P, T) data. mu0: (P, K) initial knot values. s0: (P,) initial scale.
        Returns knots (P, K), s (P,), converged (P,) bool, niter (P,) int.
        If profile=True, also returns a dict of cumulative wall-clock seconds spent in
        each phase (gemm_solve, scale_solve_exact, scale_solve_approx), so you can see
        which part actually dominates instead of guessing from FLOP counts.
        """
        P = Q.shape[0]
        knots = mu0.astype(float).copy()
        s = s0.astype(float).copy()
        converged = np.zeros(P, dtype=bool)
        niter = np.zeros(P, dtype=int)
        active = np.ones(P, dtype=bool)

        timings = {'gemm_solve': 0.0, 'scale_solve_exact': 0.0, 'scale_solve_approx': 0.0}

        for it in range(max_iter):
            idx = np.where(active)[0]
            if idx.size == 0:
                break

            Qa = Q[idx]
            Ka = knots[idx]
            Sa = s[idx]

            t0 = time.monotonic() if profile else None

            R = Qa - Ka @ B.T
            Z = R / Sa[:, None]
            W = wfun(Z)

            A = (W @ G).reshape(-1, nknot, nknot) + ridge_matrix
            b_vec = (W * Qa) @ B
            try:
                # NumPy's batched solve requires b to be explicitly (..., M, 1) here --
                # a plain (P, K) array is interpreted as one (P, K) matrix RHS, not a
                # batch of P length-K vectors, and raises a core-dimension mismatch.
                knots_new = np.linalg.solve(A, b_vec[..., None])[..., 0]
            except np.linalg.LinAlgError:
                # Rare: fall back to a per-row lstsq only for the rows that fail.
                knots_new = np.empty_like(b_vec)
                for j in range(A.shape[0]):
                    try:
                        knots_new[j] = np.linalg.solve(A[j], b_vec[j])
                    except np.linalg.LinAlgError:
                        knots_new[j], *_ = np.linalg.lstsq(A[j], b_vec[j], rcond=None)

            if profile:
                t1 = time.monotonic()
                timings['gemm_solve'] += t1 - t0

            R_new = Qa - knots_new @ B.T

            if False: # it == 0:
                s_new = solve_scale_batched(R_new, Sa)
                if profile:
                    t2 = time.monotonic()
                    timings['scale_solve_exact'] += t2 - t1
            else:
                s_step = np.sqrt(np.sum(W * R_new ** 2, axis=1) / (beta * qlen))
                s_new = np.clip(s_step, 0.5 * Sa, 2.0 * Sa)
                if profile:
                    t2 = time.monotonic()
                    timings['scale_solve_approx'] += t2 - t1

            d_knots = np.max(np.abs(knots_new - Ka), axis=1) / np.maximum(1.0, np.max(np.abs(Ka), axis=1))
            d_s = np.abs(s_new - Sa) / Sa

            knots[idx] = knots_new
            s[idx] = s_new
            niter[idx] += 1

            newly_done = (d_knots < tol) & (d_s < tol)
            done_idx = idx[newly_done]
            converged[done_idx] = True
            active[done_idx] = False

            if loud:
                n_active = int(np.sum(active))
                hashes = '#' * min(32, int(32 * (P - n_active) / max(P, 1)))
                print(f"\r{prefix}{hashes:<32} (iter {it+1:2d}/{max_iter}, active {n_active}/{P})",
                      end='', flush=True)

        if profile:
            return knots, s, converged, niter, timings
        return knots, s, converged, niter

    results_mu = np.full((num_elements, nknot), np.nan)
    results_s = np.full(num_elements, np.nan)

    # Only pixels with a sane initial guess are worth solving; everything else stays NaN,
    # matching the original per-pixel skip logic.
    valid = (~np.isnan(s0_all)) & (s0_all > 0) & (~np.isnan(mu0_all))
    valid_idx = np.where(valid)[0]

    num_pix_converged = 0
    num_iter = 0
    start_time = time.monotonic()

    if valid_idx.size > 0:
        Q_valid = flat_images[valid_idx]
        mu0_valid = mu0_knots[valid_idx]
        s0_valid = s0_all[valid_idx]

        try:
            if profile:
                knots_v, s_v, converged_v, niter_v, timings = solve_all_pixels(
                    Q_valid, mu0_valid, s0_valid, profile=True)
                total = sum(timings.values())
                print("\n---- phase breakdown ----")
                for k, v in timings.items():
                    pct = 100 * v / total if total > 0 else 0.0
                    print(f"  {k:20s} {v:8.3f}s  ({pct:5.1f}%)")
                print(f"  {'sum':20s} {total:8.3f}s")
                print("--------------------------")
            else:
                knots_v, s_v, converged_v, niter_v = solve_all_pixels(Q_valid, mu0_valid, s0_valid)
            results_mu[valid_idx, :] = knots_v
            results_s[valid_idx] = s_v
            num_pix_converged = int(np.sum(converged_v))
            num_iter = int(np.sum(niter_v))
        except Exception:
            import traceback
            traceback.print_exc()

    if loud:
        end_time = time.monotonic()
        hashes = '#' * 32
        print(f"\r{prefix}{hashes} (Converged: {num_pix_converged}, Time: {end_time - start_time:.2f}s)", flush=True)

    res_mu = results_mu.reshape(np.append(grid_shape, nknot)) if grid_shape else results_mu[0]
    res_s = results_s.reshape(grid_shape) if grid_shape else results_s[0]
    tknot_ns = (tknot_s * 1e9)
    return tknot_ns, res_mu, res_s

def apply_precomputed_spline_pedestal_correction(camera_images, tknot_ns, knots):
    """
    Applies a precomputed spline pedestal correction to each pixel's image values.
 
    Args:
        camera_images (CameraImages): The image container.
        tknot_ns (np.ndarray): Time of the knots in nanoseconds.
        knots (np.ndarray): 3D array of shape (32, 32, Nknot) containing knot values.

    Returns:
        CameraImages: A new object with adjusted images.
    """
    from scipy.interpolate import CubicSpline
    t_sec = camera_images.gti_pcap_times / 1e9
    tknot_s = tknot_ns / 1e9
    
    cs = CubicSpline(tknot_s, knots, axis=-1)
    res = cs(t_sec)
    
    return CameraImages(
        images=camera_images.images - res,
        event_times=camera_images.event_times,
        pcap_times=camera_images.pcap_times,
        gti_indexes=camera_images.gti_indexes,
        gti_pcap_times=camera_images.gti_pcap_times,
        quabo_masks=camera_images.quabo_masks,
        gtis=camera_images.gtis,
        events=camera_images.events,
        filter=dict(camera_images.filter),
        source=dict(camera_images.source),
        quabo_pcap_time=camera_images.quabo_pcap_time,
        quabo_event_time=camera_images.quabo_event_time
    )

def apply_spline_pedestal_correction(camera_images, dtknot=600, nknot=None, loss='huber',
                                     max_iter=20, tol=1e-6, ridge=1e-8, loud=False, loss_parameter=1.0):
    """
    Calculates a spline pedestal model for each pixel and subtracts it.
    The fit is performed independently for each GTI.

    Args:
        camera_images (CameraImages): The image container.
        dtknot (float): maximum time between knots
        nknot (int): force number of knots, overriding tknot if not None
        loss (str): The loss function to use ('huber', 'soft_l1', 'cauchy', 'linear').
                    Defaults to 'huber'.
        max_iter (int): maximum number of IRLS iterations per pixel.
        tol (float): relative convergence tolerance.
        ridge (float): small Tikhonov regularization.

    Returns:
        CameraImages: A new container with pedestal-subtracted images.
    """
    def _correct_gti_images(gti_images):
        tknot_ns, locations, scales = calculate_spline_location_and_scale(
            gti_images, dtknot=dtknot, nknot=nknot, loss=loss,
            max_iter=max_iter, tol=tol, ridge=ridge, loud=loud, loss_parameter=loss_parameter
        )
        return apply_precomputed_spline_pedestal_correction(gti_images, tknot_ns, locations)

    return camera_images.map_gtis(_correct_gti_images)