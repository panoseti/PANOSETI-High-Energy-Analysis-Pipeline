#!/usr/bin/env python3

# event_rate.py: Generate event rate and delta-t plots for PANOSETI

# Author: Stephen Fegan <sfegan@llr.in2p3.fr> (2026-05-17)
# Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime
import math
from scipy.optimize import curve_fit

from ..eventbuilder import CameraImages

def plot_event_rate(camera_images, bin_width_min=1.0, subplots=False, figsize=(10, 6), uttime=False, clip=False, gti_labels=None, fig=None, axes=None, **kwargs):
    """
    Plots the event rate in Hz as a function of GTI event time.

    The rate is plotted as a broken histogram (horizontal lines) with Poisson error bars.
    Multiple GTIs are plotted overlaid with a legend, or in separate subplots.

    Args:
        camera_images (CameraImages): The images and metadata container.
        bin_width_min (float): Bin width in minutes for rate calculation.
        subplots (bool): If True, each GTI is plotted in its own subplot with a shared X-axis.
        figsize (tuple): Size of the figure (width, height).
        uttime (bool): If True, plot against absolute UT time (HH:MM) instead of seconds
                       since GTI start.
        clip (bool): If True, set Y-axis limits based on 5th and 95th percentiles of all rates
                     (expanded by 2.0x) and mark out-of-range points with triangles.
        **kwargs: Additional keyword arguments passed to hlines and errorbar (e.g., alpha, linewidth).

    Returns:
        tuple: (matplotlib.figure.Figure, list of matplotlib.axes.Axes)
    """
    unique_gtis = camera_images.unique_gti_indexes
    num_gtis = len(unique_gtis)

    if num_gtis == 0:
        print("No GTIs found in the provided camera images.")
        return None, None

    bin_width_sec = bin_width_min * 60.0
    # Format bin width to remove .0 if it's an integer
    bw_str = f"{bin_width_min:g}"
    ylabel = f"Rate (Hz) [{bw_str}-min bins]"

    if fig is None and axes is None:
        if subplots:
            fig, axes = plt.subplots(num_gtis, 1, sharex=True,
                                     figsize=figsize,
                                     squeeze=False)
            axes = axes.flatten()
            fig.tight_layout(rect=[0.05, 0, 1, 1])
            fig.subplots_adjust(hspace=0)
        else:
            fig, ax = plt.subplots(figsize=figsize)
            axes = [ax] * num_gtis
            fig.tight_layout()
    else:
        if not isinstance(axes, (list, tuple, np.ndarray)):
            axes = [axes] * num_gtis
        elif len(axes) == 1 and num_gtis > 1:
            axes = list(axes) * num_gtis
        if fig is None and len(axes) > 0:
            fig = axes[0].figure

    # Pre-calculate all rates if clipping is requested
    all_values_for_clipping = []
    gti_data_list = []

    for i, gti_idx in enumerate(unique_gtis):
        mask = (camera_images.gti_indexes == gti_idx)

        if uttime:
            # Use pcap_times for plotting since event_time is unreliable.
            # Convert nanoseconds to seconds.
            times = (camera_images.pcap_times[mask] / 1e9) % 86400
        else:
            # gti_pcap_times is now in nanoseconds, convert to seconds for plotting
            times = camera_images.gti_pcap_times[mask] / 1e9

        if len(times) == 0:
            gti_data_list.append(None)
            continue

        # Determine bins aligned to bin_width_sec
        min_t = np.min(times)
        max_t = np.max(times)
        bin_start = np.floor(min_t / bin_width_sec) * bin_width_sec
        bins = np.arange(bin_start, max_t + bin_width_sec, bin_width_sec)

        # Calculate histogram and rates
        counts, bin_edges = np.histogram(times, bins=bins)
        
        # Calculate exposure for each bin based on the first and last event time
        exposures = np.clip(bin_edges[1:], None, max_t) - np.clip(bin_edges[:-1], min_t, None)
        
        # Avoid division by zero
        valid = exposures > 0
        rates = np.zeros_like(counts, dtype=float)
        errors = np.zeros_like(counts, dtype=float)
        rates[valid] = counts[valid] / exposures[valid]
        errors[valid] = np.sqrt(counts[valid]) / exposures[valid]

        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0

        if clip:
            all_values_for_clipping.extend(rates + errors)
            all_values_for_clipping.extend(np.maximum(0, rates - errors))

        gti_data_list.append({
            'rates': rates, 'errors': errors, 'bin_edges': bin_edges,
            'bin_centers': bin_centers, 'gti_idx': gti_idx
        })

    y_min, y_max = None, None
    if clip and all_values_for_clipping:
        p5, p95 = np.percentile(all_values_for_clipping, [5, 95])
        mid = (p5 + p95) / 2
        half_width = (p95 - p5) * 2.0 / 2
        y_min, y_max = max(0.0, mid - half_width), mid + half_width
        # Ensure range is non-zero
        if y_max <= y_min:
            y_max = y_min + 1.0

    # Get colors from default cycle
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']

    # Default plotting parameters
    plot_kwargs = {'fmt': '.', 'capsize': 2}
    plot_kwargs.update(kwargs)

    for i, data in enumerate(gti_data_list):
        if data is None: continue
        cur_ax = axes[i]
        # Use a single color for subplots, cycle for overlaid
        color = colors[0] if subplots else colors[i % len(colors)]
        rates, errors, bin_edges, bin_centers = data['rates'], data['errors'], data['bin_edges'], data['bin_centers']
        gti_idx = data['gti_idx']

        # Construct label
        label = f"GTI {gti_idx}"
        if camera_images.gtis:
            try:
                if gti_idx in camera_images.gtis:
                    gti_info = camera_images.gtis[gti_idx]
                    start = gti_info.get('start', gti_info.get('stop', 'Unknown'))

                    dt = None
                    if isinstance(start, str):
                        try:
                            # Try ISO format
                            dt = datetime.datetime.fromisoformat(start.replace('Z', '+00:00'))
                        except ValueError:
                            pass
                    elif isinstance(start, (int, float)) and math.isfinite(start):
                        dt = datetime.datetime.fromtimestamp(start, tz=datetime.timezone.utc)

                    if dt:
                        fmt_str = "%Y-%m-%d" if uttime else "%Y-%m-%d %H:%M:%S"
                        label = f"{dt.strftime(fmt_str)} ({label})"
                    else:
                        label = f"{start} ({label})"
            except (AttributeError, IndexError, TypeError, ValueError):
                pass

        if gti_labels is not None and gti_idx in gti_labels:
            label = gti_labels[gti_idx]

        if uttime:
            # Convert bin edges and centers to datetime objects (time of day)
            # We use UTC timestamps % 86400 which are seconds since midnight 1970-01-01
            x_centers = [datetime.datetime.fromtimestamp(t, tz=datetime.timezone.utc) for t in bin_centers]
        else:
            x_centers = bin_centers

        # Plot error bars with markers. Passing color if not in kwargs.
        current_plot_kwargs = plot_kwargs.copy()
        if 'color' not in current_plot_kwargs and 'ecolor' not in current_plot_kwargs:
            current_plot_kwargs['color'] = color

        cur_ax.errorbar(x_centers, rates, yerr=errors,
                        label=label if not subplots else None,
                        **current_plot_kwargs)

        if clip:
            cur_ax.set_ylim(y_min, y_max)
            # Mark out-of-range points
            hi_mask = rates > y_max
            lo_mask = rates < y_min
            if np.any(hi_mask):
                cur_ax.scatter(np.array(x_centers)[hi_mask], [y_max*0.95 + y_min*0.05]*np.sum(hi_mask),
                               marker='^', color='red', s=20, zorder=5)
            if np.any(lo_mask):
                cur_ax.scatter(np.array(x_centers)[lo_mask], [y_min*0.95 + y_max*0.05]*np.sum(lo_mask),
                               marker='v', color='red', s=20, zorder=5)

        if subplots:
            # Place label inside the axis
            cur_ax.text(0.01, 0.95, label, transform=cur_ax.transAxes,
                        va='top', ha='left', fontsize=9, fontweight='bold',
                        bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1))

    if not subplots:
        if num_gtis > 1:
            axes[0].legend()
        axes[0].set_ylabel(ylabel)
    else:
        # Single common Y-axis label for subplots
        if fig is not None:
            fig.supylabel(ylabel)
        else:
            axes[len(axes)//2].set_ylabel(ylabel)

    if uttime:
        axes[-1].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        axes[-1].set_xlabel("UT Time (HH:MM)")
    else:
        axes[-1].set_xlabel("Time since GTI start (s)")

    return fig, axes

def plot_delta_t(camera_images, combine_gtis=False, semilog=False, density=True, fit=False, num_bins=100, figsize=(10, 6), time_type='event', gti_labels=None,  **kwargs):
    """
    Plots the distribution of times between consecutive events (delta_t).

    The distribution is plotted as log-log (default) or semilog-y.
    By default, it bins by log(delta_t) (equally spaced in log space).
    The 'semilog' option bins by delta_t (equally spaced in linear space).
    Optionally fits an exponential model to the data.

    Args:
        camera_images (CameraImages): The images and metadata container.
        combine_gtis (bool): If True, combine delta_t from all GTIs into one distribution.
        semilog (bool): If True, bin by dt linearly and plot semilog-y.
                        If False (default), bin by log(dt) and plot log-log.
        density (bool): If True, normalize the distribution (integral = 1).
        fit (bool): If True, fit an exponential model to each distribution.
        num_bins (int): Number of bins for the histogram.
        figsize (tuple): Size of the figure (width, height).
        time_type (str): Type of time to use for intervals: 'pcap' (default) or 'event'.
        **kwargs: Additional keyword arguments passed to axes.stairs.

    Returns:
        tuple: (matplotlib.figure.Figure, matplotlib.axes.Axes)
    """
    unique_gtis = camera_images.unique_gti_indexes

    gti_dts = {}
    all_dts_list = []

    for gti_idx in unique_gtis:
        mask = (camera_images.gti_indexes == gti_idx)
        
        if time_type.lower() == 'event':
            # Maintain int64 precision through sort and diff
            times_ns = np.sort(camera_images.event_times[mask])
        else:
            # Maintain int64 precision through sort and diff
            times_ns = np.sort(camera_images.pcap_times[mask])

        if len(times_ns) < 2:
            continue
        
        # Calculate diff in nanoseconds first, then convert to float seconds
        dt = np.diff(times_ns) / 1e9
        dt = dt[dt > 0] # Ensure positive intervals
        if len(dt) > 0:
            gti_dts[gti_idx] = dt
            all_dts_list.append(dt)

    if not gti_dts:
        print("No delta_t values to plot.")
        return None, None

    all_dts = np.concatenate(all_dts_list)

    # Calculate global binning range
    min_dt, max_dt = np.min(all_dts), np.max(all_dts)
    if semilog:
        bins = np.linspace(min_dt, max_dt, num_bins + 1)
        ylabel = r"dN/dt [s$^{-1}$]" if density else "Counts"
    else:
        log_min, log_max = np.log(min_dt), np.log(max_dt)
        bins = np.linspace(log_min, log_max, num_bins + 1)
        ylabel = r"dN/dlog($\Delta t$) [1]" if density else "Counts"

    if combine_gtis:
        data_to_plot = {'Combined': all_dts}
    else:
        data_to_plot = gti_dts

    fig, ax = plt.subplots(figsize=figsize)

    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']

    all_y_values = []
    min_normalized_y = 1.0 # Default fallback

    for i, (gti_idx, dts) in enumerate(data_to_plot.items()):
        color = colors[i % len(colors)]

        if semilog:
            counts, bin_edges = np.histogram(dts, bins=bins)
            x_plot = bin_edges
            x_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0
            bw = np.diff(bins)
        else:
            counts, bin_edges = np.histogram(np.log(dts), bins=bins)
            x_plot = np.exp(bin_edges)
            x_centers = np.exp((bin_edges[:-1] + bin_edges[1:]) / 2.0)
            bw = np.diff(bins)

        y_plot = counts.astype(float)
        # Calculate normalization factor: 1.0 / (total * bin_width)
        # This is used to find the value of "0.5 events" in the density space
        norm_factor = 1.0 / (np.sum(counts) * bw) if density else 1.0

        if density:
            y_plot = y_plot * norm_factor

        # Set min_normalized_y to the value of 0.5 events for this GTI
        # We take the minimum across all GTIs to ensure we see everything
        current_gti_min = 0.5 * np.min(norm_factor)
        if i == 0:
            min_normalized_y = current_gti_min
        else:
            min_normalized_y = min(min_normalized_y, current_gti_min)

        all_y_values.extend(y_plot[counts > 0])

        # Construct base label
        label = f"GTI {gti_idx}" if gti_idx != 'Combined' else 'Combined GTIs'
        if gti_idx != 'Combined' and camera_images.gtis:
            try:
                if gti_idx in camera_images.gtis:
                    gti_info = camera_images.gtis[gti_idx]
                    start = gti_info.get('start', gti_info.get('stop', 'Unknown'))
                    dt_start = None
                    if isinstance(start, str):
                        try:
                            dt_start = datetime.datetime.fromisoformat(start.replace('Z', '+00:00'))
                        except ValueError: pass
                    elif isinstance(start, (int, float)) and math.isfinite(start):
                        dt_start = datetime.datetime.fromtimestamp(start, tz=datetime.timezone.utc)

                    if dt_start:
                        label = f"{dt_start.strftime('%Y-%m-%d %H:%M:%S')} ({label})"
                    else:
                        label = f"{start} ({label})"
            except: pass

        if gti_labels is not None and gti_idx in gti_labels:
            label = gti_labels[gti_idx]

        if fit:
            # Fit exponential model. Ignore empty bins.
            mask = counts > 0
            if np.sum(mask) > 2:
                xx = x_centers[mask]
                yy = y_plot[mask]
                # Poisson errors: sigma = counts / (total * bw) / sqrt(counts) = yy / sqrt(counts)
                sigmas = yy / np.sqrt(counts[mask])

                if semilog:
                    def model(x, A, lam): return A * np.exp(-lam * x)
                    p0 = [yy[0], 1.0 / np.mean(dts)]
                else:
                    def model(x, A, lam): return A * x * np.exp(-lam * x)
                    p0 = [np.max(yy) * np.exp(1) * (1.0/np.mean(dts)), 1.0 / np.mean(dts)]

                try:
                    popt, _ = curve_fit(model, xx, yy, p0=p0, sigma=sigmas)
                    label = f"{label} (Rate: {popt[1]:.2f} Hz)"
                    x_smooth = np.logspace(np.log10(min_dt), np.log10(max_dt), 200) if not semilog else np.linspace(min_dt, max_dt, 200)
                    y_smooth = model(x_smooth, *popt)
                    # Clip smooth curve to sensible limits for plotting
                    y_smooth = np.maximum(y_smooth, min_normalized_y * 0.1)
                    ax.plot(x_smooth, y_smooth, color=color, alpha=0.8, linestyle='--')
                except Exception as e:
                    print(f"Fit failed for {label}: {e}")

        current_kwargs = kwargs.copy()
        if 'color' not in current_kwargs:
            current_kwargs['color'] = color
        if 'label' not in current_kwargs:
            current_kwargs['label'] = label

        ax.stairs(y_plot, x_plot, **current_kwargs)

    ax.set_yscale('log')
    if all_y_values:
        y_max = np.max(all_y_values) * 2.0
        ax.set_ylim(min_normalized_y, y_max)

    ax.set_xlabel(rf'$\Delta t_{{{time_type.lower()}}}$ (s)')
    ax.set_ylabel(ylabel)

    if semilog:
        ax.set_xscale('linear')
    else:
        ax.set_xscale('log')
        def safe_reciprocal(x):
            with np.errstate(divide='ignore', invalid='ignore'):
                return 1.0 / x
        try:
            secax = ax.secondary_xaxis('top', functions=(safe_reciprocal, safe_reciprocal))
            secax.set_xlabel('Frequency (Hz)')
        except (AttributeError, ValueError): pass

    if len(data_to_plot) > 1 or fit:
        ax.legend()

    fig.tight_layout()
    return fig, ax
