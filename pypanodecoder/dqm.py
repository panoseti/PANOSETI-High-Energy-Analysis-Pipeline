#!/usr/bin/env python3

# dqm.py: Data Quality Monitoring tools for PANOSETI

# Author: Gemini CLI (2026-05-17)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime
import os
import sys

# Ensure we can import from the current directory
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.append(current_dir)

try:
    from eventbuilder import PanosetiCameraImages
except ImportError:
    # Fallback if imported from outside the package
    from .eventbuilder import PanosetiCameraImages

def plot_event_rate(camera_images, bin_width_min=1.0, gtis=None, subplots=False, figsize=(10, 6), uttime=False, clip=False, **kwargs):
    """
    Plots the event rate in Hz as a function of GTI event time.
    
    The rate is plotted as a broken histogram (horizontal lines) with Poisson error bars.
    Multiple GTIs are plotted overlaid with a legend, or in separate subplots.

    Args:
        camera_images (PanosetiCameraImages): The images and metadata container.
        bin_width_min (float): Bin width in minutes for rate calculation.
        gtis (list, optional): Original list of GTI dictionaries. If provided, 
                               GTI start times will be included in the legend/titles.
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

    # Pre-calculate all rates if clipping is requested
    all_values_for_clipping = []
    gti_data_list = []

    for i, gti_idx in enumerate(unique_gtis):
        mask = (camera_images.gti_indexes == gti_idx)
        
        if uttime:
            times = camera_images.event_times[mask] % 86400
        else:
            times = camera_images.gti_event_times[mask]
        
        if len(times) == 0:
            gti_data_list.append(None)
            continue
            
        # Determine bins
        min_t = np.min(times)
        max_t = np.max(times)
        bins = np.arange(min_t, max_t + bin_width_sec, bin_width_sec)
        
        # Calculate histogram and rates
        counts, bin_edges = np.histogram(times, bins=bins)
        rates = counts / bin_width_sec
        errors = np.sqrt(counts) / bin_width_sec
        
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
        if gtis is not None:
            try:
                # Use gti_idx for lookup if it's a valid index for the list
                if isinstance(gti_idx, (int, np.integer)) and 0 <= gti_idx < len(gtis):
                    gti_info = gtis[gti_idx]
                    start = gti_info.get('start', gti_info.get('stop', 'Unknown'))
                    
                    dt = None
                    if isinstance(start, str):
                        try:
                            # Try ISO format
                            dt = datetime.datetime.fromisoformat(start.replace('Z', '+00:00'))
                        except ValueError:
                            pass
                    elif isinstance(start, (int, float)):
                        dt = datetime.datetime.fromtimestamp(start, tz=datetime.timezone.utc)
                    
                    if dt:
                        fmt_str = "%Y-%m-%d" if uttime else "%Y-%m-%d %H:%M:%S"
                        label = f"{dt.strftime(fmt_str)} ({label})"
                    else:
                        label = f"{start} ({label})"
            except (AttributeError, IndexError, TypeError, ValueError):
                pass
        
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
                cur_ax.scatter(np.array(x_centers)[hi_mask], [y_max]*np.sum(hi_mask), 
                               marker='^', color='red', s=20, zorder=5)
            if np.any(lo_mask):
                cur_ax.scatter(np.array(x_centers)[lo_mask], [y_min]*np.sum(lo_mask), 
                               marker='v', color='red', s=20, zorder=5)

        if subplots:
            # Place label inside the axis
            cur_ax.text(0.01, 0.95, label, transform=cur_ax.transAxes, 
                        va='top', ha='left', fontsize=9, fontweight='bold',
                        bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1))
            
    if not subplots:
        axes[0].legend()
        axes[0].set_ylabel(ylabel)
    else:
        # Single common Y-axis label for subplots
        fig.supylabel(ylabel)
        
    if uttime:
        axes[-1].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        axes[-1].set_xlabel("UT Time (HH:MM)")
    else:
        axes[-1].set_xlabel("Time since GTI start (s)")
    
    return fig, axes

if __name__ == "__main__":
    # Minimal self-test with dummy data
    print("Running self-test for plot_event_rate...")
    
    # Create dummy PanosetiCameraImages
    # GTI 0: 100 events in 100 seconds
    # GTI 1: 50 events in 100 seconds
    times0 = np.random.uniform(0, 100, 100)
    times1 = np.random.uniform(0, 100, 50)
    
    dummy_images = PanosetiCameraImages(
        images=np.zeros((32, 32, 150)),
        event_times=np.concatenate([times0, times1 + 1000]),
        gti_indexes=np.array([0]*100 + [1]*50),
        gti_event_times=np.concatenate([times0, times1]),
        quabo_masks=np.zeros(150, dtype=int)
    )
    
    # Test overlaid
    fig1, _ = plot_event_rate(dummy_images, bin_width_min=0.2)
    plt.close(fig1)
    print("Overlaid plot test passed.")
    
    # Test subplots
    fig2, _ = plot_event_rate(dummy_images, bin_width_min=0.2, subplots=True)
    plt.close(fig2)
    print("Subplots test passed.")
    
    # Test clip and uttime
    fig3, _ = plot_event_rate(dummy_images, uttime=True, clip=True)
    plt.close(fig3)
    print("UT Time and Clip test passed.")
    
    print("dqm.py self-test complete.")
