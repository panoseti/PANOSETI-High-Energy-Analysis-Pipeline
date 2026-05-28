#!/usr/bin/env python3

# charge.py: Diagnostic functions for charge spectra

# Author: Gemini CLI
# Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

import numpy as np
import matplotlib.pyplot as plt
from ..pedestals import ChargeHistogram

def plot_charge_spectra(qspecs, fig=None, title=None, pixel=None, **kwargs):
    """
    Plots charge spectra in three different ways: linear, log-y, and log-log.
    
    Args:
        qspecs: A single ChargeHistogram or a dict of {label: ChargeHistogram}.
        fig (matplotlib.figure.Figure, optional): Figure to plot into.
        title (str, optional): Title for the plot.
        pixel (tuple, optional): (j, i) coordinates if qspecs contains a grid.
        **kwargs: Additional arguments (e.g., color, alpha) passed to plot functions.
        
    Returns:
        tuple: (matplotlib.figure.Figure, list of matplotlib.axes.Axes)
    """
    if isinstance(qspecs, ChargeHistogram):
        qspecs_dict = {None: qspecs}
    else:
        qspecs_dict = qspecs

    if fig is None:
        fig = plt.figure(figsize=(12, 4))
    
    # Check if axes already exist in the figure
    if fig.axes:
        axes = fig.axes
        if len(axes) < 3:
            # Not enough axes, clear and recreate
            fig.clf()
            axes = [fig.add_subplot(1, 3, i+1) for i in range(3)]
    else:
        axes = [fig.add_subplot(1, 3, i+1) for i in range(3)]

    ax1, ax2, ax3 = axes

    has_legend = len(qspecs_dict) > 1 or (None not in qspecs_dict and len(qspecs_dict) > 0)
    
    first_qspec = None

    for label, qspec in qspecs_dict.items():
        if pixel is not None and len(qspec.shape) >= 2:
            q_to_plot = qspec.extract(pixel[0], pixel[1])
        else:
            q_to_plot = qspec

        if first_qspec is None:
            first_qspec = q_to_plot

        x = q_to_plot.qcenter
        y = q_to_plot.qhist
        
        # If still multi-dimensional (e.g. didn't provide pixel but it's a grid), 
        # we might need to handle it. For now assume it's 1D after extract or by input.
        if y.ndim > 1:
            y = y.reshape(-1, y.shape[-1])[0] # Take first one as fallback

        # Plotting parameters
        plot_params = kwargs.copy()
        if label is not None:
            plot_params['label'] = label

        ax1.plot(x, y, **plot_params)
        ax2.semilogy(x, y, **plot_params)
        ax3.loglog(x, y, **plot_params)

    if first_qspec is not None:
        try:
            plinlo, plinhi = first_qspec.quantiles([0.0, 0.95])
            ploglo, ploghi = first_qspec.quantiles([0.0, 0.98])
            
            ax1.set_xlim(plinlo, plinhi)
            ax2.set_xlim(ploglo, ploghi)
        except (ValueError, IndexError):
            pass

    ax1.set_ylabel('Probability density [1/DC]')
    ax1.grid(True)

    ax2.set_xlabel('ADC value [DC]')
    ax2.grid(True)
    # y-axis is already log from semilogy

    # Log-log plot: full positive range (1 to 4096 as per prototype)
    # ax3.set_xlim(1, max(4096, np.max(x)))
    ax3.grid(True)
    
    if has_legend:
        ax3.legend(loc=1)

    if title is not None:
        if pixel is not None:
            title += f": pixel[{pixel[0]},{pixel[1]}]"
        fig.suptitle(title)

    fig.tight_layout()

    return fig, axes

def calculate_max123(camera_images):
    """
    Calculates the number of times each pixel was the 1st, 2nd, or 3rd largest 
    in the camera for each event.
    
    Args:
        camera_images (CameraImages): The images and metadata container.
        
    Returns:
        tuple: (max1, max2, max3) as 32x32 numpy arrays of counts.
    """
    # images shape: (32, 32, n_events)
    images = camera_images.images
    n_events = images.shape[2]
    
    if n_events == 0:
        return np.zeros((32, 32), dtype=int), np.zeros((32, 32), dtype=int), np.zeros((32, 32), dtype=int)
    
    # Flatten pixels for each event: (1024, n_events)
    flat_images = images.reshape(1024, -1)
    
    # Find indices of 3 largest values for each event.
    # argsort is ascending, so we take the last 3.
    # top3_indices shape: (3, n_events)
    top3_indices = np.argsort(flat_images, axis=0)[-3:]
    
    # max1 = 1st largest, max2 = 2nd largest, max3 = 3rd largest
    max1_idx = top3_indices[2, :]
    max2_idx = top3_indices[1, :]
    max3_idx = top3_indices[0, :]
    
    # Count occurrences using bincount
    max1 = np.bincount(max1_idx, minlength=1024).reshape(32, 32)
    max2 = np.bincount(max2_idx, minlength=1024).reshape(32, 32)
    max3 = np.bincount(max3_idx, minlength=1024).reshape(32, 32)
    
    return max1, max2, max3
