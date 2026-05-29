#!/usr/bin/env python3

# plotting.py: Standard plotting functions for PANOSETI data.

# Author: Stephen Fegan <sfegan@llr.in2p3.fr> (2026-05-28)
# Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

import numpy as np
import matplotlib.pyplot as plt

def plot_image(image, transpose=True, ax=None, fig=None, colorbar_label=None, show_colorbar=True, cmap='viridis', plate_scale=None, **kwargs):
    """
    Plots a PANOSETI module image (pixels, SiPMs, or Quabos).
    
    The function supports three standard input shapes:
    - 32x32 array: Plots 1024 individual pixel values, delineating the 4x4 SiPMs
      (each 8x8 pixels) with white lines.
    - 4x4 array: Plots 16 SiPM values, delineating the 2x2 Quabos with white lines.
    - 2x2 array: Plots 4 Quabo values.
    
    Args:
        image (array_like): 2D array of values to plot. Must be of shape (32, 32), (4, 4), or (2, 2).
        transpose (bool): If True (default), transposes the image array before plotting.
        ax (matplotlib.axes.Axes, optional): Pre-existing axes to plot on.
        fig (matplotlib.figure.Figure, optional): Pre-existing figure.
        colorbar_label (str, optional): Label for the colorbar.
        show_colorbar (bool): If True (default), adds a colorbar to the plot.
        cmap (str): Colormap to use (default: 'viridis').
        plate_scale (float, optional): If provided, the axes will be scaled by this value
                                       and centered at (0, 0). For example, 1.0 gives
                                       centered pixel indices, while 0.31 could give
                                       degrees on the sky.
        **kwargs: Additional keyword arguments passed to ax.imshow().
        
    Returns:
        tuple: (matplotlib.figure.Figure, matplotlib.axes.Axes, matplotlib.image.AxesImage)
    """
    image = np.asarray(image)
    if image.ndim != 2:
        raise ValueError("Input image must be a 2-dimensional array.")
        
    shape = image.shape
    if shape not in [(32, 32), (4, 4), (2, 2)]:
        raise ValueError(
            f"Unsupported image shape {shape}. "
            f"Supported shapes are (32, 32) [pixels], (4, 4) [SiPMs], and (2, 2) [Quabos]."
        )
        
    if transpose:
        image = image.T
        
    if ax is None:
        if fig is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        else:
            ax = fig.gca()
    else:
        if fig is None:
            fig = ax.get_figure()

    if 'origin' not in kwargs:
        kwargs['origin'] = 'lower'
    if 'cmap' not in kwargs:
        kwargs['cmap'] = cmap

    # Handle coordinate scaling
    ny, nx = image.shape
    if plate_scale is not None:
        # Effective scale depends on the resolution of the image
        # plate_scale refers to a single pixel (1/32 of the module)
        if shape == (32, 32):
            eff_scale = plate_scale
        elif shape == (4, 4):
            eff_scale = 8 * plate_scale
        elif shape == (2, 2):
            eff_scale = 16 * plate_scale
        else:
            eff_scale = plate_scale

        # Calculate extent to center the image at (0,0)
        half_width = 0.5 * nx * eff_scale
        half_height = 0.5 * ny * eff_scale
        extent = [-half_width, half_width, -half_height, half_height]
        kwargs['extent'] = extent
        
        # Line positions in scaled coordinates
        x_min, x_max = -half_width, half_width
        y_min, y_max = -half_height, half_height
        
        # Step size for lines in scaled units
        if shape == (32, 32):
            # Lines every 8 pixels (SiPM boundaries)
            line_step = 8 * plate_scale
            line_start = -half_width + line_step
            lines = np.arange(line_start, half_width - 1e-9, line_step)
        elif shape == (4, 4):
            # Lines every 2 SiPMs (Quabo boundaries)
            line_step = 2 * eff_scale
            line_start = -half_width + line_step
            lines = np.arange(line_start, half_width - 1e-9, line_step)
        else:
            lines = []
    else:
        x_min, x_max = -0.5, nx - 0.5
        y_min, y_max = -0.5, ny - 0.5
        if shape == (32, 32):
            lines = np.arange(7.5, 32, 8)
        elif shape == (4, 4):
            lines = [1.5]
        else:
            lines = []

    pc = ax.imshow(image, **kwargs)

    # Delineate individual components with white lines
    if len(lines) > 0:
        ax.hlines(lines, x_min, x_max, color='w', lw=0.5)
        ax.vlines(lines, y_min, y_max, color='w', lw=0.5)

    # Set integer ticks for SiPM and Quabo images in index mode
    if plate_scale is None:
        if shape == (4, 4):
            ax.set_xticks(np.arange(4))
            ax.set_yticks(np.arange(4))
        elif shape == (2, 2):
            ax.set_xticks(np.arange(2))
            ax.set_yticks(np.arange(2))

    if show_colorbar or colorbar_label is not None:
        cbar = fig.colorbar(pc, ax=ax)
        if colorbar_label is not None:
            cbar.set_label(colorbar_label)

    return fig, ax, pc
