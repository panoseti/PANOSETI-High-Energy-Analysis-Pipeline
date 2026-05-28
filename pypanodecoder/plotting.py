#!/usr/bin/env python3

# plotting.py: Standard plotting functions for PANOSETI data.

# Author: Stephen Fegan <sfegan@llr.in2p3.fr> (2026-05-28)
# Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

import numpy as np
import matplotlib.pyplot as plt

def plot_image(image, transpose=True, ax=None, fig=None, colorbar_label=None, show_colorbar=True, cmap='viridis', **kwargs):
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
            fig = plt.figure(dpi=200)
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

    pc = ax.imshow(image, **kwargs)

    # Delineate individual components with white lines
    if shape == (32, 32):
        # 32x32 pixels: 4x4 SiPMs of 8x8 pixels each.
        # Lines at boundaries 7.5, 15.5, 23.5 spanning from -0.5 to 31.5.
        ax.hlines(np.arange(7.5, 32, 8), -0.5, 31.5, color='w', lw=0.5)
        ax.vlines(np.arange(7.5, 32, 8), -0.5, 31.5, color='w', lw=0.5)
    elif shape == (4, 4):
        # 4x4 SiPMs: 2x2 Quabos of 2x2 SiPMs each.
        # Lines at boundary 1.5 spanning from -0.5 to 3.5.
        ax.hlines([1.5], -0.5, 3.5, color='w', lw=0.5)
        ax.vlines([1.5], -0.5, 3.5, color='w', lw=0.5)
    elif shape == (2, 2):
        # 2x2 Quabos: no internal divisions.
        pass

    if show_colorbar or colorbar_label is not None:
        cbar = fig.colorbar(pc, ax=ax)
        if colorbar_label is not None:
            cbar.set_label(colorbar_label)

    return fig, ax, pc
