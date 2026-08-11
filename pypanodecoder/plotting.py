#!/usr/bin/env python3

# plotting.py: Standard plotting functions for PANOSETI data.

# Author: Stephen Fegan <sfegan@llr.in2p3.fr> (2026-05-28)
# Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

import numpy as np
import matplotlib.pyplot as plt
import math

def plot_image(image, transpose=False, ax=None, fig=None, colorbar_label=None, show_colorbar=True, cmap='viridis', plate_scale=None, **kwargs):
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

def plot_star_field(source=None, ra=None, dec=None, label=None, radius_deg=8.0, max_magnitude=5.0, ax=None, east_on_left=True, theta=0.0, **kwargs):
    """
    Plots a star field around a specific source or set of coordinates.

    Args:
        source (str, optional): Name of the source to look up.
        ra (float, optional): RA of the center (degrees).
        dec (float, optional): Dec of the center (degrees).
        label (str, optional): Label for the center source. If None, no cross or label is plotted.
        radius_deg (float): Radius to fetch stars for (default 8.0).
        max_magnitude (float): Maximum magnitude of stars to display (default 6.0).
        ax (matplotlib.axes.Axes, optional): Axes to plot on.
        east_on_left (bool): If True, East is to the left (standard astronomical view).
        theta (float): Rotation angle in degrees (default 0.0).
        **kwargs: Additional arguments passed to overlay_stars().

    Returns:
        PointingSolution: The pointing solution used for the plot.
    """
    from . import pointing
    from . import stars

    if source is not None:
        coords = pointing.get_target_coordinates(source)
        if coords is None:
            raise ValueError(f"Unknown source: {source}")
        ra0, dec0 = coords
    elif ra is not None and dec is not None:
        ra0, dec0 = ra, dec
    else:
        raise ValueError("Must provide either 'source' or both 'ra' and 'dec'.")

    if ax is None:
        ax = plt.gca()

    # Create a pointing solution centered at (ra0, dec0) at image (0,0)
    # We use a unit plate scale so that image coordinates are degrees from center
    ps = pointing.PointingSolution(ra0, dec0, plate_scale=1.0, theta=theta, east_on_left=east_on_left)

    # Fetch stars
    star_list = stars.get_bright_stars(ra0, dec0, radius_deg=1.5*radius_deg, max_magnitude=max_magnitude)

    ax.set_aspect('equal')
    # Set limits slightly larger than the visual aids
    limit = radius_deg if radius_deg > 5.5 else 6.0
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)

    # Plot center marker and label if requested
    if label is not None:
        ax.plot(0, 0, 'kx')
        ax.annotate(label, xy=(0, 0), xytext=(2.5, 2.5), textcoords='offset points',
                    ha='left', va='bottom', fontsize=8)

    # Default overlay options
    overlay_kwargs = {
        'color': 'b',
        'show_mags': True,
        'clip_to_axes': True,
        'fontsize': 4,
        'auto_align': True,
        'east_on_left': east_on_left,
        'clip_to_axes': True
    }
    overlay_kwargs.update(kwargs)

    # Plot stars
    overlay_stars(star_list, ps=ps, ax=ax, **overlay_kwargs)

    # Visual aids: 5-degree circle and square
    theta = np.linspace(0, 2 * np.pi, 360)
    ax.plot(5.0 * np.cos(theta), 5.0 * np.sin(theta), 'k--', lw=0.5, alpha=0.5)

    # Square boundary at +/- 5 degrees
    ax.vlines([-5.0, 5.0], -5.0, 5.0, colors='k', linestyles='--', lw=0.75, alpha=0.5)
    ax.hlines([-5.0, 5.0], -5.0, 5.0, colors='k', linestyles='--', lw=0.75, alpha=0.5)

    ax.grid(True, alpha=0.3)
    ax.set_xlabel('West [deg]' if east_on_left else 'East [deg]')
    ax.set_ylabel('North [deg]')

    return ps

def plot_north_west_guide(ps, ax=None, loc=1, color='r', zorder=None, **kwargs):
    """
    Draws a North/West L-guide in a specified corner of the plot.

    Args:
        ps (PointingSolution): The pointing solution for the field.
        ax (matplotlib.axes.Axes, optional): Axes to plot on.
        loc (int): Corner to display guide: 1 (top right), 2 (top left),
                   3 (bottom left), 4 (bottom right).
        color (str): Color for the guide lines and labels (default: 'r').
        zorder (float, optional): Layer ordering for the guide elements.
        **kwargs: Additional arguments passed to ax.text() for label styling (e.g., fontsize).
    """
    if loc not in [1, 2, 3, 4]:
        return

    if ax is None:
        ax = plt.gca()

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    # Handle potentially inverted axes for placement
    x_left, x_right = (xlim[1], xlim[0]) if ax.xaxis_inverted() else (xlim[0], xlim[1])
    y_bottom, y_top = (ylim[1], ylim[0]) if ax.yaxis_inverted() else (ylim[0], ylim[1])

    dx_span = x_right - x_left
    dy_span = y_top - y_bottom
    scale_ref = min(dx_span, dy_span)

    arm_len = 0.08 * scale_ref
    margin_x = 0.10 * dx_span
    margin_y = 0.10 * dy_span

    # Center of the imaginary square of the L
    if loc == 1:    # top right
        xc = x_right - margin_x
        yc = y_top - margin_y
    elif loc == 2:  # top left
        xc = x_left + margin_x
        yc = y_top - margin_y
    elif loc == 3:  # bottom left
        xc = x_left + margin_x
        yc = y_bottom + margin_y
    elif loc == 4:  # bottom right
        xc = x_right - margin_x
        yc = y_bottom + margin_y

    theta_rad = ps.theta
    cos_th = math.cos(theta_rad)
    sin_th = math.sin(theta_rad)
    f = -1.0 if ps.east_on_left else 1.0

    u_N = (f * sin_th, cos_th)
    u_W = (-f * cos_th, sin_th)

    # L origin (intersection point of the two arms)
    x0 = xc - 0.5 * arm_len * (u_N[0] + u_W[0])
    y0 = yc - 0.5 * arm_len * (u_N[1] + u_W[1])

    # Draw L guide lines
    plot_kwargs = {'color': color, 'lw': 1.5}
    if zorder is not None:
        plot_kwargs['zorder'] = zorder
    ax.plot([x0, x0 + arm_len * u_N[0]], [y0, y0 + arm_len * u_N[1]], **plot_kwargs)
    ax.plot([x0, x0 + arm_len * u_W[0]], [y0, y0 + arm_len * u_W[1]], **plot_kwargs)

    # Draw N and W labels
    theta_deg = math.degrees(theta_rad)
    arm_rotation_deg = -f * theta_deg
    text_rot = (arm_rotation_deg + 45) % 90 - 45

    label_offset = 0.25 * arm_len
    x_N = x0 + (arm_len + label_offset) * u_N[0]
    y_N = y0 + (arm_len + label_offset) * u_N[1]
    x_W = x0 + (arm_len + label_offset) * u_W[0]
    y_W = y0 + (arm_len + label_offset) * u_W[1]

    text_kwargs = {'ha': 'center', 'va': 'center', 'fontsize': 8}
    text_kwargs.update(kwargs)
    if zorder is not None:
        text_kwargs['zorder'] = zorder

    ax.text(x_N, y_N, 'N', color=color, rotation=text_rot, **text_kwargs)
    ax.text(x_W, y_W, 'W', color=color, rotation=text_rot, **text_kwargs)

def overlay_stars(stars, p1=None, p2=None, east_on_left=True, ax=None, use_index=False, clip_to_axes=False, color='r', plate_scale=None, ps=None, show_mags=False, show_cross=False, auto_align=False, loc=0, zorder=None, **kwargs):
    """
    Overlays stars on a field and optionally calculates the pointing solution.

    Args:
        stars (list): List of star dictionaries from stars.get_bright_stars().
        p1 (tuple, optional): (index, x, y) for the first reference star.
        p2 (tuple, optional): (index, x, y) for the second reference star.
        east_on_left (bool): If True, East is to the left (standard astronomical view).
        ax (matplotlib.axes.Axes, optional): Axes to plot on.
        use_index (bool): If True, label stars with their index in the stars list.
        clip_to_axes (bool): If True, only show stars within the current axes limits.
        color (str): Color for the star markers and labels (default: 'r').
        plate_scale (float, optional): If provided, force this plate scale (deg/pix).
        ps (PointingSolution, optional): Pre-calculated pointing solution.
        show_mags (bool): If True, show star magnitudes in parentheses after labels.
        auto_align (bool): If True, automatically set horizontal alignment based on
                           position to avoid clipping at edges.
        loc (int): Corner to display North-West guide (L): 0 to disable (default),
                   1 for top right, 2 for top left, 3 for bottom left, 4 for bottom right.
        **kwargs: Additional arguments passed to ax.text() for label styling (e.g., fontsize).

    Returns:
        PointingSolution: The pointing solution used.
    """
    from .pointing import PointingSolution

    if ps is None:
        if p1 is None or p2 is None:
            raise ValueError("Must provide either a PointingSolution (ps) or two reference stars (p1, p2).")
        i1, x1, y1 = p1
        i2, x2, y2 = p2
        ps = PointingSolution.solve(stars[i1], (x1, y1), stars[i2], (x2, y2), east_on_left=east_on_left, plate_scale=plate_scale)

    if ax is None:
        ax = plt.gca()

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    # Default text options
    text_kwargs = {'fontsize': 8, 'ha': 'left', 'va': 'bottom'}
    text_kwargs.update(kwargs)

    # Plot stars
    for i, star in enumerate(stars):
        x, y = ps.sky_to_image(star['ra_deg'], star['dec_deg'])

        if clip_to_axes:
            if not (min(xlim) <= x <= max(xlim) and min(ylim) <= y <= max(ylim)):
                continue

        # Use vmag for size
        size = max(5, (10 - star['vmag'])**2)
        scatter_kwargs = {'s': size, 'edgecolors': color, 'facecolors': 'none'}
        if zorder is not None:
            scatter_kwargs['zorder'] = zorder
        ax.scatter(x, y, **scatter_kwargs)
        if show_cross:
            cross_kwargs = {'marker': 'x', 's': 75, 'facecolors': color, 'linewidth': 0.3}
            if zorder is not None:
                cross_kwargs['zorder'] = zorder
            ax.scatter(x, y, **cross_kwargs)

        if use_index:
            label = str(i)
        else:
            name = star['name'] if star['name'] else f"HR{star['hr']}"
            # Merge multiple whitespaces
            bits = name.split()
            
            # If the first part has leading digits followed by letters, remove the digits
            if bits and bits[0]:
                stripped = bits[0].lstrip('0123456789')
                if stripped and stripped != bits[0]:
                    # Only replace if there are leftover characters (not purely digits)
                    bits[0] = stripped
            
            label = " ".join(bits)

        if show_mags:
            label = f"{label} ({star['vmag']:.1f})"

        # Copy kwargs to allow per-star modification
        current_text_kwargs = text_kwargs.copy()
        if auto_align:
            x_min, x_max = min(xlim), max(xlim)
            if (x - x_min) / (x_max - x_min) > 0.75:
                current_text_kwargs['ha'] = 'right'
            else:
                current_text_kwargs['ha'] = 'left'
            y_min, y_max = min(ylim), max(ylim)
            if (y - y_min) / (y_max - y_min) > 0.95:
                current_text_kwargs['va'] = 'top'
            else:
                current_text_kwargs['va'] = 'bottom'

        radius = math.sqrt(size) / 2.0 * 1.2
        offset_x = radius if current_text_kwargs.get('ha') == 'left' else -radius
        if zorder is not None:
            current_text_kwargs['zorder'] = zorder
        ax.annotate(label, xy=(x, y), xytext=(offset_x, 0),
                    textcoords='offset points',
                    color=color, **current_text_kwargs)

    # Draw North/West L-guide in specified corner
    if loc in [1, 2, 3, 4]:
        # Overwrite horizontal and vertical alignment for the guide
        guide_kwargs = text_kwargs.copy()
        guide_kwargs['ha'] = 'center'
        guide_kwargs['va'] = 'center'
        plot_north_west_guide(ps, ax=ax, loc=loc, color=color, zorder=zorder, **guide_kwargs)

    return ps
