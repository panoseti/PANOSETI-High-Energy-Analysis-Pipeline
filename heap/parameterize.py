"""parameterize

Functions for calculating the Hillas parameters of an image
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Ellipse

TAB10_COLORS = plt.get_cmap("tab10").colors

# 32 pixels over +-4.95 deg
CAMERA_HALF_WIDTH = 4.95 # deg
PLATE_SCALE = 2 * CAMERA_HALF_WIDTH / 32 # deg/pixel
CAMERA_EXTENT = [-CAMERA_HALF_WIDTH, CAMERA_HALF_WIDTH, -CAMERA_HALF_WIDTH, CAMERA_HALF_WIDTH] # imshow(image, origin="lower", extent=CAMERA_EXTENT) plots in degrees

# viridis from blue up, ramping from black, so empty camera pixels (0, or NaN after cleaning) read as black
CAMERA_CMAP = LinearSegmentedColormap.from_list("black_viridis", ["black", *plt.get_cmap("viridis")(np.linspace(0.3, 1, 8))])
CAMERA_CMAP.set_bad("black")

def calc_params(
        image: np.ndarray, 
        x: float=None, 
        y: float=None
    ):
    """
    Calculate the Hillas parameters

    Parameters:
        image: 2D camera image. numpy array with shape (32, 32), indexed [row, col]
        x, y: test position (degrees from camera center) from which to calculate e.g. distance. Default value is camera center.
    Returns dict with:
        - N_pix: the total number of pixels in the shower image
        - size: total intensity
        - x_c, y_c: centroid coordinates (degrees from camera center; x along columns, y along rows)
        - s_xx, s_yy, s_xy: second central moments (degrees^2)
        - length, width: rms major/minor axis (degrees)
        - miss: perpendicular distance to major axis (degrees)
        - distance: distance to test position x,y (degrees)
        - alpha: angle between image axis and distance (degrees)
        - phi: orientation angle (degrees), CCW from +x
    """

    image = np.asarray(image, dtype=float)
    N_pix = int(np.sum(np.isfinite(image)))
    size = float(np.nansum(image))

    if N_pix == 0:
        return {
            "N_pix": 0,
            "size": 0.0,
            "x_c": float("nan"),
            "y_c": float("nan"),
            "s_xx": float("nan"),
            "s_yy": float("nan"),
            "s_xy": float("nan"),
            "length": float("nan"),
            "width": float("nan"),
            "miss": float("nan"),
            "distance": float("nan"),
            "alpha": float("nan"),
            "phi": float("nan"),
        }

    H, W = image.shape
    # Set default x, y to image center if not provided
    if x is None:
        x = 0.0
    if y is None:
        y = 0.0
    # pixel centers in degrees from the camera center
    cols = (np.arange(W) - (W - 1) / 2.0) * PLATE_SCALE
    rows = (np.arange(H) - (H - 1) / 2.0) * PLATE_SCALE
    col_sums = np.nansum(image, axis=0)
    row_sums = np.nansum(image, axis=1)

    # centroid (degrees, origin at camera center)
    x_c = float(np.sum(col_sums * cols) / size)
    y_c = float(np.sum(row_sums * rows) / size)

    # second central moments
    dx = (cols - x_c)[None, :].astype(float)   
    dx = np.repeat(dx, H, axis=0)              
    dy = (rows - y_c)[:, None].astype(float)   
    dy = np.repeat(dy, W, axis=1) 

    weights = np.nan_to_num(image, nan=0.0)

    s_xx = float(np.sum(weights * dx * dx) / size)
    s_yy = float(np.sum(weights * dy * dy) / size)
    s_xy = float(np.sum(weights * dx * dy) / size)

    # length and width
    d = float(s_yy - s_xx)
    z = float(np.sqrt(d*d + 4*s_xy*s_xy))

    length = float(np.sqrt(max((s_xx + s_yy + z) / 2, 0.0)))
    width = float(np.sqrt(max((s_xx + s_yy - z) / 2, 0.0)))

    # orientation (phi)
    ac = float((d+z)*(y_c-y) + 2.0*s_xy*(x_c-x))
    bc = float(2.0*s_xy*(y_c-y) - (d-z)*(x_c-x))
    cc = float(np.sqrt(ac*ac + bc*bc))
    if cc == 0.0:
        # undefined orientation — choose sensible defaults
        cosphi = 1.0
        sinphi = 0.0
        phi = 0.0
    else:
        cosphi = float(bc / cc)
        sinphi = float(ac / cc)
        phi = float(np.arctan2(ac, bc))
        phi = float(phi % (2.0 * np.pi))

    # distance: distance from centroid to position x,y
    distance = float(np.hypot(x_c-x,y_c-y))

    # miss: perpendicular distance from point to major axis
    miss = float(abs(-sinphi*(x_c-x) + cosphi*(y_c-y)))
    if miss > distance:
        miss = distance
    
    # alpha: angle between image axis and test position x,y
    if distance == 0:
        alpha = float("nan")
    else:
        alpha = float(abs(np.arcsin(miss/distance)))

    return {
        "N_pix": N_pix,
        "size": size,
        "x_c": x_c,
        "y_c": y_c,
        "s_xx": s_xx,
        "s_yy": s_yy,
        "s_xy": s_xy,
        "length": length,
        "width": width,
        "miss": miss,
        "distance": distance,
        "alpha": float(np.rad2deg(alpha)),
        "phi": float(np.rad2deg(phi)),
    }


def draw_params(fig, params: dict, color: str="w", lw: float=2):
    """
    Draw the Hillas ellipse for a set of parameters on top of an existing figure.

    Parameters:
        fig: matplotlib Figure containing the image in degrees (e.g. from plt.imshow(image, origin="lower", extent=CAMERA_EXTENT))
        params: dict as returned by calc_params
        color: ellipse edge color
        lw: ellipse line width
    Returns the Axes the ellipse was drawn on.
    """
    ax = fig.gca()
    if params["N_pix"] == 0:
        return ax

    ellipse = Ellipse(
        (params["x_c"], params["y_c"]),
        2 * params["length"],
        2 * params["width"],
        angle=params["phi"],
        facecolor="none",
        edgecolor=color,
        lw=lw,
    )
    ax.add_patch(ellipse)
    return ax


def plot_hillas_histograms(dfs, title="Hillas Params", colors=None, pooled_df=None, pooled_label="All telescopes"):
    """Plots length/width/log10(size)/distance histograms (params_df already in degrees), one
    step-histogram line per entry in dfs overlaid on the same 4 axes. If pooled_df is given, also
    overlays it as a black alpha=0.4 filled histogram.

    dfs: {label: params_df}, e.g. one entry per telescope.
    colors: optional {label: color} so a label keeps the same color across calls; falls back to
        tab10 by position for any label not present.
    pooled_df: optional params_df pooled across labels, drawn filled alongside dfs' step lines.
    pooled_label: legend label for pooled_df (default = "All telescopes")

    Returns the Figure.
    """
    fig, axs = plt.subplots(2, 2, figsize=(12, 12))
    axs = axs.flatten()
    fig.suptitle(title)

    axs[0].set_title("length")
    axs[0].set_xlabel("degrees")
    axs[0].set_ylabel("normalized counts")

    axs[1].set_title("width")
    axs[1].set_xlabel("degrees")
    axs[1].set_ylabel("normalized counts")

    axs[2].set_title("log10(size)")
    axs[2].set_yscale("log")
    axs[2].set_xlabel("log10(ADC)")
    axs[2].set_ylabel("normalized counts")

    axs[3].set_title("distance")
    axs[3].set_xlabel("degrees")
    axs[3].set_ylabel("normalized counts")

    colors = colors or {}

    for i, (name, df) in enumerate(dfs.items()):
        color = colors.get(name, TAB10_COLORS[i % len(TAB10_COLORS)])
        label = f"{name} N={len(df)}"
        width_mean = df.width.mean()
        width_label = f"{label}, $\\mu$={width_mean:.3f}°"
        axs[0].hist(df.length, bins=80, range=(0, 2), histtype="step", density=True, label=label, color=color)
        axs[1].hist(df.width, bins=80, range=(0, 1), histtype="step", density=True, label=width_label, color=color)
        axs[1].axvline(width_mean, color=color, linestyle="--", linewidth=1.5)
        axs[2].hist(np.log10(df["size"]), bins=80, range=(0, 6), histtype="step", density=True, label=label, color=color)
        axs[3].hist(df.distance, bins=80, range=(0, 7), histtype="step", density=True, label=label, color=color)

    if pooled_df is not None:
        label = f"{pooled_label} N={len(pooled_df)}"
        width_mean = pooled_df.width.mean()
        width_label = f"{label}, $\\mu$={width_mean:.3f}°"
        axs[0].hist(pooled_df.length, bins=80, range=(0, 2), histtype="stepfilled", density=True, label=label, color="black", alpha=0.4)
        axs[1].hist(pooled_df.width, bins=80, range=(0, 1), histtype="stepfilled", density=True, label=width_label, color="black", alpha=0.4)
        axs[2].hist(np.log10(pooled_df["size"]), bins=80, range=(0, 6), histtype="stepfilled", density=True, label=label, color="black", alpha=0.4)
        axs[3].hist(pooled_df.distance, bins=80, range=(0, 7), histtype="stepfilled", density=True, label=label, color="black", alpha=0.4)

    for ax in axs:
        ax.legend(loc="upper right")

    fig.tight_layout()
    return fig
