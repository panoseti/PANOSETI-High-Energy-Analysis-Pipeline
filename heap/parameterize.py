"""parameterize

Functions for calculating the Hillas parameters of an image
"""
import numpy as np
from matplotlib.patches import Ellipse

def calc_params(
        image: np.ndarray, 
        x: float=None, 
        y: float=None
    ):
    """
    Calculate the Hillas parameters

    Parameters:
        image: 2D camera image. numpy array with shape (32, 32)
        x, y: test position (pixels) from which to calculate e.g. distance. Default value is camera center.
    Returns dict with:
        - N_pix: the total number of pixels in the shower image
        - size: total intensity
        - x_c, y_c: centroid coordinates (pixels)
        - s_xx, s_yy, s_xy: second central moments 
        - length, width: rms major/minor axis (pixels)
        - miss: perpendicular distance to major axis (pixels)
        - distance: distance to test position x,y (pixels)
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
        x = (W - 1) / 2.0
    if y is None:
        y = (H - 1) / 2.0
    cols = np.arange(W)
    rows = np.arange(H)
    col_sums = np.nansum(image, axis=0)
    row_sums = np.nansum(image, axis=1)

    # centroid (pixel coordinates, origin at camera center)
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
        fig: matplotlib Figure containing the image (e.g. from plt.imshow)
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
