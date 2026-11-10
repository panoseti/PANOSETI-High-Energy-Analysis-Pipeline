"""rotated_sipm

Fixes bad data where a quabo's wrong firmware read two of its 8x8 SiPMs out rotated by 180 deg: the
SiPMs listed in ROTATED_SIPM_FILE are rotated back in the raw frames of those nights before processing
(see heap.process_dataset.process_dataset()). See analysis_tools/fix_rotated_sipm.ipynb.

SiPMs are named [row, col] in the 4x4 grid of them in the (32, 32) [row, col] image: pixel rows
8*row to 8*row + 7, columns 8*col to 8*col + 7.
"""
from pathlib import Path

import numpy as np
import yaml
from matplotlib.patches import Rectangle

from heap.parameterize import CAMERA_HALF_WIDTH, PLATE_SCALE

ROTATED_SIPM_FILE = Path(__file__).resolve().parent.parent / "configs" / "analysis" / "rotated_sipm.yaml"

# each quabo's 16x16 [row, col] quadrant of the image, from ph1024 frames missing one quabo's packet (that
# quadrant reads 0): quabo 0, 1 and 3 seen directly, quabo 2 by elimination
QUABO_QUADRANTS = {0: (0, 0), 1: (0, 1), 2: (1, 1), 3: (1, 0)}


def rotate_sipm(frames, sipms):
    """
    A copy of frames with each 8x8 SiPM in sipms rotated by 180 deg. Doing it twice gives back the
    frames.

    Parameters:
        frames: images, shape (..., 32, 32) or (..., 1024), indexed [row, col]
        sipms: [[row, col], ...] SiPMs to rotate, each in the 4x4 grid of them

    Returns:
        the rotated frames, with the same shape as frames
    """
    frames = np.asarray(frames)
    images = frames.reshape(*frames.shape[:-1], 32, 32) if frames.shape[-1] == 1024 else frames
    images = images.copy()
    for row, col in sipms:
        rows, cols = slice(8*row, 8*row + 8), slice(8*col, 8*col + 8)
        images[..., rows, cols] = images[..., rows, cols][..., ::-1, ::-1]
    return images.reshape(frames.shape)


def load_rotated_sipm(path=ROTATED_SIPM_FILE):
    """
    ROTATED_SIPM_FILE: {telescope: [{"first": YYYYMMDD, "last": YYYYMMDD, "sipm": [[row, col], ...]}, ...]},
    the SiPMs to rotate back on every night from first to last (inclusive). {} if the file is missing or empty.
    """
    path = Path(path)
    if not path.exists():
        return {}
    with open(path) as f:
        return yaml.safe_load(f) or {}


def rotated_sipm_for(telescope, date, path=ROTATED_SIPM_FILE):
    """the [[row, col], ...] SiPMs of telescope to rotate back on night date (YYYYMMDD), sorted; [] for none"""
    sipms = set()
    for entry in load_rotated_sipm(path).get(telescope) or []:
        if str(entry["first"]) <= str(date) <= str(entry["last"]):
            sipms |= {tuple(sipm) for sipm in entry["sipm"]}
    return [list(sipm) for sipm in sorted(sipms)]


def quabo_ip(module, quabo):
    """quabo's IP address in a telescope's module (e.g. "module_253" or 253): the module's address (4 * module number) plus quabo"""
    number = 4 * int(str(module).removeprefix("module_")) + quabo
    return f"192.168.{number // 256}.{number % 256}"


def plot_layout(ax, module=None, highlight=(), highlight_label="rotated"):
    """
    The camera's quabos and SiPMs, in degrees as the images are plotted (imshow(origin="lower",
    extent=CAMERA_EXTENT)), with the pixel rows and columns on the top and right axes.

    Parameters:
        ax: matplotlib Axes
        module: optional, the telescope's module (e.g. "module_253"), to label each quabo with its IP address
        highlight: [[row, col], ...] SiPMs to fill in red
        highlight_label: their legend entry
    """
    half, width = CAMERA_HALF_WIDTH, 8 * PLATE_SCALE

    def corner(index):
        return -half + index * width

    for row, col in highlight:
        ax.add_patch(Rectangle((corner(col), corner(row)), width, width, color="C3", alpha=0.35, lw=0,
                               label=highlight_label if [row, col] == list(highlight[0]) else None))
    for k in range(5):
        lw, color = (2.5, "k") if k % 2 == 0 else (0.8, "0.5")
        ax.plot([corner(k)] * 2, [-half, half], color=color, lw=lw)
        ax.plot([-half, half], [corner(k)] * 2, color=color, lw=lw)
    for row in range(4):
        for col in range(4):
            ax.text(corner(col) + width / 2, corner(row) + width / 2, f"SiPM\n[{row}, {col}]",
                    ha="center", va="center", fontsize=8, color="0.3")
    for quabo, (qrow, qcol) in QUABO_QUADRANTS.items():
        label = f"quabo {quabo}" + (f"\n{quabo_ip(module, quabo)}" if module is not None else "")
        ax.text(corner(2*qcol + 1), corner(2*qrow + 1), label, ha="center", va="center", fontsize=9, fontweight="bold",
                bbox={"boxstyle": "round", "facecolor": "white", "edgecolor": "k"})

    ax.set(xlim=[-half, half], ylim=[-half, half], aspect="equal", xlabel="X (degrees)", ylabel="Y (degrees)")
    to_pixel = (lambda deg: (deg + half) / PLATE_SCALE, lambda pixel: pixel * PLATE_SCALE - half)
    top, right = ax.secondary_xaxis("top", functions=to_pixel), ax.secondary_yaxis("right", functions=to_pixel)
    top.set_ticks([0, 8, 16, 24, 32])
    right.set_ticks([0, 8, 16, 24, 32])
    top.set_xlabel("pixel column")
    right.set_ylabel("pixel row")
    if len(highlight):
        ax.legend(loc="upper left", bbox_to_anchor=(1.18, 1), fontsize=8)
