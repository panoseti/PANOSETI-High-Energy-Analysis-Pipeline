"""diagnostics

Tables and plots for checking heap.events.process_night()'s inputs (a night's run folders) and
outputs (<output_dir>/<date>/<telescope>/<source_slug>/{<source_slug>.npz, calibrations.npz}), for
analysis_tools/process_pff.ipynb and analysis_tools/inspect_data_products.ipynb.
"""
import json
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from heap.events import params_path
from heap.parameterize import CAMERA_CMAP, CAMERA_EXTENT, CAMERA_HALF_WIDTH, TAB10_COLORS, draw_params
from heap.process_dataset import _get_mount_hk, _mount_pointing, _run_start_epoch, discover_runs, identify_flip_side, identify_source, load_fallback_map


def run_table(raw_dir, telescopes: dict, data_product: str = "dp_ph1024"):
    """
    One night's runs as heap.process_dataset.process_dataset() sees them: each telescope's source,
    flip side and mount pointing (from hk.pff), and its data product files. Uses
    <raw_dir>/source_run_map.json as the fallback map if present (see
    heap.process_dataset.load_fallback_map()).

    Parameters:
        raw_dir: this night's raw data dir
        telescopes: {module: {"name": telescope_name, ...}}, as in the config's telescopes section
        data_product: data product tag in .pff filenames

    Returns:
        DataFrame with one row per run and telescope: Run, Start (UTC), Telescope, Source,
        FlipSide, RA, Dec (deg, median while tracking), Files, MB. Source/FlipSide are the
        error message if they couldn't be identified.
    """
    fallback_map_path = Path(raw_dir) / "source_run_map.json"
    fallback_map = load_fallback_map(fallback_map_path) if fallback_map_path.exists() else None
    rows = []
    for run_dir in discover_runs(raw_dir):
        for module, info in telescopes.items():
            name = info["name"]
            files = sorted(run_dir.glob(f"*{data_product}*{module}*.pff"))
            try:
                source = identify_source(run_dir, name, fallback_map=fallback_map)
            except ValueError as e:
                source = f"? ({e})"
            try:
                flip_side = identify_flip_side(run_dir, name, fallback_map=fallback_map)
            except ValueError as e:
                flip_side = f"? ({e})"
            mount = _get_mount_hk(run_dir, name)
            pointing = _mount_pointing(mount) if mount is not None else None
            rows.append({
                "Run": run_dir.name,
                "Start": pd.to_datetime(_run_start_epoch(run_dir), unit="s", utc=True),
                "Telescope": name,
                "Source": source,
                "FlipSide": flip_side,
                "RA": round(pointing.ra.deg, 3) if pointing is not None else np.nan,
                "Dec": round(pointing.dec.deg, 3) if pointing is not None else np.nan,
                "Files": len(files),
                "MB": round(sum(f.stat().st_size for f in files) / 1e6, 1),
            })
    return pd.DataFrame(rows)


def product_table(output_dir, dates: list = None):
    """
    Every telescope and source process_night() wrote under output_dir.

    Parameters:
        output_dir: processed data dir (holding <date>/)
        dates: nights to list, e.g. ["20260917"]; None = every night in output_dir

    Returns:
        DataFrame with one row per night, telescope and source: Date, Telescope, Source (its folder
        name, e.g. NGC_1275), Frames
        (after the packet-loss and spike cuts), Start, End (UTC), Hours, Rate (Hz, Frames over
        End - Start), Cleaned (fraction of frames with any pixel left after cleaning), GainSource,
        MB (both npz files), Settings (processing.json). A telescope that was processed but
        wrote no source (no usable frames) gets one row with Source None.
    """
    output_dir = Path(output_dir)
    dates = dates or sorted(p.name for p in output_dir.iterdir() if p.is_dir() and p.name.isdigit())
    rows = []
    for date in dates:
        if not (output_dir / date).is_dir():
            continue
        for tel_dir in sorted(p for p in (output_dir / date).iterdir() if p.is_dir()):
            settings_path = tel_dir / "processing.json"
            settings = json.loads(settings_path.read_text()) if settings_path.exists() else None
            source_dirs = sorted(p for p in tel_dir.iterdir() if (p / "calibrations.npz").exists())
            if not source_dirs:
                rows.append({"Date": date, "Telescope": tel_dir.name, "Source": None, "Frames": 0, "Settings": settings})
            for source_dir in source_dirs:
                npz = np.load(source_dir / f"{source_dir.name}.npz")
                calib = np.load(source_dir / "calibrations.npz")
                timestamps, n_pix = npz["Timestamp"], npz["N_pix"]
                start, end = timestamps.min(), timestamps.max()
                rows.append({
                    "Date": date,
                    "Telescope": tel_dir.name,
                    "Source": source_dir.name,
                    "Frames": len(timestamps),
                    "Start": pd.to_datetime(start, unit="s", utc=True).floor("s"),
                    "End": pd.to_datetime(end, unit="s", utc=True).floor("s"),
                    "Hours": round((end - start) / 3600, 2),
                    "Rate": round(len(timestamps) / (end - start), 2) if end > start else np.nan,
                    "Cleaned": round(np.mean(n_pix > 0), 3),
                    "GainSource": calib["gain_source"].item(),
                    "MB": round(sum(f.stat().st_size for f in source_dir.glob("*.npz")) / 1e6, 1),
                    "Settings": settings,
                })
    return pd.DataFrame(rows)


def load_products(output_dir, date: str, telescope: str, source: str):
    """
    One telescope's processed night of source, as written by heap.process_dataset.process_dataset().

    Returns:
        params: DataFrame of Hillas parameters (degrees, see heap.parameterize.calc_params()), Event and
            Timestamp, one row per frame, in the order frames were processed (preflip runs then
            postflip runs), not sorted by time
        images: (n, 32, 32) cleaned, gain-corrected images, same rows as params; NaN where
            cleaning removed the pixel
        calib: calibrations.npz (pedestals and pedestal_variances (n, 1024), same rows as params;
            gain (32, 32); gain_source; gain_caption)
    """
    path = params_path(output_dir, date, telescope, source)
    npz = np.load(path)
    params = pd.DataFrame({col: npz[col] for col in npz.files if col != "cleaned_images"})
    return params, npz["cleaned_images"], np.load(path.parent / "calibrations.npz")


# what each key of <source_slug>.npz and calibrations.npz holds; n = frames, same order in both files
PRODUCT_KEYS = {
    "cleaned_images": "(n, 32, 32) [frame, row, col]: (raw - pedestal) / gain, NaN where cleaning removed the pixel",
    "N_pix": "pixels left after cleaning (0 = nothing survived, Hillas parameters NaN)",
    "size": "sum of the cleaned image (ADC, gain-corrected)",
    "x_c": "centroid x, along the columns (deg from camera center)",
    "y_c": "centroid y, along the rows (deg from camera center)",
    "s_xx": "second central moment in x (deg^2)",
    "s_yy": "second central moment in y (deg^2)",
    "s_xy": "second central moment in xy (deg^2)",
    "length": "rms spread along the image axis (deg)",
    "width": "rms spread across the image axis (deg)",
    "miss": "distance from the camera center to the image axis (deg)",
    "distance": "distance from the camera center to the centroid (deg)",
    "alpha": "angle between the image axis and the line from the centroid to the camera center (deg)",
    "phi": "image axis orientation, counterclockwise from +x (deg, 0-360)",
    "Event": "row index, 0 to n-1",
    "Timestamp": "frame time (Unix seconds, UTC), not corrected between telescopes",
    "pedestals": "(n, 1024) [frame, row*32 + col]: pedestal (ADC), constant within each 600 s window, 0 if never set",
    "pedestal_variances": "(n, 1024) [frame, row*32 + col]: pedvar (ADC), the cleaning threshold unit, 0 if never set",
    "gain": "(32, 32) [row, col]: relative gain for the whole night, mean 1",
    "gain_source": "0-d string: which fallback made the gain map (own, alt, other_night:<date>, flat)",
    "gain_caption": "0-d string: the star fields the gain map came from",
}


def describe_npz(path):
    """
    Every array in an npz file process_dataset() writes (<source_slug>.npz or calibrations.npz):
    shape, dtype, range, NaN count and what it holds (see PRODUCT_KEYS).

    Returns:
        DataFrame indexed by key
    """
    rows = []
    with np.load(path) as npz:
        for key in npz.files:
            values = npz[key]
            row = {"shape": values.shape, "dtype": str(values.dtype)}
            if values.ndim == 0:
                row["value"] = values.item()
            elif np.issubdtype(values.dtype, np.number):
                finite = values[np.isfinite(values)]
                row["min"] = finite.min() if finite.size else np.nan
                row["max"] = finite.max() if finite.size else np.nan
                row["NaN"] = values.size - finite.size
            row["description"] = PRODUCT_KEYS.get(key, "")
            rows.append({"key": key, **row})
    df = pd.DataFrame(rows).set_index("key")
    return df[[col for col in ["shape", "dtype", "value", "min", "max", "NaN", "description"] if col in df]]


def plot_gain_maps(output_dir, dates: list, telescopes: list, source: str):
    """
    Each night's (rows) and telescope's (columns) gain map for source, titled with its gain_source
    (see heap.process_dataset.resolve_gain_map()), plus the other source and flip side for "alt",
    on a shared color scale.

    Returns:
        the Figure
    """
    fig, axs = plt.subplots(len(dates), len(telescopes), figsize=(3.2*len(telescopes), 3*len(dates)), squeeze=False)
    im = None
    for i, date in enumerate(dates):
        for j, telescope in enumerate(telescopes):
            ax = axs[i, j]
            path = params_path(output_dir, date, telescope, source).parent / "calibrations.npz"
            ax.set_xticks([])
            ax.set_yticks([])
            if not path.exists():
                ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)
            else:
                calib = np.load(path)
                im = ax.imshow(calib["gain"], origin="lower", vmin=0.5, vmax=1.6)
                title = calib["gain_source"].item()
                if title == "alt": # caption: "<Side>: <source> · <Side>: <other source> (different field)"
                    side, other = calib["gain_caption"].item().split(" · ")[1].removesuffix(" (different field)").split(": ", 1)
                    title = f"alt: {other} ({side.lower()})"
                ax.set_title(title, fontsize=9)
            if i == 0:
                ax.set_title(f"{telescope}\n{ax.get_title()}", fontsize=9)
            if j == 0:
                ax.set_ylabel(date)
    if im is not None:
        fig.colorbar(im, ax=axs, label="Relative gain", shrink=0.8)
    fig.suptitle(f"{source} gain maps")
    return fig


def _sorted_calibrations(params, calib):
    """calib's timestamps, pedestals and pedvars sorted by time, and the sorted row where each
    window's pedestal/pedvar (see heap.process_dataset.build_calibrations()) starts, leaving out
    frames whose pedestal was never set (all zero)."""
    order = np.argsort(params.Timestamp.to_numpy())
    timestamps = params.Timestamp.to_numpy()[order]
    pedestals = calib["pedestals"][order]
    pedvars = calib["pedestal_variances"][order]
    changed = np.r_[True, np.any(pedvars[1:] != pedvars[:-1], axis=1)]
    starts = np.flatnonzero(changed & pedvars.any(axis=1))
    return timestamps, pedestals, pedvars, starts


def plot_calibrations(params, calib, title=""):
    """
    One telescope's night of calibrations (see load_products()): the camera's median pedestal and
    pedvar over time (shaded: 16-84% of pixels), titled with how its gain map was made
    (gain_caption). Frames whose pedestal was never set (a window with too few frames and none
    before it, see heap.process_dataset.build_calibrations()) are left out.

    Returns:
        the Figure
    """
    timestamps, pedestals, pedvars, _ = _sorted_calibrations(params, calib)
    valid = pedvars.any(axis=1)
    time = pd.to_datetime(timestamps[valid], unit="s", utc=True)

    fig, axs = plt.subplots(1, 2, figsize=(14, 3.5))
    for ax, (values, label) in zip(axs, [(pedestals, "Pedestal"), (pedvars, "Pedvar")]):
        lo, mid, hi = np.percentile(values[valid], [16, 50, 84], axis=1)
        ax.fill_between(time, lo, hi, alpha=0.3, step="post")
        ax.step(time, mid, where="post")
        ax.set_ylabel(f"{label} (camera median)")
        ax.set_xlabel("Time (UTC)")
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    fig.suptitle(f"{title}\n{calib['gain_caption'].item()}", fontsize=10)
    fig.tight_layout()
    return fig


def plot_calibration_windows(params, calib, title="", n_cols: int = 6):
    """
    The pedestal and pedvar maps applied in each window they were recomputed in (every 600 s, see
    heap.process_dataset.build_calibrations()), one panel per window, each figure on one color scale
    (1-99% of all windows' pixels) so changes over the night stand out, e.g. a star raising the
    pedvar of the pixels it moves through.

    Parameters:
        params, calib: see load_products()
        title: prefix for each figure's title
        n_cols: panels per row

    Returns:
        (pedestal_fig, pedvar_fig), or (None, None) if no frame has a pedestal
    """
    timestamps, pedestals, pedvars, starts = _sorted_calibrations(params, calib)
    if not len(starts):
        return None, None

    n_rows = int(np.ceil(len(starts) / n_cols))
    figs = []
    for values, label in [(pedestals, "Pedestal"), (pedvars, "Pedvar")]:
        maps = values[starts].reshape(-1, 32, 32)
        vmin, vmax = np.percentile(maps, [1, 99])
        fig, axs = plt.subplots(n_rows, n_cols, figsize=(2.4*n_cols, 2.3*n_rows + 0.6), squeeze=False, layout="constrained")
        for ax, values_map, start in zip(axs.flat, maps, starts):
            im = ax.imshow(values_map, origin="lower", vmin=vmin, vmax=vmax)
            ax.set_title(f"{pd.to_datetime(timestamps[start], unit='s', utc=True):%H:%M:%S}", fontsize=8)
            ax.set_xticks([])
            ax.set_yticks([])
        for ax in axs.flat[len(starts):]:
            ax.axis("off")
        fig.colorbar(im, ax=axs, label=label, shrink=0.8)
        fig.suptitle(f"{title} {label.lower()} per window (start, UTC)".strip(), fontsize=10)
        figs.append(fig)
    return tuple(figs)


def plot_rates(output_dir, date: str, telescopes: list, source: str, min_npix: int = 3, bin_width: float = 60, colors: dict = None):
    """
    Each telescope's rate of frames kept after the packet-loss and spike cuts (top), and of images
    with at least min_npix pixels left after cleaning (bottom), over one night of source.

    Returns:
        the Figure
    """
    colors = colors or {}
    fig, axs = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    for i, telescope in enumerate(telescopes):
        path = params_path(output_dir, date, telescope, source)
        if not path.exists():
            continue
        npz = np.load(path)
        timestamps, n_pix = npz["Timestamp"], npz["N_pix"]
        bins = np.arange(timestamps.min(), timestamps.max() + bin_width, bin_width)
        time = pd.to_datetime(bins[:-1], unit="s", utc=True)
        color = colors.get(telescope, TAB10_COLORS[i % len(TAB10_COLORS)])
        for ax, mask in zip(axs, [np.ones_like(n_pix, dtype=bool), n_pix >= min_npix]):
            counts, _ = np.histogram(timestamps[mask], bins=bins)
            ax.step(time, counts / bin_width, where="post", color=color, label=telescope)
    axs[0].set_ylabel("Frames after cuts [Hz]")
    axs[1].set_ylabel(f"Images with N_pix >= {min_npix} [Hz]")
    axs[1].set_xlabel("Time (UTC)")
    axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[0].legend(loc="upper right")
    axs[0].set_title(f"{source} {date}, {bin_width:g} s bins")
    fig.tight_layout()
    return fig


def plot_image_grid(params, images, indices, n_cols: int = 4):
    """
    Cleaned images (see load_products()) with their Hillas ellipse (1 sigma length/width), in
    degrees from the camera center.

    Parameters:
        params, images: see load_products()
        indices: rows of params/images to show

    Returns:
        the Figure
    """
    n_rows = int(np.ceil(len(indices) / n_cols))
    fig, axs = plt.subplots(n_rows, n_cols, figsize=(3.5*n_cols, 3.2*n_rows), squeeze=False)
    for ax, idx in zip(axs.flat, indices):
        row = params.iloc[idx]
        im = ax.imshow(images[idx], origin="lower", extent=CAMERA_EXTENT, cmap=CAMERA_CMAP, vmin=0)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        fig.sca(ax)
        draw_params(fig, row, lw=1.5)
        ax.set_title(f"{idx}: N_pix={int(row.N_pix)}, size={row['size']:.0f}\n{pd.to_datetime(row.Timestamp, unit='s', utc=True):%H:%M:%S.%f} UTC", fontsize=8)
    for ax in axs.flat[len(indices):]:
        ax.axis("off")
    fig.tight_layout()
    return fig


def plot_centroids(images, title="", telescopes: list = None, bins: int = 32):
    """
    2D histogram of image centroids (x_c, y_c) in camera coordinates (degrees), one panel per
    telescope, binned by camera pixel by default.

    Parameters:
        images: rows with Telescope, x_c, y_c, e.g. from heap.events.build_array_events() or
            heap.events.apply_cuts()
        title: figure title
        telescopes: panel order (default = every telescope in images, sorted)
        bins: bins per axis across the camera (32 = one per pixel)

    Returns:
        the Figure
    """
    telescopes = telescopes or sorted(images.Telescope.unique())
    edges = np.linspace(-CAMERA_HALF_WIDTH, CAMERA_HALF_WIDTH, bins + 1)
    fig, axs = plt.subplots(1, len(telescopes), figsize=(4.2*len(telescopes), 4.2), squeeze=False, layout="constrained")
    for ax, telescope in zip(axs[0], telescopes):
        tel = images[images.Telescope == telescope].dropna(subset=["x_c", "y_c"])
        _, _, _, im = ax.hist2d(tel.x_c, tel.y_c, bins=[edges, edges], cmap=CAMERA_CMAP)
        fig.colorbar(im, ax=ax, shrink=0.8)
        ax.set(aspect="equal", xlabel="x_c (deg)", ylabel="y_c (deg)", title=f"{telescope} (N={len(tel)})")
    fig.suptitle(title)
    return fig
