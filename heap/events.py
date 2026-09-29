"""events

Builds array events (images from 2+ telescopes of the same shower) out of per-telescope Hillas
parameters (<output_dir>/<date>/<telescope>/<source_slug>/<source_slug>.npz, see process_night()
and heap.process_dataset.process_dataset()), applies the analysis cuts, and plots single events.

Hillas parameters are in degrees on the panodisplay_REALDATA.C camera (32 pixels over +-4.95 deg,
origin at camera center, see heap.parameterize.calc_params()). x_c runs along the columns and y_c
along the rows of the (32, 32) image, so x_c/y_c here are the transpose of the ROOT CSVs'
MeanX/MeanY (ROOT fills bin (i+1, j+1) from pixel [i][j]).
"""
import json
from functools import cache
from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from matplotlib.colors import to_rgba
from matplotlib.patches import Ellipse

from heap import coincidences as coinc
from heap.parameterize import CAMERA_CMAP, CAMERA_HALF_WIDTH, TAB10_COLORS, plot_hillas_histograms
from heap.significance import make_wcs
from heap.process_dataset import _get_mount_hk, _mount_pointing, _run_start_epoch, discover_runs, group_runs_by_source, load_fallback_map, process_dataset, slugify



def params_path(output_dir, date, telescope, source):
    """Path to one telescope's <source_slug>.npz for one night, see process_night()."""
    source_slug = slugify(source, sep="_")
    return Path(output_dir) / date / telescope / source_slug / f"{source_slug}.npz"


def process_night(
        raw_dir,
        output_dir,
        date: str,
        telescopes: dict,
        data_product: str = "dp_ph1024",
        image_threshold: float = 4.0,
        border_threshold: float = 2.0,
        keep_brightest_island: bool = False,
        reprocess: bool = False,
):
    """
    Run heap.process_dataset.process_dataset() for every telescope with data in one night, writing
    <output_dir>/<date>/<telescope>/<source_slug>/ for every source that night. Keeping every night
    under the same output_dir lets process_dataset() fall back to a prior night's gain map.

    Each telescope's settings (data_product, rate_cut, and the cleaning settings) are recorded in
    <output_dir>/<date>/<telescope>/processing.json; existing output made with other (or
    unrecorded) settings raises instead of being reused, unless reprocess.

    Parameters:
        raw_dir: this night's raw data dir, containing one subfolder per run
        output_dir: processed data dir (holding <date>/)
        date: night, YYYYMMDD
        telescopes: {module: {"name": telescope_name, "rate_cut": Hz}}, e.g. {"module_253": {"name": "Winter", "rate_cut": 100}}
        data_product: data product tag in .pff filenames
        image_threshold, border_threshold, keep_brightest_island: passed through to process_dataset()
        reprocess: rerun telescopes that already have output for this night
    """
    raw_dir = Path(raw_dir)
    fallback_map_path = raw_dir / "source_run_map.json"
    runs = discover_runs(raw_dir)

    for module, info in telescopes.items():
        name = info["name"]
        out_dir = Path(output_dir) / date / name
        settings = {
            "data_product": data_product, "rate_cut": info.get("rate_cut", 20),
            "image_threshold": image_threshold, "border_threshold": border_threshold,
            "keep_brightest_island": keep_brightest_island,
        }
        settings_path = out_dir / "processing.json"
        if not reprocess and any(out_dir.glob("*/calibrations.npz")):
            previous = json.loads(settings_path.read_text()) if settings_path.exists() else None
            if previous == settings:
                continue
            raise ValueError(
                f"{out_dir} was processed with {previous or 'unrecorded settings'}, not {settings}; "
                "reprocess, or use another output_dir"
            )
        module_pattern = f"{data_product}*{module}"
        if not any(any(run_dir.glob(f"*{module_pattern}*.pff")) for run_dir in runs):
            print(f"{date}: no {module_pattern} .pff files for {name}, skipping")
            continue
        try:
            process_dataset(
                raw_dir, out_dir, module_pattern, name,
                fallback_map_path=fallback_map_path if fallback_map_path.exists() else None,
                rate_cut=info.get("rate_cut", 20),
                image_threshold=image_threshold, border_threshold=border_threshold,
                keep_brightest_island=keep_brightest_island,
            )
            settings_path.write_text(json.dumps(settings, indent=2))
        except Exception as e:
            print(f"{date}: processing {name} failed: {e}")


def source_runs(raw_dir, telescope, source):
    """
    One night's runs of source, by flip side, from hk.pff (see
    heap.process_dataset.group_runs_by_source()). Uses <raw_dir>/source_run_map.json as the
    fallback map if present.

    Returns:
        {"preflip": [run_dir, ...], "postflip": [run_dir, ...]}
    """
    raw_dir = Path(raw_dir)
    fallback_map_path = raw_dir / "source_run_map.json"
    fallback_map = load_fallback_map(fallback_map_path) if fallback_map_path.exists() else None
    return group_runs_by_source(raw_dir, telescope, fallback_map=fallback_map).get(source, {"preflip": [], "postflip": []})


def _latest_run(timestamps, runs_by_flip):
    """Index into the start-sorted runs of runs_by_flip of the latest run started before each
    timestamp, and those runs as (start epoch, flip side, run_dir)."""
    starts = sorted(
        (_run_start_epoch(run_dir), flip_side, run_dir)
        for flip_side, run_dirs in runs_by_flip.items() for run_dir in run_dirs
    )
    if not starts:
        raise ValueError("No runs to take flip sides from")
    epochs = np.array([s[0] for s in starts])
    idx = np.searchsorted(epochs, timestamps, side="right") - 1
    return np.clip(idx, 0, None), starts


def flip_side_of(timestamps, runs_by_flip):
    """
    Flip side of each timestamp: that of the latest run in runs_by_flip started before it.

    Parameters:
        timestamps: Unix timestamps (s)
        runs_by_flip: {"preflip": [run_dir, ...], "postflip": [...]}, see source_runs()

    Returns:
        array of "preflip"/"postflip", same shape as timestamps
    """
    idx, starts = _latest_run(timestamps, runs_by_flip)
    return np.array([s[1] for s in starts])[idx]


def run_of(timestamps, runs_by_flip):
    """
    Run of each timestamp: the folder name of the latest run in runs_by_flip started before it.

    Parameters:
        timestamps: Unix timestamps (s)
        runs_by_flip: {"preflip": [run_dir, ...], "postflip": [...]}, see source_runs()

    Returns:
        array of run folder names, same shape as timestamps
    """
    idx, starts = _latest_run(timestamps, runs_by_flip)
    return np.array([Path(s[2]).name for s in starts])[idx]


def run_pointings(runs_by_flip, telescope):
    """
    Each run's median RA/Dec while tracking, from hk.pff's MOUNT_<TELESCOPE> table (see
    heap.process_dataset._mount_pointing()).

    Returns:
        {run folder name: SkyCoord}, leaving out runs with no MOUNT_<TELESCOPE> table or no tracking
    """
    pointings = {}
    for run_dirs in runs_by_flip.values():
        for run_dir in run_dirs:
            mount = _get_mount_hk(run_dir, telescope)
            pointing = _mount_pointing(mount) if mount is not None else None
            if pointing is not None:
                pointings[Path(run_dir).name] = pointing
    return pointings


def load_camera_frame(
        npz_path,
        telescope: str,
        runs_by_flip: dict,
        rotate_postflip: bool = True,
        rel_efficiency: float = 1.0,
        pointing_corrections: dict = None,
):
    """
    Load one telescope's Hillas parameters for one night into the camera frame (degrees), sorted by
    time.

    Parameters:
        npz_path: <source_slug>.npz written by heap.process_dataset.process_dataset()
        telescope: telescope name, stored in the Telescope column
        runs_by_flip: this night's runs of the source, see source_runs()
        rotate_postflip: rotate postflip images by 180 deg into the preflip camera frame
        rel_efficiency: relative telescope efficiency; size is divided by it
        pointing_corrections: optional {flip_side: (dx, dy)} in deg, subtracted from x_c/y_c

    Returns:
        DataFrame with ImageIndex (row in the npz's cleaned_images, which holds every run of the source that night), Telescope, Timestamp,
        FlipSide, and the npz's Hillas parameters (x_c, y_c, phi, size, N_pix, length, width, miss, distance, alpha)
    """
    npz = np.load(npz_path)
    p = pd.DataFrame({col: npz[col] for col in npz.files if col != "cleaned_images"})
    p = p.sort_values("Timestamp", ignore_index=True)

    df = pd.DataFrame({
        "ImageIndex": p.Event,
        "Telescope": telescope,
        "Timestamp": p.Timestamp,
        "FlipSide": flip_side_of(p.Timestamp.to_numpy(), runs_by_flip),
        "x_c": p.x_c,
        "y_c": p.y_c,
        "phi": p.phi,
        "size": p["size"] / rel_efficiency,
        "N_pix": p.N_pix,
        "length": p.length,
        "width": p.width,
        "miss": p.miss,
        "distance": p.distance,
        "alpha": p.alpha,
    })

    postflip = df.FlipSide == "postflip"
    if rotate_postflip:
        df.loc[postflip, ["x_c", "y_c"]] *= -1
        df.loc[postflip, "phi"] = (df.loc[postflip, "phi"] + 180) % 360

    for side, (dx, dy) in (pointing_corrections or {}).items():
        df.loc[df.FlipSide == side, "x_c"] -= dx
        df.loc[df.FlipSide == side, "y_c"] -= dy

    return df


def build_events(telescopes: dict, reference: str, window: float = 0.001, plotting: bool = False):
    """
    Match one night's images across telescopes into events.

    Each telescope's timestamps are corrected against reference (see coincidences.correct_time),
    then matched pairwise and merged (see coincidences.find_coincidences). Events where a
    telescope matched more than one image are dropped as ambiguous.

    Parameters:
        telescopes: {telescope_name: DataFrame from load_camera_frame()}
        reference: telescope other telescopes' timestamps get corrected against
        window: coincidence window (s)
        plotting: show correct_time's before/after plots

    Returns:
        DataFrame of every image in an event (telescopes' columns plus Event, numbered from 0),
        sorted by Event then Telescope, or None if reference has no data
    """
    telescopes = {name: df for name, df in telescopes.items() if df is not None and len(df)}
    if reference not in telescopes:
        return None

    ref_timestamps = telescopes[reference].Timestamp.to_numpy()
    for name, df in telescopes.items():
        if name == reference:
            continue
        df = df.copy()
        try:
            df["Timestamp"] = coinc.correct_time(df.Timestamp.to_numpy(), ref_timestamps, plot_name=f"{reference}-{name}", base_dir=None, plotting=plotting)
        except ValueError:
            print(f"No {reference}-{name} coincidences to correct timing with, leaving uncorrected")
        # timestamps in a correction bin with no coincidences come back nan
        telescopes[name] = df.dropna(subset=["Timestamp"]).sort_values("Timestamp", ignore_index=True)

    groups = coinc.find_coincidences({name: (df.Timestamp.to_numpy(), None, None) for name, df in telescopes.items()}, window=window)
    n_groups = len(groups)
    groups = [g for g in groups if len(set(g[0])) == len(g[0])]
    print(f"{len(groups)} events; dropped {n_groups - len(groups)} ambiguous groups (a telescope matched more than one image)")

    parts = []
    for name, df in telescopes.items():
        event, idx = [], []
        for e, (tel_names, event_idx, _) in enumerate(groups):
            if name in tel_names:
                event.append(e)
                idx.append(event_idx[tel_names.index(name)])
        part = df.iloc[idx].copy()
        part.insert(0, "Event", event)
        parts.append(part)

    return pd.concat(parts).sort_values(["Event", "Telescope"], ignore_index=True)


def build_array_events(
        output_dir,
        raw_dir,
        date: str,
        source: str,
        telescopes: list,
        reference: str,
        window: float = 0.001,
        rotate_postflip: bool = True,
        rel_tel_efficiency: dict = None,
        pointing_corrections: dict = None,
        plotting: bool = False,
):
    """
    load_camera_frame() for every telescope with processed data for source on date, then
    build_events(). Adds a Date column holding date and a Run column holding each event's run
    folder name (see run_of(), with reference's runs).

    Parameters:
        output_dir: processed data dir (holding <date>/), see process_night()
        raw_dir: this night's raw data dir, for flip sides (see source_runs())
        rel_tel_efficiency: optional {telescope: efficiency}, see load_camera_frame()
        pointing_corrections: optional {(date, telescope, flip_side): (dx, dy)} in deg
        window, reference, rotate_postflip, plotting: see build_events()/load_camera_frame()

    Returns:
        see build_events()
    """
    frames = {}
    runs = {}
    for name in telescopes:
        path = params_path(output_dir, date, name, source)
        if not path.exists():
            continue
        corrections = {
            side: offset for (d, tel, side), offset in (pointing_corrections or {}).items()
            if d == date and tel == name
        }
        runs[name] = source_runs(raw_dir, name, source)
        frames[name] = load_camera_frame(
            path, name, runs[name],
            rotate_postflip=rotate_postflip,
            rel_efficiency=(rel_tel_efficiency or {}).get(name, 1.0),
            pointing_corrections=corrections,
        )

    events = build_events(frames, reference, window=window, plotting=plotting)
    if events is not None:
        events.insert(1, "Date", date)
        # each event's run, from its earliest image and reference's runs, whose mount gives the pointing
        events.insert(2, "Run", run_of(events.groupby("Event").Timestamp.transform("min").to_numpy(), runs[reference]))
    return events


def apply_cut(array, column, mode, threshold):
    """
    Keep images with column below the cut.

    Parameters:
        array: images, with Telescope and column
        mode: "nsigma" (per-telescope mean + threshold*std) or "value" (threshold, same for all telescopes)
        threshold: nsigma, or the cut value
    """
    if mode == "nsigma":
        cut_stats = array.groupby('Telescope')[column].agg(['mean', 'std'])
        cut_stats['cut'] = cut_stats['mean'] + threshold*cut_stats['std']
        array = array.merge(cut_stats[['cut']], on='Telescope')
        array = array[array[column] < array.cut]
        return array.drop(columns=['cut'])
    elif mode == "value":
        return array[array[column] < threshold]
    else:
        raise ValueError(f"Unknown cut mode for {column}: {mode!r}")


def apply_cuts(
        df,
        min_tel: int = 2,
        min_npix: int = 3,
        telescopes: list = None,
        cuts: dict = None,
):
    """
    Analysis cuts, applied in the order listed below.

    Parameters:
        df: events from build_events()/build_array_events()
        min_tel: minimum number of telescopes in an event (applied before the N_pix cut)
        min_npix: minimum number of pixels in an image
        telescopes: restrict to these telescopes (default = all)
        cuts: optional {column: (mode, threshold)} image cuts, applied in order, see apply_cut(),
            e.g. {"width": ("nsigma", 0.5), "length": ("value", 1.2)}

    Returns:
        the images passing the cuts
    """
    # remove nans and duplicates
    array = df.dropna(subset=["length", "width", "miss", "distance", "alpha"])
    array = array.drop_duplicates(subset=array.drop(["Telescope", "Event", "ImageIndex"], axis=1))
    array = array[array["width"] > 0]

    # cut by number of telescopes
    array = array.groupby('Event', group_keys=False).filter(lambda x: len(x) > (min_tel-1))

    # cut by minimum number of pixels in image
    array = array[array.N_pix >= min_npix]

    # cut by telescope
    if telescopes is not None:
        array = array[array.Telescope.isin(telescopes)]

    for column, (mode, threshold) in (cuts or {}).items():
        array = apply_cut(array, column, mode, threshold)

    return array


def plot_hillas(images, title, colors=None):
    """
    heap.parameterize.plot_hillas_histograms() of images (load_camera_frame() columns), one line per
    telescope plus all telescopes pooled.

    Returns:
        the Figure
    """
    tel_dfs = {name: g for name, g in images.groupby("Telescope")}
    return plot_hillas_histograms(tel_dfs, title=title, colors=colors, pooled_df=images)


@cache
def _cleaned_images(npz_path):
    return np.load(npz_path)["cleaned_images"]


def draw_hillas(ax, tel, color, fill_alpha=0.0):
    """Draw one image's Hillas ellipse (1 sigma length/width) and image axis, in camera degrees."""
    half = CAMERA_HALF_WIDTH
    ax.add_patch(Ellipse((tel.x_c, tel.y_c), 2*tel.length, 2*tel.width, angle=tel.phi,
                         facecolor=to_rgba(color, fill_alpha), edgecolor=color, lw=1.5, label=tel.Telescope))
    t = np.array([-2*half, 2*half])
    phi_rad = np.deg2rad(tel.phi)
    ax.plot(tel.x_c + t*np.cos(phi_rad), tel.y_c + t*np.sin(phi_rad), color=color, lw=0.8, ls="--")


def draw_positions(ax, direction, color, source_xy=None, source=None):
    """Mark the reconstructed direction, and source_xy (camera degrees) if given."""
    ax.scatter(direction.Xoffset, direction.Yoffset, marker="x", color=color, s=80, zorder=3, label="Reconstructed")
    if source_xy is not None:
        ax.scatter(*source_xy, marker="*", facecolor="none", edgecolor=color, lw=1.2, s=200, zorder=3, label=source)


def plot_event(
        event,
        direction,
        output_dir,
        source: str,
        source_position=None,
        pointings: dict = None,
        colors: dict = None,
        rotate_postflip: bool = True,
        pointing_corrections: dict = None,
):
    """
    Each cleaned image of one event in the camera frame, then every image's Hillas ellipse on one
    camera plane, all with the reconstructed direction and the source position.

    Parameters:
        event: one event's images, see apply_cuts()
        direction: that event's row from heap.reconstruction.reconstruct_directions()
        output_dir: processed data dir (holding <date>/), see process_night()
        source: source name, for params_path() and the legend
        source_position: optional SkyCoord of the source; needs pointings to be drawn
        pointings: optional {Run: SkyCoord} each run's pointing, see run_pointings()
        colors: optional {telescope: color} for the combined Hillas panel
        rotate_postflip, pointing_corrections: as passed to build_array_events()

    Returns:
        the Figure
    """
    half = CAMERA_HALF_WIDTH
    colors = colors or {}
    pointing_corrections = pointing_corrections or {}

    # source_position in this run's camera frame (deg from its pointing)
    source_xy = None
    if source_position is not None and pointings is not None and pointings.get(direction.Run) is not None:
        source_xy = make_wcs(pointings[direction.Run]).wcs_world2pix(source_position.ra.deg, source_position.dec.deg, 1)

    fig, axs = plt.subplots(1, len(event) + 1, figsize=(4*(len(event) + 1), 4), squeeze=False)
    for ax, (_, tel) in zip(axs[0], event.iterrows()):
        img = _cleaned_images(params_path(output_dir, tel.Date, tel.Telescope, source))[tel.ImageIndex]
        if rotate_postflip and tel.FlipSide == "postflip":
            img = np.rot90(img, k=2)
        dx, dy = pointing_corrections.get((tel.Date, tel.Telescope, tel.FlipSide), (0, 0))
        im = ax.imshow(img, origin="lower", extent=[-half-dx, half-dx, -half-dy, half-dy], cmap=CAMERA_CMAP, vmin=0)
        fig.colorbar(im, ax=ax, label="ADU", fraction=0.046, pad=0.04)
        draw_hillas(ax, tel, "white")
        draw_positions(ax, direction, "white", source_xy, source)
        ax.set(xlim=[-half, half], ylim=[-half, half], xlabel="X (degrees)", ylabel="Y (degrees)")
        ax.set_title(f"{tel.Telescope} ({tel.FlipSide})", fontsize=9)

    ax = axs[0, -1]
    for i, (_, tel) in enumerate(event.iterrows()):
        draw_hillas(ax, tel, colors.get(tel.Telescope, TAB10_COLORS[i % len(TAB10_COLORS)]), fill_alpha=0.3)
    draw_positions(ax, direction, "k", source_xy, source)
    ax.set(xlim=[-half, half], ylim=[-half, half], aspect="equal", xlabel="X (degrees)", ylabel="Y (degrees)")
    ax.set_title("Hillas parameterization", fontsize=9)
    ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.02, 1), borderaxespad=0)
    fig.suptitle(f"Event {direction.Event}, {pd.to_datetime(event.Timestamp.iloc[0], unit='s', utc=True):%Y-%m-%d %H:%M:%S.%f} UTC", fontsize=10)
    fig.tight_layout()
    return fig
