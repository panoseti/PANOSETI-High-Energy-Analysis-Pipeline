"""events

Builds array events (images from 2+ telescopes of the same shower) out of per-telescope Hillas
parameters (<output_dir>/<date>/<telescope>/<source_slug>/<source_slug>.npz, see process_night()
and heap.process_dataset.process_dataset()), applies the analysis cuts, and plots single events.

Hillas parameters are in degrees on the camera (32 pixels over +-4.95 deg, origin at camera
center, see heap.parameterize.calc_params()). x_c runs along the columns and y_c
along the rows of the (32, 32) image, so x_c/y_c here are the transpose of the ROOT CSVs'
MeanX/MeanY (ROOT fills bin (i+1, j+1) from pixel [i][j]).

Every telescope in a run points at the same target (Ekos aligns each one to it before the run, then
guides), so each run has one pointing, shared by every telescope: the mount positions at the run's
start (run_pointings()). Later mount positions drift away from it as they follow the guide
corrections. pointing_corrections correct each camera's offset from that pointing
(load_camera_frame()). Only timing uses a reference telescope (build_events()).
"""
import json
from functools import cache
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from matplotlib.colors import to_rgba
from matplotlib.patches import Ellipse

from heap import coincidences as coinc
from heap.parameterize import CAMERA_CMAP, CAMERA_HALF_WIDTH, TAB10_COLORS, plot_hillas_histograms
from heap.significance import make_wcs
from heap.process_dataset import _mount_pointing, _run_start_epoch, discover_runs, group_runs_by_source, load_fallback_map, process_dataset, slugify



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
    under the same output_dir lets process_dataset() fall back to another night's gain map.

    Each telescope's settings (data_product, rate_cut, and the cleaning settings) are recorded in
    <output_dir>/<date>/<telescope>/processing.json; existing output made with other (or
    unrecorded) settings raises instead of being reused, unless reprocess.

    Parameters:
        raw_dir: this night's raw data dir, containing one subfolder per run
        output_dir: processed data dir (holding <date>/)
        date: night, YYYYMMDD
        telescopes: {module: {"name": telescope_name, "rate_cut": multiple of median rate}}, e.g. {"module_253": {"name": "Winter", "rate_cut": 3}}
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
            "data_product": data_product, "rate_cut": info.get("rate_cut", 3),
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
                rate_cut=info.get("rate_cut", 3),
                image_threshold=image_threshold, border_threshold=border_threshold,
                keep_brightest_island=keep_brightest_island,
            )
            settings_path.write_text(json.dumps(settings, indent=2))
        except Exception as e:
            print(f"{date}: processing {name} failed: {e}")


def source_runs(raw_dir, telescope, source):
    """
    One night's runs of source, by flip side, from the mount logs (see
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


def run_pointings(raw_dir, source, telescopes: list, pointing_overrides: dict = None):
    """
    Each run of source on one night's pointing, shared by every telescope: the mount positions at the
    run's start (heap.process_dataset._mount_pointing()), averaged over the telescopes. Each run
    starts right after Ekos aligns every telescope to the target, and guiding then holds it there.

    pointing_overrides replace this for runs that didn't start with an alignment (e.g. a DAQ restart
    mid-observation) or have no mount data; runs with neither are left out, with a printout.

    Parameters:
        raw_dir: this night's raw data dir, see source_runs()
        source: source name
        telescopes: telescope names
        pointing_overrides: optional {run folder name: SkyCoord}

    Returns:
        {run folder name: SkyCoord}
    """
    starts = {} # {run: each telescope's mount position at its start}
    for telescope in telescopes:
        for run_dirs in source_runs(raw_dir, telescope, source).values():
            for run_dir in run_dirs:
                positions = starts.setdefault(Path(run_dir).name, [])
                position = _mount_pointing(run_dir, telescope)
                if position is not None:
                    positions.append(position)

    pointings = {}
    for run, positions in sorted(starts.items()):
        if run in (pointing_overrides or {}):
            pointings[run] = pointing_overrides[run]
        elif positions:
            mean = SkyCoord(SkyCoord(positions).cartesian.mean(), frame="icrs")
            pointings[run] = SkyCoord(mean.ra, mean.dec)
        else:
            print(f"{run}: no mount position at its start and no pointing_overrides, leaving it out")
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
    Load one telescope's Hillas parameters for one night, in camera coordinates (degrees), sorted by
    time.

    Parameters:
        npz_path: <source_slug>.npz written by heap.process_dataset.process_dataset()
        telescope: telescope name, stored in the Telescope column
        runs_by_flip: this night's runs of the source, see source_runs()
        rotate_postflip: rotate postflip images by 180 deg to the preflip orientation
        rel_efficiency: relative telescope efficiency; size is divided by it
        pointing_corrections: optional {flip_side: (dx, dy)} in deg, corrects this camera's
            misalignment, subtracted from x_c/y_c: where a sky position appears in this camera minus
            where it should appear given the run's pointing (run_pointings(); +x east, +y south).
            E.g. pointed at the Crab but it appears at (1, 1): (dx, dy) = (1, 1). Postflip offsets
            are subtracted after the 180 deg rotation and are NOT rotated themselves, so they must
            already be in rotated coordinates: measured on rotated images, or measured on raw
            postflip images and negated by the caller. Not checked.

    Returns:
        DataFrame with ImageIndex (row in the npz's cleaned_images, which holds every run of the source that night), Telescope, Timestamp,
        FlipSide, Run, the npz's Hillas parameters (x_c, y_c, phi, size, N_pix, length, width, miss, distance, alpha),
        and x_shift/y_shift, how far pointing_corrections moved the centroid from its position in the (rotated) camera image
    """
    npz = np.load(npz_path)
    p = pd.DataFrame({col: npz[col] for col in npz.files if col != "cleaned_images"})
    p = p.sort_values("Timestamp", ignore_index=True)

    df = pd.DataFrame({
        "ImageIndex": p.Event,
        "Telescope": telescope,
        "Timestamp": p.Timestamp,
        "FlipSide": flip_side_of(p.Timestamp.to_numpy(), runs_by_flip),
        "Run": run_of(p.Timestamp.to_numpy(), runs_by_flip),
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

    x_camera, y_camera = df.x_c.to_numpy().copy(), df.y_c.to_numpy().copy()
    # correct this camera's misalignment (camera vs the run's pointing); postflip offsets are
    # assumed to already be in rotated coordinates, see pointing_corrections above
    for side, (dx, dy) in (pointing_corrections or {}).items():
        df.loc[df.FlipSide == side, "x_c"] -= dx
        df.loc[df.FlipSide == side, "y_c"] -= dy

    df["x_shift"] = df.x_c - x_camera
    df["y_shift"] = df.y_c - y_camera
    return df


TIMING_BIN_WIDTH = 120 # s, coincidences.correct_time()'s default bin_width


def _timing_segments(images: dict, order: list, bin_width: float = TIMING_BIN_WIDTH):
    """
    Split one run's images into stretches of time with the same timing reference: in each
    bin_width (s) bin, the first telescope in order with images in it. A bin with no images joins
    the stretch before it.

    Returns:
        [(reference, start, end, {telescope: images}), ...], start/end Unix times (s)
    """
    start = min(df.Timestamp.min() for df in images.values())
    bins = {name: ((df.Timestamp.to_numpy() - start) // bin_width).astype(int) for name, df in images.items()}
    n_bins = max(b.max() for b in bins.values()) + 1
    has_images = {name: np.bincount(b, minlength=n_bins) > 0 for name, b in bins.items()}
    references = []
    for i in range(n_bins):
        reference = next((name for name in order if name in has_images and has_images[name][i]), None)
        references.append(reference or references[-1]) # bin 0 always has images

    segments = []
    first = 0
    for i in range(1, n_bins + 1):
        if i == n_bins or references[i] != references[first]:
            segment = {name: df[(bins[name] >= first) & (bins[name] < i)].reset_index(drop=True) for name, df in images.items()}
            segments.append((references[first], start + first*bin_width, start + i*bin_width, {name: df for name, df in segment.items() if len(df)}))
            first = i
    return segments


def _match_segment(telescopes: dict, reference: str, window: float, plotting: bool, label: str = ""):
    """build_events() for one stretch of time, with reference as its timing reference; label is
    appended to correct_time's plot_name, e.g. the run and time range."""
    ref_timestamps = telescopes[reference].Timestamp.to_numpy()
    for name, df in list(telescopes.items()):
        if name == reference:
            continue
        df = df.copy()
        try:
            df["Timestamp"] = coinc.correct_time(df.Timestamp.to_numpy(), ref_timestamps, plot_name=f"{reference}-{name}{label}", base_dir=None, plotting=plotting, coinc_window=window)
        except ValueError:
            print(f"{name}: no coincidences with {reference} in this stretch, not matched into events here")
            del telescopes[name]
            continue
        # timestamps with no timing correction (no coincidences in their bin) come back nan
        telescopes[name] = df.dropna(subset=["Timestamp"]).sort_values("Timestamp", ignore_index=True)
    if len(telescopes) < 2:
        return None

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


def build_events(telescopes: dict, reference: list, window: float = 0.001, plotting: bool = False):
    """
    Match one night's images across telescopes into events, run by run.

    Each run is split into TIMING_BIN_WIDTH bins, and each bin's timing reference is the first
    telescope in reference with images in that bin, else the first other telescope that has some
    (in the order of telescopes). So a run keeps its events when the preferred timing reference
    has data for only part of it. In each stretch of bins with the same timing reference, the other
    telescopes' timestamps are corrected against it (see coincidences.correct_time), leaving out a
    telescope with no coincidences with it, then matched pairwise and merged (see
    coincidences.find_coincidences). Events where a telescope matched more than one image are
    dropped as ambiguous.

    Parameters:
        telescopes: {telescope_name: DataFrame from load_camera_frame()}
        reference: list of timing reference telescopes, in order of preference
        window: coincidence window (s)
        plotting: show correct_time's before/after plots

    Returns:
        DataFrame of every image in an event (telescopes' columns plus Event, numbered from 0),
        sorted by Event then Telescope, or None if no run has images from 2+ telescopes
    """
    telescopes = {name: df for name, df in telescopes.items() if df is not None and len(df)}
    order = list(reference)
    order += [name for name in telescopes if name not in order]

    parts = []
    for run in sorted(set().union(*(df.Run for df in telescopes.values()))):
        run_images = {name: df[df.Run == run].reset_index(drop=True) for name, df in telescopes.items() if (df.Run == run).any()}
        if len(run_images) < 2:
            continue
        for segment_reference, start, end, images in _timing_segments(run_images, order):
            if len(images) < 2:
                continue
            span = f"{pd.to_datetime(start, unit='s'):%H:%M:%S}-{pd.to_datetime(end, unit='s'):%H:%M:%S} UTC"
            skipped = order[:order.index(segment_reference)] # no images in this stretch
            why = f" (no {', '.join(skipped)} images)" if skipped else ""
            counts = ", ".join(f"{name} {len(df)}" for name, df in images.items())
            print(f"{run}, {span}: timing reference {segment_reference}{why}; {counts} images")
            label = f", run {pd.to_datetime(_run_start_epoch(Path(run)), unit='s'):%H:%M:%S}, {span}"
            events = _match_segment(images, segment_reference, window, plotting, label)
            print() # blank line between stretches
            if events is None or not len(events):
                continue
            events["Event"] += parts[-1].Event.max() + 1 if parts else 0
            parts.append(events)

    return pd.concat(parts, ignore_index=True) if parts else None


def load_camera_frames(
        output_dir,
        raw_dir,
        date: str,
        source: str,
        telescopes: list,
        rotate_postflip: bool = True,
        rel_tel_efficiency: dict = None,
        pointing_corrections: dict = None,
        pointing_overrides: dict = None,
):
    """
    load_camera_frame() for every telescope with processed data for source on date, keeping only
    runs with a pointing (run_pointings()).

    Parameters:
        output_dir: processed data dir (holding <date>/), see process_night()
        raw_dir: this night's raw data dir, for flip sides and mount positions (see source_runs())
        pointing_overrides: see run_pointings()
        rel_tel_efficiency: optional {telescope: efficiency}, see load_camera_frame()
        pointing_corrections: optional {(date, telescope, flip_side): (dx, dy)} in deg, each camera's
            offset from the run's pointing, see load_camera_frame()
        rotate_postflip: see load_camera_frame()

    Returns:
        {telescope: DataFrame from load_camera_frame()}, {run folder name: SkyCoord} each run's pointing
    """
    pointings = run_pointings(raw_dir, source, telescopes, pointing_overrides)
    frames = {}
    for name in telescopes:
        path = params_path(output_dir, date, name, source)
        if not path.exists():
            continue
        corrections = {
            side: offset for (d, tel, side), offset in (pointing_corrections or {}).items()
            if d == date and tel == name
        }
        frame = load_camera_frame(
            path, name, source_runs(raw_dir, name, source),
            rotate_postflip=rotate_postflip,
            rel_efficiency=(rel_tel_efficiency or {}).get(name, 1.0),
            pointing_corrections=corrections,
        )
        frames[name] = frame[frame.Run.isin(pointings)].reset_index(drop=True) # runs left out by run_pointings() say so there
    return frames, pointings


def build_array_events(
        output_dir,
        raw_dir,
        date: str,
        source: str,
        telescopes: list,
        reference: list,
        window: float = 0.001,
        rotate_postflip: bool = True,
        rel_tel_efficiency: dict = None,
        pointing_corrections: dict = None,
        pointing_overrides: dict = None,
        plotting: bool = False,
):
    """
    load_camera_frames() then build_events(). Adds a Date column holding date; x_c/y_c are relative
    to each event's run's pointing (see run_pointings()).

    Parameters:
        output_dir, raw_dir, rotate_postflip, rel_tel_efficiency, pointing_corrections, pointing_overrides: see load_camera_frames()
        window, reference, plotting: see build_events()

    Returns:
        events (see build_events()), {run folder name: SkyCoord} each run's pointing (see run_pointings())
    """
    frames, pointings = load_camera_frames(
        output_dir, raw_dir, date, source, telescopes, rotate_postflip=rotate_postflip,
        rel_tel_efficiency=rel_tel_efficiency, pointing_corrections=pointing_corrections,
        pointing_overrides=pointing_overrides,
    )
    events = build_events(frames, reference, window=window, plotting=plotting)
    if events is not None:
        events.insert(1, "Date", date)
        events.insert(2, "Run", events.pop("Run")) # every image of an event is in the same run, see build_events()
    return events, pointings


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
):
    """
    Each cleaned image of one event, shifted with its centroid (x_shift/y_shift, see
    load_camera_frame()), then every image's Hillas ellipse on one camera plane, all with the
    reconstructed direction and the source position.

    Parameters:
        event: one event's images, see apply_cuts()
        direction: that event's row from heap.reconstruction.reconstruct_directions()
        output_dir: processed data dir (holding <date>/), see process_night()
        source: source name, for params_path() and the legend
        source_position: optional SkyCoord of the source; needs pointings to be drawn
        pointings: optional {Run: SkyCoord} each run's pointing, see run_pointings()
        colors: optional {telescope: color} for the combined Hillas panel
        rotate_postflip: as passed to build_array_events()

    Returns:
        the Figure
    """
    half = CAMERA_HALF_WIDTH
    colors = colors or {}

    # source_position in camera coordinates relative to this run's pointing
    source_xy = None
    if source_position is not None and pointings is not None and pointings.get(direction.Run) is not None:
        source_xy = make_wcs(pointings[direction.Run]).wcs_world2pix(source_position.ra.deg, source_position.dec.deg, 1)

    fig, axs = plt.subplots(1, len(event) + 1, figsize=(4*(len(event) + 1), 4), squeeze=False)
    for ax, (_, tel) in zip(axs[0], event.iterrows()):
        img = _cleaned_images(params_path(output_dir, tel.Date, tel.Telescope, source))[tel.ImageIndex]
        if rotate_postflip and tel.FlipSide == "postflip":
            img = np.rot90(img, k=2)
        # shifted with its centroid by the pointing correction
        im = ax.imshow(img, origin="lower", extent=[-half+tel.x_shift, half+tel.x_shift, -half+tel.y_shift, half+tel.y_shift], cmap=CAMERA_CMAP, vmin=0)
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
