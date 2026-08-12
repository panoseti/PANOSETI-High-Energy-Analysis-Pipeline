"""process_dataset

Runs process_image() over every event in a dataset, keeping every cleaned image
and collecting Hillas parameters into a table.
"""
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import pypff

from heap.coincidences import load_telescope_tv
from heap.image_cleaning import threshold_clean
from heap.make_gain_map import gain_from_pedvars
from heap.make_pedestals import calculate_pedestal_and_pedvar
from heap.parameterize import calc_params


def slugify(text, sep="-"):
    """Replaces every run of non-alphanumeric characters in text with sep, and strips
    sep from both ends (e.g. slugify("MGRO J2019+37", sep="_") -> "MGRO_J2019_37")."""
    return re.sub(r"[^A-Za-z0-9]+", sep, text).strip(sep)


def process_image(
        pixdata: np.ndarray,
        pedestal: np.ndarray,
        pedvar: np.ndarray,
        gain: np.ndarray,
        x: float = None,
        y: float = None,
        image_threshold: float = 4.0,
        border_threshold: float = 2.0,
):
    """
    Pedestal-subtract, threshold-clean, gain-correct, and parameterize a single camera image.

    signal = (pixdata - pedestal) / gain
    nsigma = (pixdata - pedestal) / pedvar

    gain is excluded from the significance test (nsigma); threshold_clean() decides which
    pixels survive using pedvar alone, and gain is applied here afterward, only to those
    surviving pixels' amplitude.

    pedestal, pedvar, and gain are all calibration products computed once for a whole
    dataset (see make_pedestals.calculate_pedestal_and_pedvar and
    make_gain_map.gain_from_pedvars) and then reused across every image in that dataset -
    they are not recomputed per image here.

    Parameters:
        pixdata: raw ADC values, shape (32, 32) or (1024,)
        pedestal: per-pixel pedestal, shape (32, 32) or (1024,)
        pedvar: per-pixel pedestal variance used for the significance test, shape (32, 32) or (1024,)
        gain: per-pixel relative gain map, shape (32, 32) or (1024,)
        x, y: test position (pixels) passed through to calc_params(), e.g. for distance/alpha.
            Default value is camera center.
        image_threshold: image pixel threshold, in units of pedvar
        border_threshold: border pixel threshold, in units of pedvar

    Returns:
        cleaned: the cleaned, gain-corrected image, shape (32, 32)
        params: dict of Hillas parameters from calc_params()
    """

    pixdata = np.asarray(pixdata, dtype=float).reshape(32, 32)
    pedestal = np.asarray(pedestal, dtype=float).reshape(32, 32)
    pedvar = np.asarray(pedvar, dtype=float).reshape(32, 32)
    gain = np.asarray(gain, dtype=float).reshape(32, 32)

    ped_subtracted = pixdata - pedestal

    cleaned, _ = threshold_clean(
        ped_subtracted,
        pedvar,
        image_threshold=image_threshold,
        border_threshold=border_threshold,
    )

    cleaned = cleaned / gain

    params = calc_params(cleaned, x, y)

    return cleaned, params


def _run_start_epoch(run_dir):
    """Parses a run folder's start time (e.g. ...start_2026-07-09T07:22:57Z...) into a Unix
    epoch (seconds), comparable against frame timestamps. Returns None if the name doesn't
    match the expected pattern (see discover_runs())."""
    match = re.search(r"start_(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2})Z", run_dir.name)
    if not match:
        return None
    return datetime.strptime(match.group(1), "%Y-%m-%dT%H:%M:%S").replace(tzinfo=timezone.utc).timestamp()


def _closest_alt_run(source, reference_epoch, groups):
    """Finds the single run, out of every other source (besides `source`) in `groups`, whose
    start time is closest to reference_epoch, across either flip side.

    Returns (run_dir, other_source, flip_side), or None if `groups` has no other source with a
    parseable run start time."""
    best = None
    best_dt = None
    for other_source, runs_by_flip in groups.items():
        if other_source == source:
            continue
        for flip_side, run_dirs in runs_by_flip.items():
            for run_dir in run_dirs:
                epoch = _run_start_epoch(run_dir)
                if epoch is None:
                    continue
                dt = abs(epoch - reference_epoch)
                if best_dt is None or dt < best_dt:
                    best_dt, best = dt, (run_dir, other_source, flip_side)
    return best


def _find_prior_night_gain(telescope_out_dir):
    """Looks for the most recent earlier night's already-computed gain map for this telescope
    (any source), to reuse when this night has no usable pair of star fields at all.

    telescope_out_dir is <output_root>/<date>/<telescope>/ (what process_dataset() calls
    out_dir); sibling nights are <output_root>/<other_date>/<telescope>/*/calibrations.npz.
    A candidate night is skipped if its own gain map was itself a "flat" fallback, so flatness
    doesn't propagate forward night after night.

    Returns (gain_map (32, 32) float32, date string), or (None, None) if nothing usable is found.
    """
    telescope_name = telescope_out_dir.name
    current_date = telescope_out_dir.parent.name
    output_root = telescope_out_dir.parent.parent

    if not output_root.is_dir():
        return None, None

    candidate_dates = sorted(
        p.name for p in output_root.iterdir()
        if p.is_dir() and p.name < current_date and (p / telescope_name).is_dir()
    )

    for date in reversed(candidate_dates):
        for calib_path in sorted((output_root / date / telescope_name).glob("*/calibrations.npz")):
            calib = np.load(calib_path)
            gain_source = calib["gain_source"].item() if "gain_source" in calib else "own"
            if gain_source == "flat":
                continue
            return calib["gain"].astype(np.float32), date

    return None, None


def resolve_gain_map(
        source: str,
        data_preflip: np.ndarray,
        timestamps_preflip: np.ndarray,
        data_postflip: np.ndarray,
        timestamps_postflip: np.ndarray,
        groups: dict,
        module_pattern: str,
        telescope_out_dir,
        rate_cut: float,
        nsig: float = 5.0,
        fit_gaussian: bool = True,
):
    """
    Determines a source's per-pixel relative gain map, trying progressively weaker fallbacks:

      1. own night, own source: compare this source's own preflip vs postflip pooled pedvar
         (two star fields taken through the same optics minimize the impact of any single
         bright star; see make_gain_map.gain_from_pedvars).
      2. own night, other source: if only one flip side is available for this source, pair it
         against the same night's single temporally-closest run belonging to another source
         (either flip side) as a stand-in second star field - the pre/post-flip labeling
         doesn't matter to gain_from_pedvars, it just needs two different pointings to compare.
      3. another night: if this source is the only one in the whole night (no other source to
         pair with), reuse the most recent earlier night's already-computed gain map for this
         telescope (see _find_prior_night_gain()).
      4. flat ones: if none of the above apply (e.g. the first night this telescope is processed).

    Returns:
        gain_map: (32, 32) float32
        gain_source: a short machine-checkable tag for which tier was used - "own", "alt",
            "prior_night:<date>", or "flat" (the exact string "flat" is checked elsewhere to keep
            flatness from propagating forward night after night; see _find_prior_night_gain()).
        gain_caption: a human-readable sentence naming the actual preflip/postflip frames (or
            explaining why there weren't two comparable ones), saved into calibrations.npz
            alongside gain_source so the gallery can show it as a caption.
    """
    if len(data_preflip) > 0 and len(data_postflip) > 0:
        _, pedvar_preflip = calculate_pedestal_and_pedvar(data_preflip, nsig=nsig, fit_gaussian=fit_gaussian)
        _, pedvar_postflip = calculate_pedestal_and_pedvar(data_postflip, nsig=nsig, fit_gaussian=fit_gaussian)
        gain_map = gain_from_pedvars(pedvar_preflip, pedvar_postflip).astype(np.float32)
        return gain_map, "own", f"Preflip: {source} · Postflip: {source}"

    own_data = data_preflip if len(data_preflip) > 0 else data_postflip
    own_timestamps = timestamps_preflip if len(data_preflip) > 0 else timestamps_postflip
    own_flip_side = "preflip" if len(data_preflip) > 0 else "postflip"

    if len(own_data) > 0:
        _, own_pedvar = calculate_pedestal_and_pedvar(own_data, nsig=nsig, fit_gaussian=fit_gaussian)

        alt = _closest_alt_run(source, own_timestamps.mean(), groups)
        if alt is not None:
            alt_run, alt_source, alt_flip_side = alt
            alt_data, _ = _load_runs([alt_run], module_pattern, rate_cut)
            if len(alt_data) > 0:
                _, alt_pedvar = calculate_pedestal_and_pedvar(alt_data, nsig=nsig, fit_gaussian=fit_gaussian)
                gain_map = gain_from_pedvars(own_pedvar, alt_pedvar).astype(np.float32)
                caption = (
                    f"{own_flip_side.capitalize()}: {source} · "
                    f"{alt_flip_side.capitalize()}: {alt_source} (different field)"
                )
                return gain_map, "alt", caption

    prior_gain, prior_date = _find_prior_night_gain(telescope_out_dir)
    if prior_gain is not None:
        return prior_gain, f"prior_night:{prior_date}", f"Reused from prior night {prior_date} (no usable frame pair this night)"

    return np.ones((32, 32), dtype=np.float32), "flat", "Flat gain (1.0) — no usable frame pair this night or a prior night"


def build_calibrations(
        data_preflip: np.ndarray,
        timestamps_preflip: np.ndarray,
        data_postflip: np.ndarray,
        timestamps_postflip: np.ndarray,
        out_dir,
        gain_map: np.ndarray,
        gain_source: str,
        gain_caption: str,
        time_window: float = 600,
        max_gap: float = 300,
        nsig: float = 5.0,
        fit_gaussian: bool = True,
):
    """
    Build pedestal and pedvar calibrations for a source's night, alongside an already-resolved
    gain map (see resolve_gain_map()).

    Pedestal and pedvar are computed per time_window-second chunk, respecting continuous
    observing segments (a gap larger than max_gap starts a new segment; see
    make_pedestals.calculate_pedestal_and_pedvar for the per-chunk calculation). A chunk with
    too few frames (< 100) reuses the previous chunk's pedestal/pedvar within its segment; if
    there is no previous chunk yet, those frames' pedestals are left at zero.

    Writes <out_dir>/calibrations.npz holding pedestals, pedestal_variances (both
    shape (n_pre + n_post, 1024)), gain (shape (1024,), since it is constant across every frame
    in a source's night - it is not repeated per frame on disk), gain_source, and gain_caption.

    Parameters:
        data_preflip, timestamps_preflip: pre-flip run data frames (n_pre, 1024) and their Unix
            timestamps in seconds (n_pre,). Pass empty arrays if no pre-flip run is available.
        data_postflip, timestamps_postflip: post-flip run data frames (n_post, 1024) and their
            Unix timestamps in seconds (n_post,). Pass empty arrays if no post-flip run is available.
        out_dir: directory to write calibrations.npz to (created if missing)
        gain_map: (32, 32) relative gain map, as returned by resolve_gain_map()
        gain_source: label for how gain_map was resolved, as returned by resolve_gain_map()
        gain_caption: human-readable caption for how gain_map was resolved, as returned by
            resolve_gain_map()
        time_window: width in seconds of each pedestal chunk (default = 600)
        max_gap: gap in seconds between frames that starts a new observing segment (default = 300)
        nsig: number of sigmas above mean to define the pedvar outlier mask (default = 5)
        fit_gaussian: whether to fit a Gaussian for the pedvar (default = True), see
            make_pedestals.calculate_pedestal_and_pedvar

    Returns:
        data: pre-flip and post-flip data frames concatenated, shape (n_pre + n_post, 1024)
        timestamps: pre-flip and post-flip timestamps concatenated, shape (n_pre + n_post,)
        pedestals: shape (n_pre + n_post, 1024)
        pedestal_variances: shape (n_pre + n_post, 1024)
        gain: shape (n_pre + n_post, 1024)
    """

    gain_map = np.asarray(gain_map, dtype=np.float32).reshape(32, 32)

    data = np.concatenate([data_preflip, data_postflip], axis=0)
    timestamps = np.concatenate([timestamps_preflip, timestamps_postflip], axis=0)

    gain = np.broadcast_to(gain_map.reshape(1, 1024), (len(data), 1024))

    sort_idx = np.argsort(timestamps)
    sorted_times = timestamps[sort_idx]
    sorted_data = data[sort_idx]

    time_diffs = np.diff(sorted_times)
    gap_indices = np.where(time_diffs > max_gap)[0]

    segment_starts = [0] + (gap_indices + 1).tolist()
    segment_ends = gap_indices.tolist() + [len(sorted_times)]

    print(f"  Found {len(segment_starts)} continuous observing segments")

    pedestals = np.zeros((len(data), 1024), dtype=np.float32)
    pedestal_variances = np.zeros((len(data), 1024), dtype=np.float32)

    for seg_start, seg_end in zip(segment_starts, segment_ends):
        seg_times = sorted_times[seg_start:seg_end]
        seg_data = sorted_data[seg_start:seg_end]

        if len(seg_data) < 100:
            print(f"    Skipping short segment with only {len(seg_data)} frames")
            continue

        seg_date = datetime.fromtimestamp(seg_times[0]).strftime('%Y%m%d')
        print(f"    Processing segment from {seg_date}, {len(seg_data)} frames")

        start_time = seg_times[0]
        end_time = seg_times[-1]

        last_ped, last_pedvar = None, None
        current_time = start_time
        while current_time < end_time:
            mask = (seg_times >= current_time) & (seg_times < current_time + time_window)
            chunk_data = seg_data[mask]
            chunk_indices = np.where(mask)[0] + seg_start  # Global indices
            window_str = datetime.fromtimestamp(current_time).strftime('%H:%M:%S')

            if len(chunk_data) > 100:  # Need enough frames to make a good pedestal
                ped, ped_var = calculate_pedestal_and_pedvar(chunk_data, nsig=nsig, fit_gaussian=fit_gaussian)
                last_ped, last_pedvar = ped, ped_var
                print(f"      Created pedestal at {window_str} using {len(chunk_data)} frames")
            elif last_ped is not None:
                ped, ped_var = last_ped, last_pedvar
                print(f"      Carrying forward previous pedestal for window at {window_str}, only {len(chunk_data)} frames")
            else:
                print(f"      Skipping window at {window_str}, only {len(chunk_data)} frames, no previous pedestal to carry forward")
                current_time += time_window
                continue

            for idx in chunk_indices:
                pedestals[sort_idx[idx]] = ped
                pedestal_variances[sort_idx[idx]] = ped_var

            current_time += time_window

    print(f"  Assigned pedestals to {len(data)} frames")

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    np.savez(
        out_dir / "calibrations.npz",
        pedestals=pedestals,
        pedestal_variances=pedestal_variances,
        gain=gain_map,
        gain_source=gain_source,
        gain_caption=gain_caption,
    )

    return data, timestamps, pedestals, pedestal_variances, gain


def discover_runs(raw_dir):
    """
    List a night's run subfolders, sorted chronologically by start time.

    Run folder names encode their start time (e.g.
    obs_Palomar.start_2026-07-09T07:22:57Z.runtype_obs-test.pffd), so a lexicographic sort
    is already a chronological sort.
    """
    return sorted(p for p in Path(raw_dir).iterdir() if p.is_dir())


def _get_mount_hk(run_dir, telescope):
    """Returns hk.pff's MOUNT_<TELESCOPE> table for this run, or None if it's not present
    (e.g. telescope has no tracked mount or it was never written)."""
    return pypff.PanosetiRun(run_dir).get_hk().get(f"MOUNT_{telescope.upper()}")


def identify_source(run_dir, telescope, fallback_map=None):
    """
    Identify which astronomical source a run was tracking.

    Tries target_name from hk.pff's MOUNT_<TELESCOPE> table first. 
    Falls back to fallback_map when it's missing or blank.

    Parameters:
        run_dir: path to one run folder
        telescope: telescope name, used to look up hk.pff's MOUNT_<TELESCOPE> table
        fallback_map: optional {run_dir_name: {"source": ..., "flip_side": ...}} dict (see
            load_fallback_map()), used when hk.pff's target_name isn't usable for this run

    Returns:
        source_name
    """
    mount = _get_mount_hk(run_dir, telescope)
    if mount is not None:
        names = {n for n in mount["target_name"] if n}
        if len(names) == 1:
            return names.pop()

    if fallback_map is not None and run_dir.name in fallback_map:
        return fallback_map[run_dir.name]["source"]

    raise ValueError(f"Could not identify source for {run_dir}: no usable hk target_name and no fallback_map entry")


def identify_flip_side(run_dir, telescope, fallback_map=None):
    """
    Identify whether a run is pre-flip or post-flip.

    Tries side from hk.pff's MOUNT_<TELESCOPE> table first. 
    Falls back to fallback_map when it's missing or ambiguous.

    Parameters:
        run_dir: path to one run folder
        telescope: telescope name, used to look up hk.pff's MOUNT_<TELESCOPE> table
        fallback_map: optional {run_dir_name: {"source": ..., "flip_side": ...}} dict (see
            load_fallback_map()), used when hk.pff's side isn't usable for this run

    Returns:
        "preflip" or "postflip"
    """
    mount = _get_mount_hk(run_dir, telescope)
    if mount is not None:
        sides = {s for s in mount["side"] if s}
        if len(sides) == 1:
            side = sides.pop().upper()
            if side == "WEST":
                return "preflip"
            if side == "EAST":
                return "postflip"

    if fallback_map is not None and run_dir.name in fallback_map:
        return fallback_map[run_dir.name]["flip_side"]

    raise ValueError(f"Could not identify flip side for {run_dir}: no usable hk side and no fallback_map entry")


def load_fallback_map(path):
    """
    Load a user-maintained fallback map for identify_source()/identify_flip_side(), for runs
    where hk.pff's mount data is missing.

    Expects a JSON file: {run_dir_name: {"source": source_name, "flip_side": "preflip" | "postflip"}}
    """
    with open(path) as f:
        return json.load(f)


def group_runs_by_source(raw_dir, telescope, fallback_map=None):
    """
    Group a night's runs by source and flip side.

    Parameters:
        raw_dir: path to a night's raw data folder, containing one subfolder per run
        telescope: telescope name, passed through to identify_source()/identify_flip_side()
        fallback_map: optional {run_dir_name: {"source": ..., "flip_side": ...}} dict (see
            load_fallback_map()), used when hk.pff's mount data is missing

    Returns:
        {source_name: {"preflip": [run_dir, ...], "postflip": [run_dir, ...]}}
    """
    groups = {}
    for run_dir in discover_runs(raw_dir):
        source = identify_source(run_dir, telescope, fallback_map=fallback_map)
        flip_side = identify_flip_side(run_dir, telescope, fallback_map=fallback_map)
        groups.setdefault(source, {"preflip": [], "postflip": []})[flip_side].append(run_dir)
    return groups


def _load_runs(run_dirs, module_pattern, rate_cut):
    """Load + clean one telescope's data across a list of run folders (see
    coincidences.load_telescope_tv for the per-run read_pff -> cut_pkt_loss_old ->
    spike_cut pipeline)."""

    if not run_dirs:
        return np.empty((0, 1024)), np.empty((0,))

    all_data, all_timestamps = [], []
    for run_dir in run_dirs:
        data, timestamps = load_telescope_tv(module_pattern, run_dir, rate_cut)
        all_data.append(data)
        all_timestamps.append(timestamps)

    return np.concatenate(all_data, axis=0), np.concatenate(all_timestamps, axis=0)


def process_dataset(
        raw_dir,
        out_dir,
        module_pattern: str,
        telescope: str,
        fallback_map_path=None,
        rate_cut: float = 20,
        time_window: float = 600,
        max_gap: float = 300,
        nsig: float = 5.0,
        fit_gaussian: bool = True,
        x: float = None,
        y: float = None,
        image_threshold: float = 4.0,
        border_threshold: float = 2.0,
):
    """
    Run the full pipeline for one telescope's data in a night's raw_dir.

    Discovers the night's runs and groups them by source and flip side (see
    group_runs_by_source()), resolves each source's gain map (see resolve_gain_map()) and
    builds its pedestal/pedvar calibrations (see build_calibrations()), then runs
    process_image() over every event, keeping every cleaned image and collecting Hillas
    parameters into a table. Writes one <out_dir>/<source_slug>/<source_slug>.npz per source
    (source_slug is source with non-alphanumeric characters, e.g. spaces, replaced by
    underscores; see slugify()), holding cleaned_images plus one array per params_df column,
    alongside that same source's calibrations.npz (see build_calibrations()).

    Parameters:
        raw_dir: path to a night's raw data folder, containing one subfolder per run
        out_dir: directory to write each source's <source>.npz to (created if missing)
        module_pattern: glob pattern identifying this telescope's files within a run folder,
            e.g. "dp_ph1024*module_252"
        telescope: telescope name, used to look up hk.pff's MOUNT_<TELESCOPE> table (see
            identify_source()/identify_flip_side())
        fallback_map_path: optional path to a fallback map file (see load_fallback_map()),
            used when hk.pff's mount data is missing
        rate_cut: spike_cut's trigger-rate threshold (Hz), see coincidences.load_telescope_tv
        time_window, max_gap, nsig, fit_gaussian: passed through to build_calibrations()
        x, y, image_threshold, border_threshold: passed through to process_image()

    Returns:
        {source_name: (cleaned_images, params_df, raw_images, data_preflip, timestamps_preflip,
        data_postflip, timestamps_postflip)}, one entry per source found in raw_dir - see
        build_calibrations() and process_image() for what cleaned_images/params_df hold.
        raw_images is the pre-cleaning frame for each row of cleaned_images/params_df, aligned
        by position (it's build_calibrations()'s sorted/merged `data`, not a copy - nothing
        extra is read from raw_dir or written to disk for it). The last four arrays are the same
        ones passed to build_calibrations(), returned so callers can build plots straight from
        this source's own frames (e.g. pedestal/pedvar-over-time) without re-reading raw_dir.
    """

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    fallback_map = load_fallback_map(fallback_map_path) if fallback_map_path is not None else None

    results = {}

    groups = group_runs_by_source(raw_dir, telescope, fallback_map=fallback_map)

    for source, runs_by_flip in groups.items():
        source_slug = slugify(source, sep="_")
        source_dir = out_dir / source_slug

        data_preflip, timestamps_preflip = _load_runs(runs_by_flip["preflip"], module_pattern, rate_cut)
        data_postflip, timestamps_postflip = _load_runs(runs_by_flip["postflip"], module_pattern, rate_cut)

        gain_map, gain_source, gain_caption = resolve_gain_map(
            source, data_preflip, timestamps_preflip, data_postflip, timestamps_postflip,
            groups, module_pattern, out_dir, rate_cut,
            nsig=nsig, fit_gaussian=fit_gaussian,
        )
        print(f"  Gain map for {source}: {gain_caption}")

        data, timestamps, pedestals, pedestal_variances, gain = build_calibrations(
            data_preflip, timestamps_preflip, data_postflip, timestamps_postflip,
            source_dir, gain_map, gain_source, gain_caption,
            time_window=time_window, max_gap=max_gap, nsig=nsig, fit_gaussian=fit_gaussian,
        )

        cleaned_images = []
        all_params = []

        for i, (event, pedestal, pedvar, gain_row, timestamp) in enumerate(
                zip(data, pedestals, pedestal_variances, gain, timestamps)):
            cleaned, params = process_image(
                event, pedestal, pedvar, gain_row,
                x=x, y=y,
                image_threshold=image_threshold,
                border_threshold=border_threshold,
            )

            cleaned_images.append(cleaned)

            params['Event'] = i
            params['Timestamp'] = timestamp
            all_params.append(params)

        cleaned_images = np.stack(cleaned_images, axis=0).astype(np.float32)
        params_df = pd.DataFrame(all_params)

        np.savez(
            source_dir / f"{source_slug}.npz",
            cleaned_images=cleaned_images,
            **{col: params_df[col].values for col in params_df.columns},
        )

        results[source] = (
            cleaned_images, params_df, data,
            data_preflip, timestamps_preflip, data_postflip, timestamps_postflip,
        )

    return results
