"""
Runs heap analysis tools for one night of data and packages the
resulting plots into a browsable directory structure plus a manifest.json
that a database ingestion step (or anything else) can read.

A night's raw data is split into runs (separate acquisition sessions,
e.g. either side of a meridian flip or a DAQ restart), each in its own
subfolder named after the run's start time. Reads raw .pff files from
<raw_data_dir>/<date>/<run_id>/ and writes plots + manifest.json +
index.html to <out_data_dir>/<date>/, with each run's plots under
<run_id>/ (raw_data_dir and output_dir come from a config file, see
nextday_config.yaml, and can point at the same tree or different ones).

For every run within the night this:
  1. Loads + cleans every telescope's data (heap.Coincidences.load_telescope_tv),
     which also writes each telescope's spike_cut.png into <run_id>/module_<n>/.
  2. Corrects every other telescope's timestamps against a reference
     telescope and computes pairwise coincidence rates, saving plots under
     <run_id>/pairs/<ref>-<name>/.
  3. Computes triple coincidence rates for every pair of non-reference
     telescopes against the reference, saving plots under
     <run_id>/triples/<ref>-<a>-<b>/.

Once every run is loaded, process_telescope_datasets() builds each telescope's
pedestal/pedvar/gain calibrations (grouped by source across the whole night, see
heap.process_dataset.build_calibrations) and 30-min-interval pedestal/pedvar heatmap
grids per flip side (see heap.make_pedestals.plot_pedestal_pedvar_intervals), and
build_hillas_plots() collects every telescope's Hillas parameters (length/width/
log10(size)/distance) into a per-telescope histogram plus one pooling every telescope
together. add_calibration_plots() reuses those per-source/per-telescope files to add
pedestal/pedvar/gain/Hillas plots for every telescope back into each run it belongs to
(rather than recomputing them per run).

Once every run is processed this writes manifest.json (per run: telescopes
present and plot path/type/telescopes for each plot) and index.html (a
static gallery of everything generated, grouped by run, then by source
telescope(s)/plot type as a table).

Usage:
    python run_nextday_pipeline.py --date 20260116
    python run_nextday_pipeline.py --date 20260116 --config my_config.yaml
    python run_nextday_pipeline.py --date 20260116 --raw-dir /data/incoming/20260116 --out-dir /srv/www/heap/20260116
"""
import argparse
import html
import json
import sys
import matplotlib.pyplot as plt
plt.switch_backend("Agg")
import matplotlib.dates as mdates
import numpy as np
import pandas as pd
import yaml
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CONFIG = Path(__file__).resolve().parent / "nextday_config.yaml"
sys.path.insert(0, str(PROJECT_ROOT))

from heap import coincidences as coinc
from heap.make_pedestals import plot_pedestal_mean_over_interval, plot_pedestal_pedvar_intervals, plot_pedvar_histogram
from heap.process_dataset import identify_source, load_fallback_map, process_dataset, slugify


# A run directory can hold multiple data products per module (dp_ph1024, dp_ph256,
# dp_img8, dp_img16, ...); this pipeline's coincidence analysis only wants the
# photon-timing product. Filenames are start_<ts>.<data_product>.bpp_<n>.<module>.seqno_<n>.pff,
# so the product tag always precedes the module tag.
DATA_PRODUCT = "dp_ph1024"

TAB10_COLORS = plt.get_cmap("tab10").colors

PLOT_TYPE_LABELS = {
    "spike_cut": "Event rate and spike cut",
    "coincidence_rate": "Coincidence rates",
    "timing_offset_correction": "Timing offset correction",
    "pedestal": "Pedestal",
    "pedestal_interval": "Pedestal (30-min intervals)",
    "pedestal_mean_over_interval": "Pedestal mean over interval",
    "pedvar": "Pedestal variance",
    "pedvar_interval": "Pedestal variance (30-min intervals)",
    "pedvar_histogram": "Pedestal variance histogram",
    "gain": "Gain map",
    "hillas_params": "Hillas params",
}

# Plot types rendered fresh per run. Everything else in PLOT_TYPE_LABELS is rendered once per
# source across the whole night (see process_telescope_datasets()/build_hillas_plots()). 
PER_RUN_PLOT_TYPES = {"spike_cut", "timing_offset_correction", "coincidence_rate", "triple_coincidence_rate"}


def resolve(path_str):
    path = Path(path_str)
    return path if path.is_absolute() else (PROJECT_ROOT / path)


def load_config(config_path):
    with open(config_path) as f:
        return yaml.safe_load(f)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--date", required=True, help="Night to process, YYYYMMDD.")
    parser.add_argument("--config", default=DEFAULT_CONFIG, type=Path, help=f"Path to config YAML (default: {DEFAULT_CONFIG}).")
    parser.add_argument("--raw-dir", help="Overrides <raw_data_dir>/<date> from the config.")
    parser.add_argument("--out-dir", help="Overrides <output_dir>/<date> from the config.")
    parser.add_argument("--reference", help="Overrides the reference telescope from the config.")
    return parser.parse_args()


def discover_runs(raw_dir):
    """A night's raw_dir contains one subfolder per run, named after the run's start
    time; that name is unique and used as-is for the run_id. Returns them sorted
    chronologically (their names sort lexicographically since they're timestamps).
    """
    return sorted(p for p in raw_dir.iterdir() if p.is_dir())


def telescope_color_map(telescope_map):
    """Assigns each configured telescope a stable color from the tab10 colormap (by name,
    alphabetically), so a given telescope is always the same color across the Hillas plots and
    the combined event-rate plot, regardless of which telescopes are actually present in a given
    run or dataset. Adding/removing telescopes from the config will shift these assignments
    around, but a given config produces the same colors every run.
    """
    names = sorted(info["name"] for info in telescope_map.values())
    return {name: TAB10_COLORS[i % len(TAB10_COLORS)] for i, name in enumerate(names)}


def save_plot(path):
    """Saves the current matplotlib figure to path (creating parent dirs) and closes it.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(path, dpi=300)
    plt.close()


def save_fig(fig, path):
    """Makes fig the current figure and saves it to path (see save_plot())."""
    plt.figure(fig.number)
    save_plot(path)


def save_pedvar_histogram(pedvar, path, color="tab:blue"):
    """Saves a histogram of the 1024 per-pixel pedvar values (e.g. the mean-over-integral array)
    to path, colored per color (default = "tab:blue"; see heap.make_pedestals.plot_pedvar_histogram)."""
    save_fig(plot_pedvar_histogram(pedvar, color=color), path)


def save_calibration_heatmap(values, colorbar_label, path):
    """Saves a 32x32 per-pixel calibration array (pedestal or pedvar) as a heatmap to path."""
    plt.figure(figsize=(4, 3.5))
    im = plt.imshow(values.reshape(32, 32), origin="lower", cmap="viridis")
    plt.colorbar(im, label=colorbar_label)
    plt.xlabel("X")
    plt.ylabel("Y")
    save_plot(path)


def save_pedestal_pedvar_interval_plots(data, timestamps, telescope_dir):
    """Saves 30-min-interval pedestal/pedvar heatmap grids for a telescope's whole night (both
    flip sides pooled together - the pointing changes but the interval grid is about calibration
    drift over time, not flip side) to <telescope_dir>/pedestal_interval.png and
    pedvar_interval.png (see heap.make_pedestals.plot_pedestal_pedvar_intervals). Does nothing if
    data is empty."""
    pedestal_fig, pedvar_fig = plot_pedestal_pedvar_intervals(data, timestamps)
    if pedestal_fig is None:
        return

    save_fig(pedestal_fig, telescope_dir / "pedestal_interval.png")
    save_fig(pedvar_fig, telescope_dir / "pedvar_interval.png")


def save_pedestal_mean_over_interval_plot(data, timestamps, path):
    """Saves the pedestal-mean-over-time plot for a sparse grid of representative pixels to path
    (see heap.make_pedestals.plot_pedestal_mean_over_interval). Does nothing if data is empty."""
    fig = plot_pedestal_mean_over_interval(data, timestamps)
    if fig is None:
        return
    save_fig(fig, path)


def add_plot(plots, path, plot_type, telescopes, date_out_dir, caption=None):
    """Appends a manifest record for a plot file to plots, if it was actually written.
    Path is stored relative to the night's out_dir (not the run's), so the gallery
    can link to it as <run_id>/.... caption is optional extra text shown under the
    plot in the gallery (e.g. to flag a gain map resolved via a fallback tier).
    """
    if not path.exists():
        return
    entry = {
        "path": str(path.relative_to(date_out_dir)),
        "type": plot_type,
        "telescopes": telescopes,
    }
    if caption:
        entry["caption"] = caption
    plots.append(entry)


def plot_combined_event_rate(telescopes, colors, bin_width=30):
    """Overlays each telescope's post-spike-cut trigger rate (re-binned the same way as the
    "After Cut" panel of heap.pre_cleaning.spike_cut) on one set of axes, so a run's telescopes
    can be compared directly. colors: {telescope_name: color}, shared with the Hillas plots (see
    telescope_color_map()).
    """
    fig, ax = plt.subplots(figsize=(10, 4))
    for name, (_, timestamps) in telescopes.items():
        if len(timestamps) == 0:
            continue
        bins = np.arange(timestamps.min(), timestamps.max() + bin_width, bin_width)
        counts, _ = np.histogram(timestamps, bins=bins)
        rate = counts / bin_width
        bin_centers = bins[:-1] + bin_width / 2
        time = pd.to_datetime(bin_centers, unit="s", utc=True).tz_convert("America/Los_Angeles")
        ax.step(time, rate, where="mid", label=name, color=colors.get(name))

    ax.set_title("Event rate (all telescopes)")
    ax.set_xlabel("Time")
    ax.set_ylabel("Trigger Rate [Hz]")
    ax.legend(loc="upper right")
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))
    fig.tight_layout()
    return fig


def load_all_telescopes(telescope_map, raw_dir, run_out_dir, date_out_dir, manifest):
    """Loads + cleans data for every configured telescope found in raw_dir (one run).
    spike_cut plots land in run_out_dir/<DATA_PRODUCT>*<module>/; raw data is read from raw_dir.
    """
    telescopes = {}
    for module, info in telescope_map.items():
        name = info["name"]
        module_pattern = f"{DATA_PRODUCT}*{module}"
        if not sorted(raw_dir.glob(f"*{module_pattern}*.pff")):
            print(f"Skipping {name} ({module}): no {DATA_PRODUCT} .pff files found in {raw_dir}")
            continue
        data, timestamps = coinc.load_telescope_tv(module_pattern, raw_dir, info.get("rate_cut", 20), plotting=True)
        telescopes[name] = (data, timestamps)
        spike_cut_path = run_out_dir / module_pattern / "spike_cut.png"
        save_plot(spike_cut_path)
        add_plot(manifest["plots"], spike_cut_path, "spike_cut", [name], date_out_dir)
    return telescopes


def process_pairs(telescopes, reference, run_out_dir, date_out_dir, manifest):
    """Corrects + matches every other telescope against `reference`. Returns corrected timestamps by name."""
    if reference not in telescopes:
        print(f"Reference telescope '{reference}' not present in this run's data; skipping pair/triple coincidences")
        return {}

    ref_data, ref_timestamps = telescopes[reference]
    corrected = {reference: ref_timestamps}

    for name, (data, timestamps) in telescopes.items():
        if name == reference:
            continue
        pair_name = f"{reference}-{name}"
        pair_dir = run_out_dir / "pairs" / pair_name
        pair_dir.mkdir(parents=True, exist_ok=True)

        timestamps_corr = coinc.correct_time(timestamps, ref_timestamps, plot_name=pair_name, base_dir=pair_dir, plotting=True)
        corrected[name] = timestamps_corr
        dt_corr_path = pair_dir / f"{pair_name}_dt_corr.png"
        save_plot(dt_corr_path)
        add_plot(manifest["plots"], dt_corr_path, "timing_offset_correction", [reference, name], date_out_dir)

        time_coinc, _, _, _ = coinc.match_coinc(timestamps_corr, data, ref_timestamps, ref_data)
        if len(time_coinc):
            coinc.coinc_rate(time_coinc, pair_name, pair_dir, bin_width=300, plotting=True)
            coinc_rate_path = pair_dir / f"{pair_name}_coinc_rate.png"
            save_plot(coinc_rate_path)
            add_plot(manifest["plots"], coinc_rate_path, "coincidence_rate", [reference, name], date_out_dir)

    return corrected


def process_triples(telescopes, corrected, reference, run_out_dir, date_out_dir, manifest):
    """Runs 3-way coincidence for every pair of non-reference telescopes against the reference."""
    others = [name for name in telescopes if name != reference and name in corrected]

    for i in range(len(others)):
        for j in range(i + 1, len(others)):
            a, b = others[i], others[j]
            triple_name = f"{reference}-{a}-{b}"
            triple_dir = run_out_dir / "triples" / triple_name
            triple_dir.mkdir(parents=True, exist_ok=True)

            try:
                time_coinc_ref, _, _, _, _, _ = coinc.match_coinc3(
                    corrected[reference], telescopes[reference][0],
                    corrected[a], telescopes[a][0],
                    corrected[b], telescopes[b][0],
                )
            except NameError:
                # heap.Coincidences.match_coinc3 hits a bug (undefined `tz`) when there
                # are zero triple coincidences; treat that the same as "no matches".
                print(f"No triple coincidences for {triple_name}")
                continue

            if len(time_coinc_ref):
                coinc.coinc_rate(time_coinc_ref, triple_name, triple_dir, bin_width=600, plotting=True)
                triple_coinc_rate_path = triple_dir / f"{triple_name}_coinc_rate.png"
                save_plot(triple_coinc_rate_path)
                add_plot(manifest["plots"], triple_coinc_rate_path, "triple_coincidence_rate", [reference, a, b], date_out_dir)


def process_run(run_dir, telescope_map, reference, date_out_dir, run_manifest, colors):
    """Runs the full single-run pipeline (load, pairs, triples, combined event rate) into
    date_out_dir/<run_id>/. colors: {telescope_name: color}, shared with the Hillas plots (see
    telescope_color_map())."""
    run_out_dir = date_out_dir / run_manifest["run_id"]
    run_out_dir.mkdir(parents=True, exist_ok=True)

    telescopes = load_all_telescopes(telescope_map, run_dir, run_out_dir, date_out_dir, run_manifest)
    run_manifest["telescopes"] = list(telescopes)

    corrected = process_pairs(telescopes, reference, run_out_dir, date_out_dir, run_manifest)
    process_triples(telescopes, corrected, reference, run_out_dir, date_out_dir, run_manifest)

    if len(telescopes) > 1:
        fig = plot_combined_event_rate(telescopes, colors)
        combined_path = run_out_dir / "event_rate_combined.png"
        save_fig(fig, combined_path)
        add_plot(run_manifest["plots"], combined_path, "spike_cut", sorted(telescopes), date_out_dir)


def _sorted_by_time(data_parts, timestamps_parts):
    """Concatenates chunks of (data, timestamps) - one chunk per source - and sorts the result by
    timestamp. results.items() (see process_telescope_datasets()) groups process_dataset()'s
    output by source rather than by time, so pooling multiple sources' frames for one telescope
    needs an explicit re-sort before any interval-based plot, which assumes chronological order."""
    if not data_parts:
        return np.empty((0, 1024)), np.empty((0,))
    data = np.concatenate(data_parts, axis=0)
    timestamps = np.concatenate(timestamps_parts, axis=0)
    sort_idx = np.argsort(timestamps)
    return data[sort_idx], timestamps[sort_idx]


def process_telescope_datasets(telescope_map, raw_dir, date_out_dir, colors, image_threshold, border_threshold, event_preview_min_pixels):
    """Runs process_dataset() once per configured telescope for the whole night: builds
    pedestal/pedvar/gain calibrations, cleans every event's image, and parameterizes it into
    Hillas parameters (see heap.process_dataset.process_dataset()). It groups the night's runs
    by source itself, unlike process_run() below which works one run at a time. Writes
    <date_out_dir>/<name>/<source_slug>/<source_slug>.npz plus that same folder's
    calibrations.npz (source_slug is source with spaces/other non-alphanumerics replaced by
    underscores; see heap.process_dataset.slugify()).

    Also renders each source's pedestal/pedvar/gain heatmaps and pedvar histogram here (once,
    since calibrations.npz already covers the whole night) rather than leaving add_calibration_plots()
    to redraw them from scratch for every run that shares the source. The pedestal/pedvar interval
    grids and the pedestal-mean-over-interval plot aren't source- or flip-side-specific, so those
    pool every source's and flip side's frames into one whole-night plot per telescope instead (at
    <date_out_dir>/<name>/, not nested under a source_slug).

    colors: {telescope_name: color}, used to color each telescope's pedvar histogram (see
        telescope_color_map()).

    image_threshold, border_threshold: threshold_clean()'s per-pixel cut, in units of pedvar (see
        heap.image_cleaning.threshold_clean, called from heap.process_dataset.process_image());
        from nextday_config.yaml's pipeline section.

    event_preview_min_pixels: an event needs more than this many surviving (post-cleaning) pixels
        to be included in the returned telescope_events.

    Looks for <raw_dir>/source_run_map.json as a fallback map for identify_source()/
    identify_flip_side(), used when hk.pff's mount tracking is missing or ambiguous (see
    process_dataset.load_fallback_map()).

    Returns {telescope_name: (timestamps, cleaned_images, raw_images)}, one entry per telescope
    with usable data this night, each array concatenated across every source in chronological
    order - build_event_preview() samples coincident events from these. Held in memory only
    (process_dataset() doesn't write raw_images to disk); nothing extra is read from raw_dir or
    saved for this.
    """
    fallback_map_path = raw_dir / "source_run_map.json"
    if not fallback_map_path.exists():
        fallback_map_path = None

    telescope_events = {}

    for module, info in telescope_map.items():
        name = info["name"]
        module_pattern = f"{DATA_PRODUCT}*{module}"
        if not sorted(raw_dir.glob(f"*/*{module_pattern}*.pff")):
            print(f"Skipping {name} ({module}) dataset processing: no {DATA_PRODUCT} .pff files found in {raw_dir}")
            continue

        try:
            results = process_dataset(
                raw_dir, date_out_dir / name, module_pattern, name,
                fallback_map_path=fallback_map_path, rate_cut=info.get("rate_cut", 20),
                image_threshold=image_threshold, border_threshold=border_threshold,
            )
        except Exception as e:
            print(f"Dataset processing for {name} failed: {e}")
            continue

        timestamps_parts, cleaned_parts, raw_parts = [], [], []
        combined_data_parts, combined_ts_parts = [], []

        for source, (cleaned_images, params_df, raw_images, data_preflip, timestamps_preflip, data_postflip, timestamps_postflip) in results.items():
            source_dir = date_out_dir / name / slugify(source, sep="_")

            calib = np.load(source_dir / "calibrations.npz")
            mean_pedvar = calib["pedestal_variances"].mean(axis=0)
            save_calibration_heatmap(calib["pedestals"].mean(axis=0), "Mean", source_dir / "pedestal.png")
            save_calibration_heatmap(mean_pedvar, r"$\sigma$", source_dir / "pedvar.png")
            save_pedvar_histogram(mean_pedvar, source_dir / "pedvar_hist.png", color=colors.get(name, "tab:blue"))
            save_calibration_heatmap(calib["gain"], "Relative gain", source_dir / "gain.png")

            above_pixel_cut = params_df["N_pix"].values > event_preview_min_pixels
            timestamps_parts.append(params_df["Timestamp"].values[above_pixel_cut])
            cleaned_parts.append(cleaned_images[above_pixel_cut])
            raw_parts.append(raw_images[above_pixel_cut])

            combined_data_parts += [data_preflip, data_postflip]
            combined_ts_parts += [timestamps_preflip, timestamps_postflip]

        print(f"{name}: wrote {len(results)} source(s) to {date_out_dir / name}: {', '.join(sorted(results))}")

        # Pedestal/pedvar interval grids and the mean-over-interval plot aren't source- or
        # flip-side-specific (unlike pedestal.png/pedvar.png/gain.png above, which are one
        # calibration per source) - pool every source's and flip side's frames for this telescope
        # into one whole-night pair of plots, sorted back into chronological order since
        # results.items() groups by source rather than by time.
        telescope_dir = date_out_dir / name
        data_all, timestamps_all = _sorted_by_time(combined_data_parts, combined_ts_parts)
        save_pedestal_pedvar_interval_plots(data_all, timestamps_all, telescope_dir)
        save_pedestal_mean_over_interval_plot(data_all, timestamps_all, telescope_dir / "pedestal_mean_over_interval.png")

        if timestamps_parts:
            telescope_events[name] = (
                np.concatenate(timestamps_parts),
                np.concatenate(cleaned_parts, axis=0),
                np.concatenate(raw_parts, axis=0),
            )

    return telescope_events


def add_calibration_plots(run_dir, run_manifest, date_out_dir, fallback_map, hillas_telescopes=()):
    """Adds pedestal/pedvar/gain heatmaps for every telescope present in this run, plus each
    telescope's own Hillas params histogram and, once, a column pooling every telescope in
    hillas_telescopes together.

    Reuses the per-source calibrations.npz, pedestal/pedvar/gain heatmaps, and
    per-telescope hillas_params.png already written by
    process_telescope_datasets()/build_hillas_plots() (built once across the whole night) instead
    of recomputing them from this run's frames - must run after those so those files exist.
    """
    for name in run_manifest["telescopes"]:
        try:
            source = identify_source(run_dir, name, fallback_map=fallback_map)
        except ValueError as e:
            print(f"Skipping calibration plots for {name} in {run_dir.name}: {e}")
            continue

        source_dir = date_out_dir / name / slugify(source, sep="_")
        calib_path = source_dir / "calibrations.npz"
        if not calib_path.exists():
            continue

        calib = np.load(calib_path)

        add_plot(run_manifest["plots"], source_dir / "pedestal.png", "pedestal", [name], date_out_dir, caption=source)
        add_plot(run_manifest["plots"], source_dir / "pedvar.png", "pedvar", [name], date_out_dir, caption=source)
        add_plot(run_manifest["plots"], source_dir / "pedvar_hist.png", "pedvar_histogram", [name], date_out_dir, caption=source)

        gain_caption = calib["gain_caption"].item() if "gain_caption" in calib else None
        add_plot(run_manifest["plots"], source_dir / "gain.png", "gain", [name], date_out_dir, caption=gain_caption)

        add_plot(run_manifest["plots"], date_out_dir / name / "hillas_params.png", "hillas_params", [name], date_out_dir)

        add_plot(
            run_manifest["plots"], date_out_dir / name / "pedestal_mean_over_interval.png",
            "pedestal_mean_over_interval", [name], date_out_dir,
        )
        add_plot(run_manifest["plots"], date_out_dir / name / "pedestal_interval.png", "pedestal_interval", [name], date_out_dir)
        add_plot(run_manifest["plots"], date_out_dir / name / "pedvar_interval.png", "pedvar_interval", [name], date_out_dir)

    if hillas_telescopes:
        add_plot(
            run_manifest["plots"], date_out_dir / "hillas_params_combined.png", "hillas_params",
            list(hillas_telescopes), date_out_dir,
        )


def find_coincidences(telescope_events, window=0.001):
    """Finds every 2-way coincidence between every pair of telescopes in telescope_events (on
    raw, uncorrected timestamps - see heap.coincidences.match_coinc), then unions overlapping
    matches into groups via union-find, so a group spans however many telescopes ended up linked
    together (2 or more).

    telescope_events: {telescope_name: (timestamps, cleaned_images, raw_images)}, see
        process_telescope_datasets().

    Returns a list of (telescope_names, event_indices, event_times) tuples, one per coincidence
    group, with telescope_names/event_indices/event_times all ordered the same way (by
    telescope_events' iteration order).
    """
    names = list(telescope_events)
    parent = {}

    def find(node):
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for name, (timestamps, _, _) in telescope_events.items():
        for idx in range(len(timestamps)):
            parent[(name, idx)] = (name, idx)

    node_time = {}
    for a, b in combinations(names, 2):
        timestamps_a = telescope_events[a][0]
        timestamps_b = telescope_events[b][0]
        t1, i1, t2, i2 = coinc.match_coinc(
            timestamps_a, np.arange(len(timestamps_a)), timestamps_b, np.arange(len(timestamps_b)), window=window,
        )
        for ta, ia, tb, ib in zip(t1, i1, t2, i2):
            node_a, node_b = (a, ia), (b, ib)
            node_time[node_a] = ta
            node_time[node_b] = tb
            union(node_a, node_b)

    groups = {}
    for node in node_time:
        groups.setdefault(find(node), []).append(node)

    coincidences = []
    for nodes in groups.values():
        if len({n[0] for n in nodes}) < 2:
            continue
        nodes = sorted(nodes, key=lambda n: names.index(n[0]))
        coincidences.append(([n[0] for n in nodes], [n[1] for n in nodes], [node_time[n] for n in nodes]))

    return coincidences


def plot_preview_event(tel_names, event_idx, telescope_events, all_names):
    """One figure for a single coincident event: 2 rows (raw, cleaned) x len(all_names) columns
    (one per telescope with data this night), matching next_day_plots_pff.ipynb's plot_event().
    Telescopes not part of this coincidence get a blank (but bordered) column, so events with
    different numbers of participating telescopes stay directly comparable. Each image is scaled
    to its own min/max (no shared color scale)."""
    fig, axs = plt.subplots(2, len(all_names), figsize=(3.5 * len(all_names), 7), squeeze=False)

    for col, name in enumerate(all_names):
        if name not in tel_names:
            for ax in (axs[0, col], axs[1, col]):
                ax.axis("off")
            continue

        i = tel_names.index(name)
        idx = event_idx[i]
        timestamps, cleaned_images, raw_images = telescope_events[name]
        ts_str = pd.to_datetime(timestamps[idx], unit="s", utc=True).tz_convert("America/Los_Angeles").strftime("%Y-%m-%d %H:%M:%S.%f")

        for row, (label, img) in enumerate([("raw", raw_images[idx].reshape(32, 32)), ("cleaned", cleaned_images[idx])]):
            ax = axs[row, col]
            im = ax.imshow(img, aspect="equal", origin="lower")
            plt.colorbar(im, ax=ax, label="ADC")
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(f"{name} {label} event {idx}\n{ts_str}", fontsize=9)

    fig.tight_layout()
    return fig


def build_event_preview(telescope_events, out_dir, manifest, image_threshold, border_threshold, n_events=12, window=0.001):
    """Finds every cross-telescope coincidence across the whole night (see find_coincidences()),
    randomly samples up to n_events of them, and saves one raw+cleaned image grid per sampled
    event to <out_dir>/event_preview/event_<n>.png (see plot_preview_event()). Adds an
    "event_preview" entry to manifest: a description string plus one {path, telescopes,
    timestamp} record per sampled event, for write_gallery() to render as a bottom-of-page
    section. Does nothing if fewer than 2 telescopes have data, or no coincidences are found.

    image_threshold, border_threshold: only used to quote threshold_clean()'s actual cut in the
        description text - the images themselves were already cleaned with these same values by
        process_telescope_datasets().
    """
    if len(telescope_events) < 2:
        return

    coincidences = find_coincidences(telescope_events, window=window)
    if not coincidences:
        print("Event preview: no cross-telescope coincidences found")
        return

    all_names = list(telescope_events)
    n_sample = min(n_events, len(coincidences))
    sample = np.random.default_rng().choice(len(coincidences), size=n_sample, replace=False)

    events = []
    for n, i in enumerate(sample):
        tel_names, event_idx, event_time = coincidences[i]
        fig = plot_preview_event(tel_names, event_idx, telescope_events, all_names)
        path = out_dir / "event_preview" / f"event_{n}.png"
        save_fig(fig, path)
        ts_str = pd.to_datetime(min(event_time), unit="s", utc=True).tz_convert("America/Los_Angeles").strftime("%Y-%m-%d %H:%M:%S")
        events.append({"path": str(path.relative_to(out_dir)), "telescopes": tel_names, "timestamp": ts_str})

    manifest["event_preview"] = {
        "description": (
            f"{n_sample} events randomly sampled from the {len(coincidences)} coincidences "
            "found across the whole night. Each event is shown as raw (top row) and cleaned (bottom row) images "
            "for each participating telescope. Threshold cleaning is used with image cuts at "
            f"{image_threshold:g}σ and border cuts at {border_threshold:g}σ (see heap.image_cleaning.threshold_clean)."
        ),
        "events": events,
    }


def convert_units(df):
    """Converts a params_df's pixel-unit columns (see heap.parameterize.calc_params()) to
    degrees, using the camera's plate scale (0.31 deg/pixel), and re-centers x_c/y_c on the
    camera's optical center. Returns a copy; df is left untouched."""
    df = df.copy()
    columns = ["x_c", "y_c", "s_xx", "s_yy", "s_xy", "length", "width", "miss", "distance"]
    for c in columns:
        df[c] = df[c] * 0.31
        if c in ["x_c", "y_c"]:
            df[c] = df[c] - 4.96 + 0.31
    return df


def postprocess_df(df):
    """Drops events with no well-defined Hillas ellipse (missing length/width/miss/distance/alpha,
    zero width, or too few surviving pixels) and duplicate rows, ahead of plotting."""
    df = df.dropna(subset=["length", "width", "miss", "distance", "alpha"])
    df = df[df["width"] > 0]
    df = df[df.N_pix > 3]
    df = df.drop_duplicates(subset=df.drop(["Event"], axis=1))
    return df


def plot_hillas_histograms(dfs, title="Hillas Params", colors=None, pooled_df=None, pooled_label="All telescopes"):
    """
    Plots length/width/log10(size)/distance histograms (params_df must already be in degrees,
    see convert_units()), one histtype="step" line per entry in dfs, overlaid on the same 4 axes.
    If pooled_df is given, also overlays it as a black, alpha=0.4 filled histogram (e.g. every
    telescope's events pooled together, alongside each telescope's own step histogram in dfs).

    Parameters:
        dfs: {label: params_df}, e.g. one entry per telescope for a per-telescope comparison.
        title: figure title (default = "Hillas Params")
        colors: optional {label: color}, so a given telescope keeps the same color across calls
            (e.g. the per-telescope plot and the combined plot). Falls back to the tab10
            colormap, by position in dfs, for any label not present.
        pooled_df: optional params_df of every telescope's events pooled together, drawn as a
            filled black histogram alongside the per-entry step histograms in dfs.
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
        axs[3].hist(df.distance, bins=80, range=(0, 6), histtype="step", density=True, label=label, color=color)

    if pooled_df is not None:
        label = f"{pooled_label} N={len(pooled_df)}"
        width_mean = pooled_df.width.mean()
        width_label = f"{label}, $\\mu$={width_mean:.3f}°"
        axs[0].hist(pooled_df.length, bins=80, range=(0, 2), histtype="stepfilled", density=True, label=label, color="black", alpha=0.4)
        axs[1].hist(pooled_df.width, bins=80, range=(0, 1), histtype="stepfilled", density=True, label=width_label, color="black", alpha=0.4)
        axs[2].hist(np.log10(pooled_df["size"]), bins=80, range=(0, 6), histtype="stepfilled", density=True, label=label, color="black", alpha=0.4)
        axs[3].hist(pooled_df.distance, bins=80, range=(0, 6), histtype="stepfilled", density=True, label=label, color="black", alpha=0.4)

    for ax in axs:
        ax.legend(loc="upper right")

    fig.tight_layout()
    return fig


def load_telescope_params_df(telescope_dir):
    """Loads and concatenates every source's params_df for one telescope (see
    heap.process_dataset.process_dataset(), which writes one <source_slug>.npz per source under
    telescope_dir), converted to degrees and postprocessed (see convert_units()/postprocess_df()).
    Returns None if telescope_dir has no usable source data.
    """
    dfs = []
    for npz_path in sorted(telescope_dir.glob("*/*.npz")):
        if npz_path.name == "calibrations.npz":
            continue
        npz = np.load(npz_path)
        dfs.append(pd.DataFrame({col: npz[col] for col in npz.files if col != "cleaned_images"}))

    if not dfs:
        return None

    df = postprocess_df(convert_units(pd.concat(dfs, ignore_index=True)))
    return df if len(df) else None


def build_hillas_plots(telescope_map, out_dir, colors):
    """Builds each telescope's own Hillas parameter histogram (length/width/log10(size)/distance,
    over its whole night of data across every source) at <out_dir>/<name>/hillas_params.png, plus
    one at <out_dir>/hillas_params_combined.png overlaying every telescope's own step histogram
    (colored per colors) alongside every telescope's events pooled into one filled histogram (see
    plot_hillas_histograms()). Reused across every run via add_calibration_plots(), the same way
    pedestal/pedvar/gain are.

    colors: {telescope_name: color}, shared with the combined event-rate plot so a given
        telescope keeps the same color everywhere (see telescope_color_map()).

    Returns the sorted list of telescopes with usable Hillas data (i.e. those making up the
    combined plot), or [] if none had any.
    """
    telescope_dfs = {}
    for info in telescope_map.values():
        name = info["name"]
        df = load_telescope_params_df(out_dir / name)
        if df is None:
            continue
        telescope_dfs[name] = df

        fig = plot_hillas_histograms({name: df}, title=f"Hillas Params — {name}", colors=colors)
        save_fig(fig, out_dir / name / "hillas_params.png")

    if not telescope_dfs:
        return []

    pooled_df = pd.concat(telescope_dfs.values(), ignore_index=True)
    fig = plot_hillas_histograms(
        telescope_dfs, title="Hillas Params (combined)", colors=colors, pooled_df=pooled_df,
    )
    save_fig(fig, out_dir / "hillas_params_combined.png")

    return sorted(telescope_dfs)


def render_plot_table(plots):
    """Renders one <table class="plots"> for a list of manifest plot entries. Columns are
    sources (the telescope or telescope-combo a plot belongs to), in the order first
    encountered; rows are plot types, in PLOT_TYPE_LABELS order."""
    source_columns = []
    present_types = set()
    plots_by_cell = {}
    for plot in plots:
        source_key = tuple(sorted(plot["telescopes"]))
        if source_key not in source_columns:
            source_columns.append(source_key)
        row_type = "coincidence_rate" if plot["type"] == "triple_coincidence_rate" else plot["type"]
        present_types.add(row_type)
        plots_by_cell.setdefault((source_key, row_type), []).append(plot)

    type_rows = [plot_type for plot_type in PLOT_TYPE_LABELS if plot_type in present_types]

    header_cells = "".join(f'<th>{html.escape(" + ".join(s))}</th>' for s in source_columns)
    parts = [f'<div class="table-scroll"><table class="plots"><thead><tr><th>Plot</th>{header_cells}</tr></thead><tbody>']
    for plot_type in type_rows:
        label = html.escape(PLOT_TYPE_LABELS.get(plot_type, plot_type))
        row_cells = []
        for source_key in source_columns:
            cell_plots = plots_by_cell.get((source_key, plot_type), [])
            imgs = "".join(
                f'<a href="{html.escape(p["path"])}" target="_blank" rel="noopener">'
                f'<img src="{html.escape(p["path"])}" loading="lazy"></a>'
                + (f'<div class="caption">{html.escape(p["caption"])}</div>' if p.get("caption") else "")
                for p in cell_plots
            )
            row_cells.append(f"<td>{imgs}</td>")
        parts.append(f'<tr><th scope="row">{label}</th>{"".join(row_cells)}</tr>')
    parts.append("</tbody></table></div>")
    return "".join(parts)


def write_gallery(manifest, out_dir):
    nav_items = []
    sections = []
    all_telescopes = set()

    # Whole-night plots are referenced by every run that shares their source; keep only the
    # first occurrence of each file for the calibration section (telescopes with more than one
    # source this night still get one image per source, since each source has its own file).
    calibration_plots = []
    seen_paths = set()
    for run in manifest["runs"]:
        for plot in run["plots"]:
            if plot["type"] in PER_RUN_PLOT_TYPES:
                continue
            if plot["path"] in seen_paths:
                continue
            seen_paths.add(plot["path"])
            calibration_plots.append(plot)

    if calibration_plots:
        nav_items.append('<a href="#calibration">Calibration</a>')
        sections.append('<h2 id="calibration">Calibration (whole night)</h2>')
        sections.append(render_plot_table(calibration_plots))

    for run in manifest["runs"]:
        run_anchor = "run-" + slugify(run["run_id"])
        all_telescopes.update(run["telescopes"])
        nav_items.append(f'<a href="#{run_anchor}">{html.escape(run["run_id"])}</a>')
        sections.append(f'<h2 id="{run_anchor}">Run {html.escape(run["run_id"])}</h2>')

        if run.get("error"):
            sections.append(f'<p class="error">Failed: {html.escape(run["error"])}</p>')

        run_plots = [p for p in run["plots"] if p["type"] in PER_RUN_PLOT_TYPES]
        if not run_plots:
            sections.append('<p class="empty">No plots for this run.</p>')
            continue

        sections.append(render_plot_table(run_plots))

    if manifest.get("event_preview"):
        preview = manifest["event_preview"]
        nav_items.append('<a href="#event-preview">Event Preview</a>')
        sections.append('<h2 id="event-preview">Event Preview</h2>')
        sections.append(f'<p>{html.escape(preview["description"])}</p>')
        sections.append('<div class="event-grid">')
        for event in preview["events"]:
            caption = html.escape(f'{" + ".join(event["telescopes"])} · {event["timestamp"]}')
            sections.append(
                '<div class="event-card">'
                f'<a href="{html.escape(event["path"])}" target="_blank" rel="noopener">'
                f'<img src="{html.escape(event["path"])}" loading="lazy"></a>'
                f'<div class="caption">{caption}</div>'
                '</div>'
            )
        sections.append('</div>')

    telescope_labels = [
        f"{name} (reference)" if name == manifest["reference_telescope"] else name
        for name in sorted(all_telescopes)
    ]

    page = f"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>HEAP Next Day — {html.escape(manifest['date'])}</title>
<style>
  :root {{ color-scheme: light dark; }}
  body {{ font-family: system-ui, sans-serif; margin: 2rem; }}
  h1 {{ margin-bottom: 0.2rem; }}
  .meta {{ color: #888; margin-bottom: 2rem; }}
  nav {{ display: flex; flex-wrap: wrap; gap: 0.6rem; margin-bottom: 2rem; }}
  nav a {{
    border: 1px solid #8884; border-radius: 999px; padding: 0.25rem 0.8rem;
    font-size: 0.9rem; text-decoration: none; color: inherit;
  }}
  h2 {{ margin-top: 2.5rem; border-bottom: 1px solid #8884; padding-bottom: 0.3rem; }}
  .empty {{ color: #888; font-style: italic; }}
  .error {{ color: #c0392b; }}
  .table-scroll {{ margin-bottom: 1.5rem; }}
  table.plots {{ border-collapse: collapse; }}
  table.plots th, table.plots td {{ text-align: left; padding: 0.5rem 0.8rem; border-bottom: 1px solid #8884; vertical-align: middle; }}
  table.plots thead th {{ color: #888; font-weight: 600; font-size: 0.9rem; white-space: nowrap; }}
  table.plots th[scope="row"] {{ white-space: nowrap; font-size: 0.9rem; }}
  table.plots img {{ display: block; height: 160px; width: 220px; object-fit: contain; border-radius: 6px; border: 1px solid #8884; }}
  table.plots .caption {{ color: #888; font-size: 0.75rem; max-width: 200px; }}
  .event-grid {{ display: flex; flex-wrap: wrap; gap: 1.2rem; margin-bottom: 1.5rem; }}
  .event-card img {{ display: block; max-width: 480px; width: 100%; border-radius: 6px; border: 1px solid #8884; }}
  .event-card .caption {{ color: #888; font-size: 0.8rem; margin-top: 0.3rem; }}
</style>
</head>
<body>
<h1>HEAP Next Day — {html.escape(manifest['date'])}</h1>
<div class="meta">Telescopes: {html.escape(', '.join(telescope_labels))} &middot;
Generated {html.escape(manifest['generated_at'])}</div>
<nav>{''.join(nav_items)}</nav>
{''.join(sections)}
</body>
</html>
"""
    (out_dir / "index.html").write_text(page)


def main():
    args = parse_args()
    config = load_config(args.config)

    telescope_map = config["telescopes"]
    reference = args.reference or config["pipeline"]["reference"]
    image_threshold = config["pipeline"].get("image_threshold", 4.0)
    border_threshold = config["pipeline"].get("border_threshold", 2.0)
    event_preview_min_pixels = config["pipeline"].get("event_preview_min_pixels", 3)
    colors = telescope_color_map(telescope_map)

    raw_dir = Path(args.raw_dir) if args.raw_dir else resolve(config["paths"]["raw_data_dir"]) / args.date
    out_dir = Path(args.out_dir) if args.out_dir else resolve(config["paths"]["output_dir"]) / args.date

    if not raw_dir.exists():
        sys.exit(f"Raw data directory not found: {raw_dir}")
    out_dir.mkdir(parents=True, exist_ok=True)

    runs = discover_runs(raw_dir)
    if not runs:
        sys.exit(f"No run subfolders found in {raw_dir}")

    manifest = {
        "date": args.date,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "reference_telescope": reference,
        "runs": [],
    }

    failed_runs = 0
    for run_dir in runs:
        run_manifest = {"run_id": run_dir.name, "telescopes": [], "plots": []}
        try:
            process_run(run_dir, telescope_map, reference, out_dir, run_manifest, colors)
        except Exception as e:
            print(f"Run {run_manifest['run_id']} failed: {e}")
            run_manifest["error"] = str(e)
            failed_runs += 1
        manifest["runs"].append(run_manifest)

    telescope_events = process_telescope_datasets(
        telescope_map, raw_dir, out_dir, colors, image_threshold, border_threshold, event_preview_min_pixels,
    )
    hillas_telescopes = build_hillas_plots(telescope_map, out_dir, colors)
    build_event_preview(telescope_events, out_dir, manifest, image_threshold, border_threshold)

    fallback_map_path = raw_dir / "source_run_map.json"
    fallback_map = load_fallback_map(fallback_map_path) if fallback_map_path.exists() else None
    for run_dir, run_manifest in zip(runs, manifest["runs"]):
        if run_manifest.get("error"):
            continue
        add_calibration_plots(run_dir, run_manifest, out_dir, fallback_map, hillas_telescopes)

    total_plots = sum(len(run_manifest["plots"]) for run_manifest in manifest["runs"])
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    write_gallery(manifest, out_dir)

    ok_runs = len(runs) - failed_runs
    print(f"\nWrote {total_plots} plots across {ok_runs}/{len(runs)} run(s) (failed_runs={failed_runs}), manifest.json and index.html to {out_dir}")


if __name__ == "__main__":
    main()
