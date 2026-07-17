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

Once every run is processed this writes manifest.json (per run: telescopes
present and plot path/type/telescopes for each plot) and index.html (a
static gallery of everything generated, grouped by run, then by source
telescope(s), then by plot type).

Usage:
    python run_nextday_pipeline.py --date 20260116
    python run_nextday_pipeline.py --date 20260116 --config my_config.yaml
    python run_nextday_pipeline.py --date 20260116 --raw-dir /data/incoming/20260116 --out-dir /srv/www/heap/20260116
"""
import argparse
import html
import json
import re
import sys
import matplotlib.pyplot as plt
import yaml
from datetime import datetime, timezone
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CONFIG = Path(__file__).resolve().parent / "nextday_config.yaml"
sys.path.insert(0, str(PROJECT_ROOT))

from heap import coincidences as coinc


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


def save_plot(path):
    """Saves the current matplotlib figure to path (creating parent dirs) and closes it.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(path, dpi=300)
    plt.close()


def plot_entry(path, plot_type, telescopes, date_out_dir):
    """Builds a manifest record for a plot file if it was actually written.
    Path is stored relative to the night's out_dir (not the run's), so the gallery
    can link to it as <run_id>/....
    """
    if not path.exists():
        return None
    return {
        "path": str(path.relative_to(date_out_dir)),
        "type": plot_type,
        "telescopes": telescopes,
    }


# A run directory can hold multiple data products per module (dp_ph1024, dp_ph256,
# dp_img8, dp_img16, ...); this pipeline's coincidence analysis only wants the
# photon-timing product. Filenames are start_<ts>.<data_product>.bpp_<n>.<module>.seqno_<n>.pff,
# so the product tag always precedes the module tag.
DATA_PRODUCT = "dp_ph1024"


def load_all_telescopes(telescope_map, raw_dir, run_out_dir, date_out_dir, manifest):
    """Loads + cleans data for every configured telescope found in raw_dir (one run).
    spike_cut plots land in run_out_dir/<DATA_PRODUCT>*<module>/; raw data is read from raw_dir.
    """
    telescopes = {}
    for module, name in telescope_map.items():
        module_pattern = f"{DATA_PRODUCT}*{module}"
        if not sorted(raw_dir.glob(f"*{module_pattern}*.pff")):
            print(f"Skipping {name} ({module}): no {DATA_PRODUCT} .pff files found in {raw_dir}")
            continue
        data, timestamps = coinc.load_telescope_tv(module_pattern, raw_dir, plotting=True)
        telescopes[name] = (data, timestamps)
        spike_cut_path = run_out_dir / module_pattern / "spike_cut.png"
        save_plot(spike_cut_path)
        entry = plot_entry(spike_cut_path, "spike_cut", [name], date_out_dir)
        if entry:
            manifest["plots"].append(entry)
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
        entry = plot_entry(dt_corr_path, "timing_offset_correction", [reference, name], date_out_dir)
        if entry:
            manifest["plots"].append(entry)

        time_coinc, _, _, _ = coinc.match_coinc(timestamps_corr, data, ref_timestamps, ref_data)
        if len(time_coinc):
            coinc.coinc_rate(time_coinc, pair_name, pair_dir, bin_width=300, plotting=True)
            coinc_rate_path = pair_dir / f"{pair_name}_coinc_rate.png"
            save_plot(coinc_rate_path)
            entry = plot_entry(coinc_rate_path, "coincidence_rate", [reference, name], date_out_dir)
            if entry:
                manifest["plots"].append(entry)

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
                entry = plot_entry(triple_coinc_rate_path, "triple_coincidence_rate", [reference, a, b], date_out_dir)
                if entry:
                    manifest["plots"].append(entry)


def process_run(run_dir, telescope_map, reference, date_out_dir, run_manifest):
    """Runs the full single-run pipeline (load, pairs, triples) into date_out_dir/<run_id>/."""
    run_out_dir = date_out_dir / run_manifest["run_id"]
    run_out_dir.mkdir(parents=True, exist_ok=True)

    telescopes = load_all_telescopes(telescope_map, run_dir, run_out_dir, date_out_dir, run_manifest)
    run_manifest["telescopes"] = sorted(telescopes)

    corrected = process_pairs(telescopes, reference, run_out_dir, date_out_dir, run_manifest)
    process_triples(telescopes, corrected, reference, run_out_dir, date_out_dir, run_manifest)


PLOT_TYPE_LABELS = {
    "spike_cut": "Spike cut (single telescope)",
    "timing_offset_correction": "Timing offset correction",
    "coincidence_rate": "Pairwise coincidence rate",
    "triple_coincidence_rate": "Triple coincidence rate",
}


def slugify(text):
    return re.sub(r"[^A-Za-z0-9]+", "-", text).strip("-")


def write_gallery(manifest, out_dir):
    nav_items = []
    sections = []
    all_telescopes = set()

    for run in manifest["runs"]:
        run_anchor = "run-" + slugify(run["run_id"])
        all_telescopes.update(run["telescopes"])
        nav_items.append(f'<a href="#{run_anchor}">{html.escape(run["run_id"])}</a>')
        sections.append(f'<h2 id="{run_anchor}">Run {html.escape(run["run_id"])}</h2>')

        if run.get("error"):
            sections.append(f'<p class="error">Failed: {html.escape(run["error"])}</p>')

        if not run["plots"]:
            sections.append('<p class="empty">No plots for this run.</p>')
            continue

        # Group by source (the telescope or telescope-combo a plot belongs to), preserving
        # the order plots were generated in (singles, then pairs, then triples).
        sources = {}
        for plot in run["plots"]:
            sources.setdefault(tuple(plot["telescopes"]), []).append(plot)

        for source_key, plots in sources.items():
            source_label = " + ".join(source_key)
            sections.append(f'<h3>{html.escape(source_label)}</h3>')

            by_type = {}
            for plot in plots:
                by_type.setdefault(plot["type"], []).append(plot)
            for plot_type, type_plots in by_type.items():
                label = PLOT_TYPE_LABELS.get(plot_type, plot_type)
                sections.append(f'<h4>{html.escape(label)}</h4><div class="grid">')
                for plot in type_plots:
                    src = html.escape(plot["path"])
                    sections.append(f'<figure><img src="{src}" loading="lazy"></figure>')
                sections.append("</div>")

    page = f"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>HEAP night report — {html.escape(manifest['date'])}</title>
<style>
  :root {{ color-scheme: light dark; }}
  body {{ font-family: system-ui, sans-serif; margin: 2rem; max-width: 1200px; }}
  h1 {{ margin-bottom: 0.2rem; }}
  .meta {{ color: #888; margin-bottom: 2rem; }}
  nav {{ display: flex; flex-wrap: wrap; gap: 0.6rem; margin-bottom: 2rem; }}
  nav a {{
    border: 1px solid #8884; border-radius: 999px; padding: 0.25rem 0.8rem;
    font-size: 0.9rem; text-decoration: none; color: inherit;
  }}
  h2 {{ margin-top: 2.5rem; border-bottom: 1px solid #8884; padding-bottom: 0.3rem; }}
  h3 {{ margin-top: 1.5rem; font-size: 1.1rem; }}
  h4 {{ margin-top: 1rem; font-size: 0.9rem; color: #888; font-weight: 600; }}
  .empty {{ color: #888; font-style: italic; }}
  .error {{ color: #c0392b; }}
  .grid {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(340px, 1fr)); gap: 1.2rem; }}
  figure {{ margin: 0; }}
  figure img {{ width: 100%; border-radius: 6px; border: 1px solid #8884; }}
</style>
</head>
<body>
<h1>HEAP night report — {html.escape(manifest['date'])}</h1>
<div class="meta">Reference telescope: {html.escape(manifest['reference_telescope'])} &middot;
Telescopes: {html.escape(', '.join(sorted(all_telescopes)))} &middot;
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

    total_plots = 0
    failed_runs = 0
    for run_dir in runs:
        run_manifest = {"run_id": run_dir.name, "telescopes": [], "plots": []}
        try:
            process_run(run_dir, telescope_map, reference, out_dir, run_manifest)
        except Exception as e:
            print(f"Run {run_manifest['run_id']} failed: {e}")
            run_manifest["error"] = str(e)
            failed_runs += 1
        manifest["runs"].append(run_manifest)
        total_plots += len(run_manifest["plots"])

    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    write_gallery(manifest, out_dir)

    ok_runs = len(runs) - failed_runs
    print(f"\nWrote {total_plots} plots across {ok_runs}/{len(runs)} run(s) (failed_runs={failed_runs}), manifest.json and index.html to {out_dir}")


if __name__ == "__main__":
    main()
