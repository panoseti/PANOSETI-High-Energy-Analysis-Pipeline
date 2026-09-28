"""config

Loads an analysis config (e.g. configs/crab.yaml, for analysis_tools/pff_analysis.ipynb) and
converts its values into the forms heap.events/heap.significance take.
"""
from pathlib import Path

import astropy.units as u
import yaml
from astropy.coordinates import SkyCoord

from heap.significance import OFF_REGION_METHODS

SECTIONS = {
    "paths": {"raw_data_dir", "output_dir"},
    "source": {"name", "dates", "pointing_overrides"},
    "telescopes": None, # module -> {name, rate_cut}
    "pipeline": {"data_product", "image_threshold", "border_threshold", "keep_brightest_island"},
    "events": {"reference", "coinc_window", "rotate_postflip", "rel_tel_efficiency", "pointing_corrections"},
    "cuts": {"min_npix", "min_tel", "telescopes", "theta", "max_distance", "off_regions", "image_cuts", "position_cuts"},
}
CUT_MODES = {"nsigma", "value"}
POSITION_CUT_COLUMNS = {"Distance", "Alpha", "Miss"}
FLIP_SIDES = {"preflip", "postflip"}
PROJECT_ROOT = Path(__file__).resolve().parent.parent


def _check_keys(section, config, required):
    missing = required - set(config)
    unknown = set(config) - required
    if missing or unknown:
        raise ValueError(f"config section {section!r}: missing {sorted(missing)}, unknown {sorted(unknown)}")


def _cuts(name, cuts, columns=None):
    """{column: [mode, threshold]} -> {column: (mode, threshold)}; None -> {}."""
    out = {}
    for column, cut in (cuts or {}).items():
        if not isinstance(cut, (list, tuple)) or len(cut) != 2 or cut[0] not in CUT_MODES:
            raise ValueError(f"cuts.{name}.{column}: expected [mode, threshold] with mode in {sorted(CUT_MODES)}, got {cut!r}")
        if columns is not None and column not in columns:
            raise ValueError(f"cuts.{name}: {column!r} not one of {sorted(columns)}")
        out[column] = (cut[0], float(cut[1]))
    return out


def _pointing(pointing):
    """"RA DEC" (hourangle, deg) string or [ra, dec] (deg) -> SkyCoord; None -> None."""
    if pointing is None:
        return None
    if isinstance(pointing, str):
        return SkyCoord(pointing, unit=(u.hourangle, u.deg))
    ra, dec = pointing
    return SkyCoord(ra, dec, unit=u.deg)


def _resolve(path):
    path = Path(path).expanduser()
    return path if path.is_absolute() else PROJECT_ROOT / path


def load_analysis_config(path):
    """
    Load and check an analysis config.

    Parameters:
        path: config YAML, see configs/crab.yaml

    Returns:
        the config's sections as nested dicts, with
            paths.*: Path (relative paths resolved against the project root; output_dir
                defaults to <raw_data_dir>/processed)
            source.dates: list of "YYYYMMDD" strings, or None
            source.pointing_overrides: {run folder name: SkyCoord}
            events.pointing_corrections: {(date, telescope, flip_side): (dx, dy)}, see
                heap.events.build_night_events()
            cuts.telescopes: list of telescope names (default = every telescope in telescopes)
            cuts.image_cuts, cuts.position_cuts: {column: (mode, threshold)}, see
                heap.events.apply_cuts() and heap.significance.OnOffCounter
    """
    path = Path(path)
    with open(path) as f:
        config = yaml.safe_load(f)

    _check_keys("top level", config, set(SECTIONS))
    for section, keys in SECTIONS.items():
        if keys is not None:
            _check_keys(section, config[section], keys)

    paths = config["paths"]
    paths["raw_data_dir"] = _resolve(paths["raw_data_dir"])
    paths["output_dir"] = _resolve(paths["output_dir"]) if paths["output_dir"] is not None else paths["raw_data_dir"] / "processed"

    source = config["source"]
    source["dates"] = [str(d) for d in source["dates"]] if source["dates"] is not None else None
    source["pointing_overrides"] = {str(run): _pointing(p) for run, p in (source["pointing_overrides"] or {}).items()}

    names = [info["name"] for info in config["telescopes"].values()]
    events = config["events"]
    if events["reference"] not in names:
        raise ValueError(f"events.reference: {events['reference']!r} not one of {names}")
    events["rel_tel_efficiency"] = events["rel_tel_efficiency"] or {}
    corrections = {}
    for date, telescopes in (events["pointing_corrections"] or {}).items():
        for telescope, sides in telescopes.items():
            for side, (dx, dy) in sides.items():
                if telescope not in names or side not in FLIP_SIDES:
                    raise ValueError(f"events.pointing_corrections.{date}.{telescope}.{side}: unknown telescope or flip side")
                corrections[(str(date), telescope, side)] = (float(dx), float(dy))
    events["pointing_corrections"] = corrections

    cuts = config["cuts"]
    cuts["telescopes"] = cuts["telescopes"] if cuts["telescopes"] is not None else names
    if cuts["off_regions"] not in OFF_REGION_METHODS:
        raise ValueError(f"cuts.off_regions: {cuts['off_regions']!r} not one of {sorted(OFF_REGION_METHODS)}")
    cuts["image_cuts"] = _cuts("image_cuts", cuts["image_cuts"])
    cuts["position_cuts"] = _cuts("position_cuts", cuts["position_cuts"], POSITION_CUT_COLUMNS)

    return config
