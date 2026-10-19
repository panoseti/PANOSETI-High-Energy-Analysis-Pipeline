"""sources

Catalog of observed sources, used to strip the wobble from a run's mount target_name (e.g.
"Crab 0.5 N" -> "Crab", see parse_target_name()), to identify a run's target from its mount pointing when
target_name is blank (see heap.process_dataset.identify_source()), and to measure each run's wobble
offset from its pointing (see wobble_offset()). Names are the source names analyses use (e.g.
pff_analysis.ipynb's SOURCE); add an entry here for each newly observed source.
"""
import re
import warnings

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, search_around_sky

SOURCES = {
    "Crab": SkyCoord("05 34 31.78 +22 01 02.6", unit=(u.hourangle, u.deg)),
    "MGRO J2019+37": SkyCoord(304.67, 36.69, unit=u.deg), # 4HWC J2018+3641 (Alfaro et al. 2026)
    "NGC 1275": SkyCoord("03 19 48.16 +41 30 42.1", unit=(u.hourangle, u.deg)),
    "Mrk 421": SkyCoord("11 04 27.31 +38 12 31.8", unit=(u.hourangle, u.deg)),
    "Mrk 501": SkyCoord("16 53 52.22 +39 45 36.6", unit=(u.hourangle, u.deg)),
    "Boomerang": SkyCoord("22 28 44 +61 10 00", unit=(u.hourangle, u.deg)),
    "1ES 1959+650": SkyCoord("19 59 59.85 +65 08 54.7", unit=(u.hourangle, u.deg)),
    "LSI +61 303": SkyCoord("02 40 31.66 +61 13 45.6", unit=(u.hourangle, u.deg)),
    "Geminga": SkyCoord("06 33 54.15 +17 46 12.9", unit=(u.hourangle, u.deg)),
}

ALIASES = {
    "M1": "Crab",
    "M 1": "Crab", # Jan 2026 mount logs
}  # other names the mount's target_name may use, mapped to their SOURCES name

MAX_OFFSET = 1.0 # deg; largest wobble offset ever used plus a margin for pointing error. A pointing matches a source only if it's the only one this close

WOBBLE_SUFFIX = re.compile(r"\s+\d+(\.\d+)?\s*[NESW]$") # e.g. " 0.5 N", for target_names of sources not in the catalog


def parse_target_name(name):
    """
    The source in a mount target_name, without its wobble suffix, e.g. "M1 0.5 N" -> "Crab".

    The source is the longest SOURCES/ALIASES name that target_name starts with (followed by a space
    or the end), so any wobble suffix works for catalog sources. Other target_names lose a trailing
    "<offset> <N/E/S/W>" (WOBBLE_SUFFIX) and are otherwise kept as they are. The wobble itself comes
    from the pointing, see wobble_offset().
    """
    name = name.strip()
    known = sorted([*SOURCES, *ALIASES], key=len, reverse=True)
    prefix = next((n for n in known if name == n or (name.startswith(n) and name[len(n)].isspace())), None)
    if prefix is None:
        m = WOBBLE_SUFFIX.search(name)
        prefix = name[:m.start()] if m else name
    return ALIASES.get(prefix, prefix)


def match_source(pointing, max_offset=MAX_OFFSET):
    """
    The catalog source within max_offset of pointing. Raises if there's more than one, since the
    pointing alone can't tell them apart.

    Parameters:
        pointing: SkyCoord
        max_offset: max separation (deg)

    Returns:
        source name, or None if no source is within max_offset
    """
    names = list(SOURCES)
    catalog = SkyCoord([SOURCES[n] for n in names])
    separations = pointing.separation(catalog).deg
    candidates = [f"{names[i]} ({separations[i]:.2f} deg)" for i in np.flatnonzero(separations <= max_offset)]
    if len(candidates) > 1:
        raise ValueError(f"{len(candidates)} catalog sources within {max_offset} deg of {pointing.to_string('hmsdms')}: {', '.join(candidates)}")
    return names[int(np.argmin(separations))] if candidates else None


def wobble_offset(source, pointing):
    """
    pointing's offset from source.

    Returns:
        (separation (deg), position angle (deg east of north: 0 = N, 90 = E, 180 = S, 270 = W)),
        or (nan, nan) if source isn't in SOURCES or pointing is None
    """
    if source not in SOURCES or pointing is None:
        return np.nan, np.nan
    return SOURCES[source].separation(pointing).deg, SOURCES[source].position_angle(pointing).deg


def close_pairs(max_separation):
    """Every pair of catalog sources closer than max_separation (deg), as (name, name, separation (deg))."""
    names = list(SOURCES)
    catalog = SkyCoord([SOURCES[n] for n in names])
    i, j, separations, _ = search_around_sky(catalog, catalog, max_separation * u.deg)
    return [(names[a], names[b], s.deg) for a, b, s in zip(i, j, separations) if a < b]


# a pointing can only be within MAX_OFFSET of two sources if they're within 2*MAX_OFFSET of each other
for a, b, separation in close_pairs(2 * MAX_OFFSET):
    warnings.warn(
        f"heap.sources: {a} and {b} are {separation:.2f} deg apart (< 2*MAX_OFFSET), so runs near them "
        "with a blank target_name can't be identified from their pointing (see match_source())"
    )
