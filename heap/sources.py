"""sources

Catalog of observed sources, used to identify a run's target from its mount pointing when
target_name is blank (see heap.process_dataset.identify_source()). Names are the source names
analyses use (e.g. pff_analysis.ipynb's SOURCE); add an entry here for each newly observed source.
"""
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord

SOURCES = {
    "Crab": SkyCoord("05 34 31.78 +22 01 02.6", unit=(u.hourangle, u.deg)),
    "MGRO J2019+37": SkyCoord(304.83, 36.83, unit=u.deg),
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
}  # other names the mount's target_name may use, mapped to their SOURCES name

MAX_OFFSET = 1.0 # deg; farthest a pointing can be from a source and still match it (wobbles are ~0.5 deg)


def match_source(pointing, max_offset=MAX_OFFSET):
    """
    The catalog source nearest pointing.

    Parameters:
        pointing: SkyCoord
        max_offset: max separation (deg)

    Returns:
        source name, or None if no source is within max_offset
    """
    names = list(SOURCES)
    catalog = SkyCoord([SOURCES[n] for n in names])
    separations = pointing.separation(catalog).deg
    i = int(np.argmin(separations))
    return names[i] if separations[i] <= max_offset else None
