"""toy_events

Synthetic reconstructed events, shaped like heap.reconstruction.reconstruct_directions()'s
directions and heap.events.apply_cuts()'s array, for checking on/off counting (heap.significance)
with known pointings, source and background, e.g. wobble vs. on-source runs.

Events are drawn directly in each run's camera coordinates: a background spread over the camera with a
Gaussian radial acceptance, plus a point source smeared by a Gaussian PSF. Each event gets one
placeholder image with its centroid at the event's direction, so position cuts (distance, alpha,
miss) and the max distance cut are not meaningful for toy events.
"""
import astropy.units as u
import numpy as np
import pandas as pd

from heap.parameterize import CAMERA_HALF_WIDTH
from heap.significance import make_wcs

WOBBLE_ANGLES = {"N": 0, "E": 90, "S": 180, "W": 270} # position angle (deg, east of north); N/S = Dec +/- offset


def wobble_pointings(source, offset, directions=("N", "S")):
    """
    Pointings offset from source in each wobble direction.

    Parameters:
        source: SkyCoord
        offset: wobble offset (deg), e.g. 0.5 for wobble N/S at Dec +/-30'
        directions: keys of WOBBLE_ANGLES

    Returns:
        {direction: SkyCoord}
    """
    return {d: source.directional_offset_by(WOBBLE_ANGLES[d]*u.deg, offset*u.deg) for d in directions}


def simulate_run(run, pointing, source, n_background, n_signal, psf, acceptance_width, rng, date="toy", first_event=0):
    """
    One run's toy events, in its camera coordinates (deg from the pointing, see heap.significance.make_wcs()).

    Parameters:
        run: Run name
        pointing: SkyCoord of the camera center
        source: SkyCoord of the point source
        n_background: background events, inside the camera (|x|, |y| < CAMERA_HALF_WIDTH)
        n_signal: source events
        psf: Gaussian PSF sigma per axis (deg)
        acceptance_width: Gaussian sigma (deg) of the background's radial acceptance
        rng: numpy Generator
        date: Date column value
        first_event: Event number of the first event

    Returns:
        (directions, array): directions with Event, Date, Run, Xoffset, Yoffset, Signal (True for
        source events); array with one placeholder image per event (Event, Date, Run, Telescope,
        x_c, y_c, phi)
    """
    # background: uniform over the camera, thinned by the acceptance
    bkg = np.empty((0, 2))
    while len(bkg) < n_background:
        xy = rng.uniform(-CAMERA_HALF_WIDTH, CAMERA_HALF_WIDTH, size=(2*n_background, 2))
        keep = rng.uniform(size=len(xy)) < np.exp(-0.5*(xy**2).sum(axis=1)/acceptance_width**2)
        bkg = np.vstack([bkg, xy[keep]])
    bkg = bkg[:n_background]

    source_xy = np.array(make_wcs(pointing).wcs_world2pix(source.ra.deg, source.dec.deg, 1))
    sig = source_xy + rng.normal(scale=psf, size=(n_signal, 2))

    xy = np.vstack([bkg, sig])
    events = first_event + np.arange(len(xy))
    directions = pd.DataFrame({
        "Event": events, "Date": date, "Run": run, "Xoffset": xy[:, 0], "Yoffset": xy[:, 1],
        "Signal": np.r_[np.zeros(len(bkg), bool), np.ones(len(sig), bool)],
    })
    array = pd.DataFrame({
        "Event": events, "Date": date, "Run": run, "Telescope": "toy",
        "x_c": xy[:, 0], "y_c": xy[:, 1], "phi": rng.uniform(0, 180, len(xy)),
    })
    return directions, array


def simulate_runs(pointings, source, n_background, n_signal, psf=0.1, acceptance_width=2.0, seed=None, date="toy"):
    """
    simulate_run() for several runs, with Event numbers unique across them.

    Parameters:
        pointings: {Run: SkyCoord}
        see simulate_run() for the rest; seed seeds the numpy Generator

    Returns:
        (directions, array), each concatenated over runs
    """
    rng = np.random.default_rng(seed)
    directions, array = [], []
    for run, pointing in pointings.items():
        first = directions[-1].Event.max() + 1 if directions else 0
        d, a = simulate_run(run, pointing, source, n_background, n_signal, psf, acceptance_width, rng, date, first)
        directions.append(d)
        array.append(a)
    return pd.concat(directions, ignore_index=True), pd.concat(array, ignore_index=True)
