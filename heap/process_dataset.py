"""process_dataset

Runs process_image() over every event in a dataset, keeping every cleaned image
and collecting Hillas parameters into a table.
"""
import numpy as np

from heap.image_cleaning import threshold_clean
from heap.parameterize import calc_params


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
