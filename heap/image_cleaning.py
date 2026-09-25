"""image_cleaning

Functions for image cleaning
"""
import numpy as np

_NEIGHBOR_DIRECTIONS = [
    (-1, -1), (-1, 0), (-1, 1),
    (0, -1),           (0, 1),
    (1, -1),  (1, 0),  (1, 1),
]

def _shifted(arr, dy, dx):
    # returns an array `out` of the same shape as arr such that out[i, j] == arr[i + dy, j + dx],
    # with out-of-bounds reads treated as 0/False.
    p = np.pad(arr, 1, constant_values=0)
    h, w = arr.shape
    return p[1 + dy: h + 1 + dy, 1 + dx: w + 1 + dx]

def _remove_pairs(mask, border_mask):
    # Remove 2-pixel image/border pairs.
    # Checked from both directions (from both pixels).
    # Returns a boolean mask.

    neighbor_arrays = [_shifted(mask, dy, dx) for dy, dx in _NEIGHBOR_DIRECTIONS]
    neighbor_count = sum(neighbor_arrays)
    border_neighbor_count = sum(_shifted(border_mask, dy, dx) for dy, dx in _NEIGHBOR_DIRECTIONS)
    image_mask = mask & (~border_mask)

    # a border pixel whose only neighbor (of either type) is this one lone pixel
    lone_border = border_mask & (neighbor_count == 1)
    # an image pixel whose only neighbor is a lone border pixel
    lone_image = image_mask & (neighbor_count == 1) & (border_neighbor_count == 1)
    pair_pixel = lone_border | lone_image

    # find the lone neighbor of each pair pixel, so it gets removed too
    pair_neighbor = np.zeros_like(mask)
    for (dy, dx), n_arr in zip(_NEIGHBOR_DIRECTIONS, neighbor_arrays):
        pair_dir = pair_pixel & (n_arr != 0)
        pair_neighbor |= _shifted(pair_dir, -dy, -dx)

    clean_mask = mask.copy()
    clean_mask[pair_pixel] = False
    clean_mask[pair_neighbor] = False

    return clean_mask

def _remove_isolated(mask):
    # Remove pixels with zero surviving neighbor/
    # Returns a boolean mask.

    neighbor_count = sum(_shifted(mask, dy, dx) for dy, dx in _NEIGHBOR_DIRECTIONS)
    isolated = mask & (neighbor_count == 0)

    clean_mask = mask.copy()
    clean_mask[isolated] = False

    return clean_mask

def threshold_clean(
        data: np.ndarray,
        pedvars: np.ndarray,
        image_threshold: float = 4.0,
        border_threshold: float = 2.0
):
    """
    Clean an image by defining a threshold for "image pixels" and some lower threshold for "border pixels"

    Parameters:
        data: data array
        pedvars: pedvar array
        image_threshold: threshold for labeling a pixel as an image pixel in units of pedvars
        border_threshold: threshold for labeling a pixel as a border pixel in units of pedvars

    Returns:
        cleaned: the cleaned image
        mask: the removed pixels
    """

    if border_threshold > image_threshold:
        raise ValueError("Border pixel threshold cannot be above image pixel threshold.")

    # read data
    arr = np.asarray(data)

    # scale thresholds to pedestal variance
    image_threshold = image_threshold * pedvars
    border_threshold = border_threshold * pedvars

    # create masks
    image_mask = arr >= image_threshold
    border_mask = (arr >= border_threshold) & (~image_mask)

    # border pixel must be next to image pixel
    p = np.pad(image_mask, 1, constant_values=0)
    neighbor_count = (
        p[:-2, :-2] + p[:-2, 1:-1] + p[:-2, 2:] +
        p[1:-1, :-2]               + p[1:-1, 2:] +
        p[2:  , :-2] + p[2:, 1:-1] + p[2:  , 2:]
    )
    border_mask = border_mask & (neighbor_count > 0)
    
    combined_mask = image_mask | border_mask
    combined_mask = _remove_pairs(combined_mask, border_mask)
    combined_mask = _remove_isolated(combined_mask)

    # apply masks
    cleaned = np.where(combined_mask, arr, np.nan)
    
    return cleaned, combined_mask