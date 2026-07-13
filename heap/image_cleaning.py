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

def _remove_islands(mask, border_mask):
    # Remove pixels whose neighbors register zero digital counts (isolated single pixels),
    # and remove 2-pixel islands (a pixel with exactly one true neighbor) whenever the
    # pixel itself is a border pixel - both pixels of the island are dropped.
    # Returns a boolean mask.

    neighbor_arrays = [_shifted(mask, dy, dx) for dy, dx in _NEIGHBOR_DIRECTIONS]
    neighbor_count = sum(neighbor_arrays)

    # isolated = pixel is True and all 8 neighbors are False
    isolated = mask & (neighbor_count == 0)

    # 2-pixel islands: pixel has exactly one true neighbor, and this pixel is a border pixel
    single_neighbor = mask & (neighbor_count == 1)
    qualifying = single_neighbor & border_mask

    # find the lone neighbor of each qualifying pixel, so it can be removed too
    neighbor_to_remove = np.zeros_like(mask)
    for (dy, dx), n_arr in zip(_NEIGHBOR_DIRECTIONS, neighbor_arrays):
        qualifying_dir = qualifying & (n_arr != 0)
        neighbor_to_remove |= _shifted(qualifying_dir, -dy, -dx)

    clean_mask = mask.copy()
    clean_mask[isolated] = False
    clean_mask[qualifying] = False
    clean_mask[neighbor_to_remove] = False

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
        image_threshold: threshold for labeling a pixel as an image pixel in units of pedvar
        border_threshold: threshold for labeling a pixel as a border pixel in units of pedvar
        
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
    combined_mask = _remove_islands(combined_mask, border_mask)

    # apply masks
    cleaned = np.where(combined_mask, arr, np.nan)
    
    return cleaned, combined_mask