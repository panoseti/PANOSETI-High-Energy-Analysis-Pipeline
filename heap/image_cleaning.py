"""image_cleaning

Functions for image cleaning
"""
import numpy as np

def _remove_islands(mask):
    # Remove pixels whose neighbors register zero digital counts. Returns a boolean mask.

    # pad mask by 1
    p = np.pad(mask, 1, constant_values=0)
    neighbor_count = (
        p[:-2, :-2] + p[:-2, 1:-1] + p[:-2, 2:] +
        p[1:-1, :-2]               + p[1:-1, 2:] +
        p[2:  , :-2] + p[2:, 1:-1] + p[2:  , 2:]
    )

    # isolated = pixel is True and all 8 neighbors are False
    isolated = (p[1:-1, 1:-1] != 0) & (neighbor_count == 0)
    clean_mask = mask.copy()
    clean_mask[isolated] = False
    
    return clean_mask

def threshold_clean(
        data: np.ndarray,
        pedvars: np.ndarray,
        image_threshold: float = 5.0,
        border_threshold: float = 2.5
):
    """
    Clean an image by defining a threshold for "image pixels" and some lower threshold for "border pixels"

    Parameters:
        data:
        image_threshold: 
        border_threshold:
        
    Returns:
        cleaned: the cleaned image
        mask: the removed pixels
    """

    if not border_threshold < image_threshold:
        raise ValueError("Border pixel threshold must be below image pixel threshold.")

    # read data
    arr = np.asarray(data)

    # scale thresholds to pedestal variance
    image_threshold = image_threshold * pedvars
    border_threshold = border_threshold * pedvars

    # create masks
    image_mask = arr >= image_threshold
    border_mask = (arr >= border_threshold) & (~image_mask)
    combined_mask = image_mask | border_mask
    combined_mask = _remove_islands(combined_mask)

    # apply masks
    cleaned = np.where(combined_mask, arr, np.nan)
    
    return cleaned, combined_mask