"""parameterize

Functions for creating image pedestals and calculating pedestal variances
"""
import numpy as np
from matplotlib import pyplot as plt

def make_pedestal_images(data: np.ndarray):
    """
    Make a map of pedestals for every pixel over many events

    Parameters:
        data: array of images to analyze pixel data. numpy array with shape (N, 1024)

    Returns
        ped: map of pedestals for every pixel. numpy array with shape (N, 1024)
        ped_var: 
    """
    
    # check dimension of data
    data = np.asarray(data, dtype=float)
    if data.ndim != 2:
        raise ValueError(f"expecting data with shape (N, P), got {data.shape}")

    # initial pedestal stats
    ped_means = np.nanmean(data, axis=0)
    ped_stds = np.nanstd(data, axis=0)

    # mask out Cherenkov events
    threshold = ped_means + 5.0 * ped_stds
    masked = np.where(data < threshold[None, :], data, np.nan)

    # pedestal variance after filtering
    ped = np.nanmean(masked, axis=0)
    ped_var = np.nanstd(masked, axis=0)

    # return arrays (not lists) for downstream numeric work
    return ped, ped_var

def plot_distribution(
        data: np.ndarray,
        pix: int=0,
        bins: int = 100
    ):
    """
    Plot the pedestal distribution of a single pixel over many events

    Parameters:
        data: array of images to analyze pixel data. numpy array with shape (N, 1024)
        pix: pixel number in image
        bins: number of histogram bins

    Returns pyplot figure
    """

    # check dimension of data
    data = np.asarray(data, dtype=float)
    if data.ndim != 2:
        raise ValueError(f"expecting data with shape (N, P), got {data.shape}")
    N, P = data.shape
    pix = int(pix)
    if not (0 <= pix < P):
        raise IndexError(f"pix index out of range 0..{P-1}: {pix}")

    # initial pedestal stats
    ped_mean = np.nanmean(data[:,pix])
    ped_std = np.nanstd(data[:,pix])

    # mask out Cherenkov events
    threshold = ped_mean + 5.0 * ped_std
    masked = np.where(data[:,pix] < threshold, data[:,pix], np.nan)

    # pedestal variance after filtering
    ped = np.nanmean(masked)
    ped_var = np.nanstd(masked)

    # plotting
    fig, axs = plt.subplots(1,2, figsize=(10,5))
    ax1, ax2 = axs
    low=ped_mean - ped_std
    high=ped_mean + 3*ped_std
    bin_edges = np.linspace(low, 2*high, bins)  
    ax1.hist(data[:,pix], bins=bin_edges)
    ax1.axvline(threshold, color='red')
    n_total = int(np.count_nonzero(np.isfinite(data[:, pix])))
    ax1.set_title("Pedestal: pixel {} (N={})".format(pix, n_total))
    ax2.hist(masked, bins=bin_edges)
    n_masked = int(np.count_nonzero(np.isfinite(masked)))
    ax2.set_title("Pedestal, Cherenkov removed: pixel {} (N={})".format(pix, n_masked))

    for ax in axs:
        ax.set_yscale('log')
        ax.set_xlabel('ADC')
        ax.set_ylabel('Counts')

    plt.tight_layout()
    plt.close()
    
    # printing
    print("Mean: ", ped_mean)
    print("Std:  ", ped_std)

    print("Updated Mean: ", ped)
    print("Updated Std:  ", ped_var)

    return fig