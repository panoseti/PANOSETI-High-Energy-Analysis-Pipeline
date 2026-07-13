'''
Functions to calculate and visualise pedestals and pedestal variances.
'''

import numpy as np
import pandas as pd
from heap import pre_cleaning as pc
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.dates as mdates


def _gaussian(x, A, mu, sigma):
    return A * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))


def calculate_pedestal_and_pedvar(data, nsig=5.0, fit_gaussian=True):
    '''
    Returns pedestal and pedestal variance per pixel.
    Pedvar is obtained from a Gaussian fit of the distribution after removing outliers (abs > nsig * sigma).

    :param data: Data frames in format (n_frames, n_pixels)
    :param nsig: Number of sigmas above mean to define mask threshold (default = 5)
    :param fit_gaussian: Whether to perform Gaussian fit (default = True). If False, uses standard deviation of masked data.
    '''

    n_frames, n_pixels = data.shape

    mean_pixels_initial = np.nanmean(data, axis=0) # initial pedestal
    sigma_pixels_initial = np.abs(np.nanstd(data, axis=0)) # initial sigma for all pixels

    threshold = mean_pixels_initial + nsig * sigma_pixels_initial
    data_masked = np.where(data < threshold[None, :], data, np.nan)
    
    mean_pixels = np.nanmean(data_masked, axis=0)
    # If no Gaussian fit, compute mean and sigma on masked data directly (vectorized)
    if not fit_gaussian:
        sigma_pixels = np.abs(np.nanstd(data_masked, axis=0))
        return mean_pixels, sigma_pixels
    
    # If Gaussian fit, recompute mean on masked data for consistency
    sigma_pixels = np.zeros(n_pixels) # pedvar to be filled in

    # If Gaussian fit, create histogram and curve fitting per pixel
    for i in range(n_pixels):
        x = data_masked[:, i]
        x_clean = x[np.isfinite(x)]  # Remove nan values

        # initial estimates
        mu0 = mean_pixels[i]
        sigma0 = sigma_pixels_initial[i]

        if x_clean.size < 5:
            sigma_pixels[i] = sigma0
            continue

        hmin = -500
        hmax = 500
        hbins = 1000
        hist, edges = np.histogram(x_clean, bins=hbins, range=(hmin, hmax))
        centers = 0.5 * (edges[1:] + edges[:-1])

        try:
            p0 = [hist.max(), mu0, sigma0]
            popt, _ = curve_fit(_gaussian, centers, hist, p0=p0, maxfev=2000)
            sigma_pixels[i] = abs(popt[2])

        except RuntimeError:
            sigma_pixels[i] = sigma0

    return mean_pixels, sigma_pixels


def select_time_interval(filename, start_time, end_time):

    '''
    Returns data frames and timestamps from given interval.
    
    :param filename: Path to pff file
    :param start_time: Start of time interval in UTC in format: YYYY-MM-DDTHH:MM:SSZ or Timestamp
    :param end_time: End of time interval in UTC in format: YYYY-MM-DDTHH:MM:SSZ or Timestamp
    '''

    rawdata, metadata = pc.read_pff(filename)
    data, timestamps = pc.cut_pkt_loss(rawdata, metadata)

    timestamps_df = pd.to_datetime(timestamps, unit="s", utc=True)
    
    # Convert string times to Timestamp if needed
    if isinstance(start_time, str):
        start_time = pd.to_datetime(start_time, utc=True)
    if isinstance(end_time, str):
        end_time = pd.to_datetime(end_time, utc=True)
    
    start_idx = timestamps_df.searchsorted(start_time, side="left")
    end_idx = timestamps_df.searchsorted(end_time, side="right")

    timestamps_cut = timestamps[start_idx:end_idx]

    data_cut = data[start_idx:end_idx, :]
    print("shape of data:", np.shape(data_cut))

    return (data_cut, timestamps_cut)


def calculate_pedestal_and_pedvar_on_interval(filename, start_time, end_time, nsig=5, fit_gaussian=True):

    '''
    Returns pedestal (mean pixel value) and pedestal variance (standard deviation of pixel values) 
    for given time interval.
    
    :param filename: Path to pff file
    :param start_time: Start of time interval in UTC in format: YYYY-MM-DDTHH:MM:SSZ
    :param end_time: End of time interval in UTC in format: YYYY-MM-DDTHH:MM:SSZ
    :param nsig: Number of sigmas above mean to define mask threshold (default = 5)
    :param fit_gaussian: Whether to perform Gaussian fit (default = True). If False, uses standard deviation of masked data.
    '''

    data_cut, timestamps_cut = select_time_interval(filename, start_time, end_time)
    pedestal, pedvar = calculate_pedestal_and_pedvar(data_cut, nsig, fit_gaussian=fit_gaussian)

    return (pedestal, pedvar, data_cut, timestamps_cut)


def calculate_pedestal_and_pedvar_on_interval_from_array(
    data,
    timestamps,
    start_time,
    end_time,
    nsig=5,
    fit_gaussian=True,
):

    '''
    Returns pedestal (mean pixel value) and pedestal variance (standard deviation of pixel values)
    for a given time interval, using in-memory data and timestamps instead of reading a file.

    :param data: Data frames in format (n_frames, n_pixels)
    :param timestamps: 1D array of Unix timestamps in seconds corresponding to each frame
    :param start_time: Start of time interval in UTC in format: YYYY-MM-DDTHH:MM:SSZ or Timestamp
    :param end_time: End of time interval in UTC in format: YYYY-MM-DDTHH:MM:SSZ or Timestamp
    :param nsig: Number of sigmas above mean to define mask threshold (default = 5)
    :param fit_gaussian: Whether to perform Gaussian fit (default = True). If False, uses standard deviation of masked data.
    '''

    # Convert timestamps to pandas DateTimeIndex in UTC
    timestamps_df = pd.to_datetime(timestamps, unit="s", utc=True)

    # Convert string times to Timestamp if needed
    if isinstance(start_time, str):
        start_time = pd.to_datetime(start_time, utc=True)
    if isinstance(end_time, str):
        end_time = pd.to_datetime(end_time, utc=True)

    # Find index range for the requested interval
    start_idx = timestamps_df.searchsorted(start_time, side="left")
    end_idx = timestamps_df.searchsorted(end_time, side="right")

    data_cut = data[start_idx:end_idx, :]
    timestamps_cut = timestamps[start_idx:end_idx]

    pedestal, pedvar = calculate_pedestal_and_pedvar(data_cut, nsig, fit_gaussian=fit_gaussian)

    return (pedestal, pedvar, data_cut, timestamps_cut)


def calculate_pedvar_diff(data_pre_flip, data_post_flip, nsig=5, fit_gaussian=True, diff_sign=1.0, squared=True, rotate_post_flip=False):

    '''
    Returns difference between pedestal variance before and after meridian flip.
    
    :param data_pre_flip: Data frames in format (n_frames, n_pixels) before flip
    :param data_post_flip:Data frames in format (n_frames, n_pixels) after flip
    :param nsig:  Number of sigmas above mean to define mask threshold (default = 5)
    :param fit_gaussian: Whether to perform Gaussian fit (default = True). If False, uses standard deviation of masked data.
    :param diff_sign: Sign of difference (default=1.0)
    :param squared: Whether the difference should be squared (default=True)
    '''

    _, pedvar_pre_flip = calculate_pedestal_and_pedvar(data_pre_flip, nsig, fit_gaussian=fit_gaussian)
    _, pedvar_post_flip = calculate_pedestal_and_pedvar(data_post_flip, nsig, fit_gaussian=fit_gaussian)

    if rotate_post_flip:
        pedvar_post_flip = np.rot90(pedvar_post_flip, k=2, axes=(0, 1))

    pedvar_diff = diff_sign * (pedvar_pre_flip - pedvar_post_flip)
    
    if squared:
        #pedvar_diff = pedvar_diff * np.abs(pedvar_diff)
        pedvar_diff = np.sign(pedvar_diff) * pedvar_diff**2

    return pedvar_diff


def calculate_pedvar_diff_on_interval(filename, start_time_pre, end_time_pre, start_time_post, end_time_post, nsig=5, fit_gaussian=True, diff_sign=1.0, squared=True, rotate_post_flip=False):

    '''
    Returns difference between pedestal variance before and after meridian flip.
    
    :param filename: Path to pff file
    :param start_time_pre: Start of time interval pre-flip in UTC in format: YYYY-MM-DDTHH:MM:SSZ
    :param end_time_pre: End of time interval pre-flip in UTC in format: YYYY-MM-DDTHH:MM:SSZ
    :param start_time_post: Start of time interval post-flip in UTC in format: YYYY-MM-DDTHH:MM:SSZ
    :param end_time_post: End of time interval post-flip in UTC in format: YYYY-MM-DDTHH:MM:SSZ
    :param nsig:  Number of sigmas above mean to define mask threshold (default = 5)
    :param fit_gaussian: Whether to perform Gaussian fit (default = True). If False, uses standard deviation of masked data.
    :param diff_sign: Sign of difference (default=1.0)
    :param squared: Whether the difference should be squared (default=True)
    '''

    _, pedvar_pre_flip, _, _ = calculate_pedestal_and_pedvar_on_interval(filename, start_time_pre, end_time_pre, nsig, fit_gaussian=fit_gaussian)
    _, pedvar_post_flip, _, _ = calculate_pedestal_and_pedvar_on_interval(filename, start_time_post, end_time_post, nsig, fit_gaussian=fit_gaussian)

    if rotate_post_flip:
        # Reshape to 2D, rotate, then flatten back
        pedvar_post_flip_2d = pedvar_post_flip.reshape((32, 32))
        pedvar_post_flip_2d = np.rot90(pedvar_post_flip_2d, k=2)
        pedvar_post_flip = pedvar_post_flip_2d.flatten()

    pedvar_diff = diff_sign * (pedvar_pre_flip - pedvar_post_flip)
    
    if squared:
        #pedvar_diff = pedvar_diff * np.abs(pedvar_diff)
        pedvar_diff = np.sign(pedvar_diff) * pedvar_diff**2

    return pedvar_diff


def calculate_pedvar_diff_on_interval_from_array(
    data,
    timestamps,
    start_time_pre,
    end_time_pre,
    start_time_post,
    end_time_post,
    nsig=5,
    fit_gaussian=True,
    diff_sign=1.0,
    squared=True,
    rotate_post_flip=False,
):

    '''
    Returns difference between pedestal variance before and after meridian flip,
    using in-memory data and timestamps instead of reading from file.

    :param data: Data frames in format (n_frames, n_pixels)
    :param timestamps: 1D array of Unix timestamps in seconds corresponding to each frame
    :param start_time_pre: Start of time interval pre-flip in UTC (string or Timestamp)
    :param end_time_pre: End of time interval pre-flip in UTC (string or Timestamp)
    :param start_time_post: Start of time interval post-flip in UTC (string or Timestamp)
    :param end_time_post: End of time interval post-flip in UTC (string or Timestamp)
    :param nsig:  Number of sigmas above mean to define mask threshold (default = 5)
    :param fit_gaussian: Whether to perform Gaussian fit (default = True). If False, uses standard deviation of masked data.
    :param diff_sign: Sign of difference (default=1.0)
    :param squared: Whether the difference should be squared (default=True)
    :param rotate_post_flip: If True, rotate post-flip pedvar map by 180° (assumes 32x32 camera)
    '''

    # Pre-flip interval
    _, pedvar_pre_flip, _, _ = calculate_pedestal_and_pedvar_on_interval_from_array(
        data,
        timestamps,
        start_time_pre,
        end_time_pre,
        nsig=nsig,
        fit_gaussian=fit_gaussian,
    )

    # Post-flip interval
    _, pedvar_post_flip, _, _ = calculate_pedestal_and_pedvar_on_interval_from_array(
        data,
        timestamps,
        start_time_post,
        end_time_post,
        nsig=nsig,
        fit_gaussian=fit_gaussian,
    )

    if rotate_post_flip:
        # Reshape to 2D, rotate, then flatten back (assumes 32x32 layout)
        pedvar_post_flip_2d = pedvar_post_flip.reshape((32, 32))
        pedvar_post_flip_2d = np.rot90(pedvar_post_flip_2d, k=2)
        pedvar_post_flip = pedvar_post_flip_2d.flatten()

    pedvar_diff = diff_sign * (pedvar_pre_flip - pedvar_post_flip)

    if squared:
        pedvar_diff = np.sign(pedvar_diff) * pedvar_diff**2

    return pedvar_diff


def plot_trigger_rate(filename, hk_filename, bin_width=1.0, start_time=None, end_time=None):
    
    '''
    Plots trigger rate over time.
    
    :param filename: Path to pff file
    :param bin_width: Width of time bins in seconds (default = 1.0)
    :param start_time: Start of time interval in UTC in format: YYYY-MM-DDTHH:MM:SSZ (default = None)
    :param end_time: End of time interval in UTC in format: YYYY-MM-DDTHH:MM:SSZ (default = None)
    '''

    rawdata, metadata = pc.read_pff(filename)
    data, timestamps = pc.cut_pkt_loss(rawdata, metadata)
    hk = pc.read_pff_hk(hk_filename)
    _, _, merflip_start, merflip_end = pc.cut_meridian_flip(data, timestamps, hk)
    merflip_start = pd.to_datetime(merflip_start, unit="s", utc=True)
    merflip_end = pd.to_datetime(merflip_end, unit="s", utc=True)
    print(f"Meridian flip start (UTC): {merflip_start}")
    print(f"Meridian flip end (UTC): {merflip_end}")

    timestamps_df = pd.to_datetime(timestamps, unit="s", utc=True)

    if start_time is not None:
        start_idx = timestamps_df.searchsorted(start_time, side="left")
        timestamps_df = timestamps_df[start_idx:]

    if end_time is not None:
        end_idx = timestamps_df.searchsorted(end_time, side="right")
        timestamps_df = timestamps_df[:end_idx]

    # Compute trigger rate
    time_bins = pd.date_range(start=timestamps_df.min(), end=timestamps_df.max(), freq=pd.Timedelta(seconds=bin_width))
    trigger_counts = np.histogram(timestamps_df, bins=time_bins)[0]
    trigger_rate = trigger_counts / bin_width

    # Plotting
    plt.figure(figsize=(10, 5))
    plt.plot(time_bins[:-1], trigger_rate, drawstyle='steps-post')
    
    # Add transparent box for meridian flip period
    ax = plt.gca()
    ax.axvspan(merflip_start, merflip_end, alpha=0.2, color='red', label='Meridian Flip')
    
    plt.xlabel('Time (UTC)')
    plt.ylabel('Trigger Rate (Hz)')
    plt.title('Trigger Rate Over Time')
    
    # Format x-axis to display time more clearly
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
    ax.xaxis.set_major_locator(mdates.AutoDateLocator())
    plt.xticks(rotation=45, ha='right')
    
    plt.legend()
    plt.grid()
    plt.tight_layout()

    return None







    