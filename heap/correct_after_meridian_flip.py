'''
Functions to rotate frames after the meridian flip by 180°.
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def divide_data_meridian_flip(data, timestamps, hk, time_after=60, plotting=False, path=None):

    '''
    Exactly the same as cut_meridian_flip from pre_cleaning.py but returns the data after cut in two chunks - before and after meridian flip cuts.
    data: data array in (n_timesteps, n_pix*n_pix) 
    timestamps: timestamp array in (n_timestamps)
    hk: housekeeping pff data as read in
    time_after: time frame in s to cut out after tracking is enabled again
    plotting: plots RA and dec (which should stay constant during tracking) zoomed in around meridian flip timeframe if true
    path: path to save the plot to
    returns: tuple with data tuple (2 elements) with data (array) before and after flip cut, timestamp tuple (2 elements) with timestamps (array) before and after cut and UTC timestamps of start and end of the flip
    '''

    # Get timestamps, differences between timestamps and tracking booleans
    time = np.asarray(hk["MOUNT_GATTINI"]["Computer_UTC"])
    time_diff = np.diff(time)
    track = np.asarray(hk["MOUNT_GATTINI"]["tracking"])

    # Filter data where tracking was off and add a point before and after flip
    no_track = np.where(track == "False")[0]
    no_track = np.insert(no_track, [0, len(no_track)], [min(no_track)-1, max(no_track)+1])

    # Also cut time when guiding was recalibrating (time_after)
    instert_val_num = int(time_after/np.median(time_diff))+1 
    no_track = np.append(no_track, np.arange(max(no_track), max(no_track)+1+instert_val_num, 1))

    # Print info on meridian flip
    time_pd = pd.to_datetime(time, unit='s', utc=True).tz_convert('America/Los_Angeles')
    print("Start of meridian flip (local time): ", time_pd[no_track[0]], ", End: ",time_pd[no_track[-1]])
    print("Duration of meridian flip cut: ", time[no_track[-1]]-time[no_track[0]])

    # NEW: Make a before/after dataset
    data_before_flip = data[(timestamps<time[no_track[0]])]
    data_after_flip = data[(timestamps>time[no_track[-1]])]
    data_post_cuts = (data_before_flip, data_after_flip)
    timestamps_before_flip = timestamps[(timestamps<time[no_track[0]])]
    timestamps_after_flip = timestamps[(timestamps>time[no_track[-1]])]
    timestamps_post_cuts = (timestamps_before_flip, timestamps_after_flip)
    print("Number of cut out events from meridian flip: ", len(data)-(len(data_before_flip)+len(data_after_flip)), " (", round((len(data)-(len(data_before_flip)+len(data_after_flip)))/len(data)*100, 2),"%)")

    # Plotting to visualise cut
    if plotting:
        ra=hk["MOUNT_GATTINI"]["ra_hours"]
        dec=hk["MOUNT_GATTINI"]["dec_deg"]
        fig, ax1 = plt.subplots()
        color = 'tab:blue'
        ax1.set_xlabel('time')
        ax1.set_ylabel('RA [h]', color=color)
        ax1.scatter(time_pd,ra,marker=".",color=color)
        ax1.tick_params(axis='y', labelcolor=color)
        
        ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis
        
        color = 'tab:red'
        ax2.set_ylabel('DEC [°]', color=color)  # we already handled the x-label with ax1
        ax2.scatter(time_pd,dec,marker=".",color=color)
        ax2.tick_params(axis='y', labelcolor=color)
        
        ax1.fill_between(time_pd[no_track],y1=13,alpha=0.3,label="Meridian flip cut",color="grey")
        ax1.set_xlim(time_pd[no_track[0]-20],time_pd[no_track[-1]+20])
        fig.legend(loc=1)
        fig.autofmt_xdate()
        fig.tight_layout()
        plt.savefig(path+"_meridian_flip.png",dpi=300)
        plt.show()

    return (data_post_cuts, timestamps_post_cuts, time[no_track[0]], time[no_track[-1]])


def rotate_after_flip(data_post_cuts, timestamps_post_cuts, n_pix=32):

    '''
    Rotates data after meridian flip by 180° and appends to data before meridian flip.
    data_post_cuts: tuple (2 elements) of data (array of dims (n_timestamps, n_pix*n_pix)) before and after flip cuts
    timestamps_post_cuts: tuple (2 elements) of timestamps (array of dims (n_timestamps)) before and after flip cuts
    returns: tupe of concatenated data and timestamps from before the flip and after the flip rotated by 180°
    '''
    
    # Unpack data before and after cut
    data_before_flip, data_after_flip = data_post_cuts
    timestamps_before_flip, timestamps_after_flip = timestamps_post_cuts

    # Reshape from (m, n) -> (m, n_pix, n_pix)
    reshaped_data_after_flip = np.reshape(data_after_flip, (len(data_after_flip), n_pix, n_pix))

    # Rotate (n_pix, n_pix) by 180°
    rotated_data_after_flip = np.rot90(reshaped_data_after_flip, k=2, axes=(1, 2))

    # Reshape from (m, n_pix, n_pix) -> (m, n)
    reshaped_rotated_data_after_flip = np.reshape(rotated_data_after_flip, (len(data_after_flip), n_pix*n_pix))

    # Append along axis 0 to xxx_before_flip
    mer_flip_corr_data = np.concatenate((data_before_flip, reshaped_rotated_data_after_flip), axis=0)
    mer_flip_corr_timestamps = np.concatenate((timestamps_before_flip, timestamps_after_flip), axis=0)

    return (mer_flip_corr_data, mer_flip_corr_timestamps)



