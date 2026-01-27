import pre_cleaning
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import os

def load_telescope_tv(module_str, base_dir, plotting=True):
    """
    Reads in all files for one telescope, applies
    read_pff -> cut_pkt_loss_old (tv-timestamps) -> spike_cut
    Returns combined timestamps and data
    """
    files = sorted(base_dir.glob(f"*{module_str}*.pff"))
    print(f"{module_str}: {len(files)} files found")

    all_timestamps = []
    all_data=[]
    os.system("mkdir -p "+ f"{base_dir}/{module_str}/")

    for f in files:
        print(f"  -> {f}")
        data, metadata = pre_cleaning.read_pff(str(f))
        data, timestamps = pre_cleaning.cut_pkt_loss_old(data, metadata)
        
        plot_path = f"{base_dir}/{module_str}/"
        data_clean, timestamps_clean = pre_cleaning.spike_cut(
            data,
            timestamps,
            plotting=plotting,
            path=plot_path
        )

        all_timestamps.append(timestamps_clean)
        all_data.append(data_clean)

    if not all_timestamps:
        raise RuntimeError(f"No files for {module_str} found.")

    # combine to one array
    data_all=np.concatenate(all_data)
    timestamps_all = np.concatenate(all_timestamps)
    return data_all, timestamps_all

def time_offset(timestamps1, timestamps2, window=0.02, plotting=False):
    """
    Determines time dependent timing offset between two telescopes t1, t2 using coincident events within a large window
    Returns timing difference dt=t1-t2 and corresponding timestamps1 and timestamps2 in unix time
    Fast version with using double loop
    """
    t1 = np.asarray(timestamps1)
    t2 = np.asarray(timestamps2)

    # finde für jedes t1 den Suchbereich in t2
    left  = np.searchsorted(t2, t1 - window, side="left")
    right = np.searchsorted(t2, t1 + window, side="right")

    idx1 = []
    idx2 = []

    for i, (l, r) in enumerate(zip(left, right)):
        if l < r:
            idx1.extend([i] * (r - l))
            idx2.extend(range(l, r))

    idx1 = np.asarray(idx1)
    idx2 = np.asarray(idx2)

    time_coinc1 = t1[idx1]
    time_coinc2 = t2[idx2]
    dt = time_coinc1 - time_coinc2

    # Pandas-Zeitachsen nur einmal erzeugen
    time_coinc_pd1 = pd.to_datetime(time_coinc1, unit='s', utc=True)\
                         .tz_convert('America/Los_Angeles')

    if plotting:
        plt.scatter(time_coinc_pd1, dt, marker=".", s=5)
        plt.xlabel("time")
        plt.ylabel("time difference [s]")
        plt.xticks(rotation=30)
        plt.tight_layout()
        plt.show()

    return time_coinc1, time_coinc2, dt

def correct_time(timestamps1, timestamps2,  plot_name, base_dir, window=0.02, bin_width=120, plotting=True):
    """
    Determines time dependant timing offset between two telescopes my matching events within a large window
    Corrects timestamp1 with the determined offset function and returns corrected timestamp 1
    """
    #get timing differences and coincident timestamps of t1
    time_coinc1,_,dt=time_offset(timestamps1,timestamps2,window=window)
    #calculating median offset over certain time range
    y,x=np.histogram(time_coinc1,bins=np.arange(min(time_coinc1),max(time_coinc1)+bin_width+1,bin_width))
    dt_median=np.zeros(len(x)-1)
    for i in range(len(x)-1):
        dt_median[i]=np.median(dt[(time_coinc1>=x[i]) & (time_coinc1<x[i+1])])
    rms=np.sqrt(np.mean(np.square(dt_median[~np.isnan(dt_median)])))
    sigma=np.std(dt_median[~np.isnan(dt_median)])
    print("Offset: RMS= ",round(rms,5),"s, std=",round(sigma, 5),"s")
    #Correcting timestamps1 with determined offset function
    t_corr=np.digitize(timestamps1,x)-1
    timestamps1_corr=timestamps1-dt_median[t_corr]
    #calculate corrected timing differences
    time_coinc1_corr,_,dt_corr=time_offset(timestamps1_corr,timestamps2,window=window)
    #check that time offset is 0 now
    y_corr,x_corr=np.histogram(time_coinc1_corr,bins=np.arange(min(time_coinc1_corr),max(time_coinc1_corr)+bin_width+1,bin_width))
    dt_median_corr=np.zeros(len(x_corr)-1)
    for i in range(len(x_corr)-1):
        dt_median_corr[i]=np.median(dt_corr[(time_coinc1_corr>=x_corr[i]) & (time_coinc1_corr<x_corr[i+1])])
    rms_corr=np.sqrt(np.mean(np.square(dt_median_corr[~np.isnan(dt_median_corr)])))
    sigma_corr=np.std(dt_median_corr[~np.isnan(dt_median_corr)])
    print("Corrected Offset: RMS= ",round(rms_corr,5),"s, std=",round(sigma_corr, 5),"s")
    if plotting:
        time_coinc_pd1 = pd.to_datetime(time_coinc1, unit='s', utc=True).tz_convert('America/Los_Angeles')
        time_coinc_pd1_corr = pd.to_datetime(time_coinc1_corr, unit='s', utc=True).tz_convert('America/Los_Angeles')
        x_pd=pd.to_datetime(x[:-1], unit='s', utc=True).tz_convert('America/Los_Angeles')
        fig,ax=plt.subplots(2,1,sharex=True)
        ax[0].scatter(time_coinc_pd1,dt,marker=".")
        ax[0].step(x_pd,dt_median,label="dt median, RMS="+str(round(rms,5))+"s",color="red")
        ax[0].set_ylabel("time offset [s]")
        ax[0].grid()
        ax[0].set_title("Before Correction")
        ax[1].scatter(time_coinc_pd1_corr,dt_corr,marker=".")
        ax[1].set_xlabel("Time")
        ax[1].set_ylabel("time offset [s]")
        ax[1].set_title("After Correction")
        ax[1].grid()
        ax[0].legend()
        fig.suptitle(str(base_dir)+", "+plot_name)
        fig.tight_layout()
        path=f"{base_dir}/"
        plt.savefig(path+plot_name+"_dt_corr.png",dpi=300)
        plt.show()
    return(timestamps1_corr)

def match_coinc(timestamps1, data1, timestamps2, data2, window=0.001, tz='America/Los_Angeles'):
    """
    Assuming corrected timestamps, looks for coincident events between 2 telescopes within smaller set window
    Returns timestamps (local time) and data of both telescopes and timing difference for coincident events
    """
    print("Number of events for telescope 1:", len(timestamps1))
    print("Number of events for telescope 2:", len(timestamps2))

    # Candidate ranges in t2 for each t1
    left  = np.searchsorted(timestamps2, timestamps1 - window, side="left")
    right = np.searchsorted(timestamps2, timestamps1 + window, side="right")

    # Build matching index pairs (i in ttimestamps1, j in timestamps2)
    idx1 = []
    idx2 = []
    for i, (l, r) in enumerate(zip(left, right)):
        if l < r:
            idx1.extend([i] * (r - l))
            idx2.extend(range(l, r))

    idx1 = np.asarray(idx1, dtype=np.int64)
    idx2 = np.asarray(idx2, dtype=np.int64)

    # Gather results
    time_coinc1 = timestamps1[idx1]
    time_coinc2 = timestamps2[idx2]
    data_coinc1 = data1[idx1]
    data_coinc2 = data2[idx2]

    # Convert to local time (only for matched timestamps)
    time_coinc_pd1 = pd.to_datetime(time_coinc1, unit='s', utc=True).tz_convert(tz)
    time_coinc_pd2 = pd.to_datetime(time_coinc2, unit='s', utc=True).tz_convert(tz)

    ncoinc = len(time_coinc1)
    p1 = (ncoinc * 100 / len(timestamps1)) if len(timestamps1) else 0.0
    p2 = (ncoinc * 100 / len(timestamps2)) if len(timestamps2) else 0.0
    print(f"Number of coincident events within {window}s: {ncoinc} "
          f"( ~{p1:.1f}% of all telescope 1 events and ~{p2:.1f}% of all telescope 2 events.)")

    return time_coinc_pd1, data_coinc1, time_coinc_pd2, data_coinc2