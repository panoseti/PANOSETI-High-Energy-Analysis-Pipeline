'''
These functions are supposed to to be the first processing scripts we apply on the data.
They include reading in the pff files and filtering out frames with package loss, spikes in the rates and frames during the meridian flip.
For these processing steps only metadata and housekeeping data is used, no need to look at the events for it.
'''


import numpy as np
import pypff
import pandas as pd
import matplotlib.pyplot as plt


def read_pff(filename):
    #Reads in pff file and returns data and metadata
    pff = pypff.io.datapff(filename)
    data, metadata = pff.readpff(metadata=True)
    print("Data read in, number of events: ",len(data))
    return(data,metadata)

def cut_pkt_loss(data,metadata):
    #cut out frames with package loss and return data and timestamps
    #get pkt number, if it is 0 the corresponding quabo will be missing in the data
    pkt_num_0=metadata["quabo_0"]["pkt_num"]
    pkt_num_1=metadata["quabo_1"]["pkt_num"]
    pkt_num_2=metadata["quabo_2"]["pkt_num"]
    pkt_num_3=metadata["quabo_3"]["pkt_num"]
    #get timestamps for every quabo
    timestamps0=metadata["quabo_0"]["tv_sec"]+metadata["quabo_0"]["tv_usec"]*10**(-6)
    timestamps1=metadata["quabo_1"]["tv_sec"]+metadata["quabo_1"]["tv_usec"]*10**(-6)
    timestamps2=metadata["quabo_2"]["tv_sec"]+metadata["quabo_2"]["tv_usec"]*10**(-6)
    timestamps3=metadata["quabo_3"]["tv_sec"]+metadata["quabo_3"]["tv_usec"]*10**(-6)
    #it is different for every quabo, I guess it is the time the paket arrived on the daq node, not the actual trigger time, taking min for now, but needs to be reviewed
    timestamps=np.min([timestamps0,timestamps1,timestamps2,timestamps3],axis=0)
    #cut out all data and timestamps with package loss for at least one quabo
    data_cut=data[(pkt_num_0!=0)&(pkt_num_1!=0)&(pkt_num_2!=0)&(pkt_num_3!=0)]
    timestamps_cut=timestamps[(pkt_num_0!=0)&(pkt_num_1!=0)&(pkt_num_2!=0)&(pkt_num_3!=0)]
    print("Number of Events with package loss: ",len(data)-len(data_cut)," (",round((len(data)-len(data_cut))/len(data)*100,2),"%)")
    return(data_cut,timestamps_cut)

def spike_cut(data,timestamps,bin_width=30,rate_cut=2,plotting=False,path=None):
    '''
    Removes spikes in trigger rate caused e.g. by planes passing through the filed of view
    bin_width: time width in s to determine the trigger rate, a plane takes ~ 10s so should be higher than that
    rate_cut: threshold at which to cut off spikes, excepted trigger rate without mood ~ 0.5-1 Hz
    plotting: plots trigger rate over time before and after cut if True
    path: path to save the plot to
    '''
    bins = np.arange(timestamps.min(), timestamps.max() + bin_width, bin_width)
    counts, _ = np.histogram(timestamps, bins=bins)
    rate = counts / bin_width  # Hz
    bad_bins = rate > rate_cut
    bin_indices = np.digitize(timestamps, bins) - 1
    #masking bins which exceed the set cut threshold
    mask = ~bad_bins[bin_indices]
    timestamps_filtered = timestamps[mask]
    data_filtered=data[mask]
    bin_centers = bins[:-1] + bin_width/2
    time = pd.to_datetime(bin_centers, unit='s', utc=True).tz_convert('America/Los_Angeles')
    if plotting:
        counts, _ = np.histogram(timestamps_filtered, bins=bins)
        rate_filtered = counts / bin_width
        fig,ax=plt.subplots(2,1,sharex=True)
        ax[0].step(time,rate)
        ax[0].set_ylabel("Trigger Rate [Hz]")
        ax[0].set_title("Before Cut")
        ax[1].step(time,rate_filtered)
        ax[1].set_xlabel("Time")
        ax[1].set_ylabel("Trigger Rate [Hz]")
        ax[1].set_title("After Cut")
        fig.tight_layout()
        plt.savefig(path+"_spike_cut.png",dpi=300)
        plt.show()
    timestamps_pd=pd.to_datetime(timestamps, unit='s', utc=True).tz_convert('America/Los_Angeles')
    print("Rate spikes at: ",time[rate>2])
    print("Number of spike cut out events: ",len(data)-len(data_filtered)," (",round((len(data)-len(data_filtered))/len(data)*100,2),"%)")
    return(data_filtered,timestamps_filtered)

def read_pff_hk(filename="hk.pff"):
    #Read in housekeeping pff file and return data, default called hk.pff, gonna need to adjust path
    hkpff = pypff.io.hkpff(filename)
    hk = hkpff.readhk()
    return(hk)

def cut_meridian_flip(data,timestamps,hk, time_after=60, plotting=False,path=None):
    '''
    Uses hk data from observation to locate meridian flip and cut out data during this time
    During this time the mount doesn't track the object so using this parameter as an indicator for it
    After meridian flip the guiding need to be calibrated again which takes some time ~ 1 min
    hk data taken independently of telescope triggers, one value every ~9 s
    hk: housekeeping pff data as read in
    time_after: time frame in s to cut out after tracking is enabled again
    plotting: plots RA and dec (which should stay constant during tracking) zoomed in around meridian flip timeframe if true
    path: path to save the plot to
    returns: data and timestamps after cut out meridian flip and the UTC timestamps of start and end of the flip
    '''
    time=np.asarray(hk["MOUNT_GATTINI"]["Computer_UTC"])
    time_diff=np.diff(time)
    track=np.asarray(hk["MOUNT_GATTINI"]["tracking"]) #False during Meridian Flip, otherwise true if pointed to the source
    no_track=np.where(track=="False")[0]
    no_track=np.insert(no_track,[0,len(no_track)],[min(no_track)-1,max(no_track)+1]) #adding one more point before and after flip
    instert_val_num=int(time_after/np.median(time_diff))+1 #number of values to cut out after meridian flip to compensate for time to readjust the guiding
    no_track=np.append(no_track,np.arange(max(no_track),max(no_track)+1+instert_val_num,1))
    time_pd=pd.to_datetime(time, unit='s', utc=True).tz_convert('America/Los_Angeles')
    print("Start of meridian flip (local time): ",time_pd[no_track[0]],", End: ",time_pd[no_track[-1]])
    print("Duration of meridian flip cut: ",time[no_track[-1]]-time[no_track[0]])
    #Cutting out timestamps and data for meridian flip intervall
    data_cut=data[(timestamps<time[no_track[0]]) | (timestamps>time[no_track[-1]])]
    timestamps_cut=timestamps[(timestamps<time[no_track[0]]) | (timestamps>time[no_track[-1]])]
    print("Number of cut out events from meridian flip: ",len(data)-len(data_cut)," (",round((len(data)-len(data_cut))/len(data)*100,2),"%)")
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
    return(data_cut,timestamps_cut,time[no_track[0]],time[no_track[-1]])

#Example how to use the functions
file="start_2025-10-20T06-50-55Z.dp_ph1024.bpp_2.module_254.seqno_0.pff"
def main():
    data,metadata=read_pff(file)
    data, timestamps=cut_pkt_loss(data,metadata)
    data_cut,timestamps_cut=spike_cut(data,timestamps,plotting=True,path="")
    hk=read_pff_hk()
    data_cut_mer,timestamps_cut_mer,mflip_start,mflip_end=cut_meridian_flip(data_cut,timestamps_cut,hk,plotting=True,path="")