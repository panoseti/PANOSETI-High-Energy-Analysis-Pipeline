"""read_pcap

Functions for reading pcap data into numpy arrays
"""

import numpy as np
import uproot

def read_tel_data(data_file):
    """
    Read the data from a .root file derived from .pcap data of a single telescope and return a numpy array

    Parameters:
        data_file: .root file containing a tree named "camdata" and branch "cam_2D_hist"
        
    Returns:
        data: numpy array with shape (N, 32, 32)
    """

    data = []
    with uproot.open(data_file) as f:
        tree = f["camdata"]
        branch = tree["cam_2D_hist"]

        for event in branch.array(library="np"):
            image = event.values(flow=False)
            data.append(image.astype(np.float32))

    data = np.stack(data)
    return data

def read_tel_peds(pedestal_file):
    """
    Read the data from a .root file derived from .pcap data of a single telescope and return a numpy array

    Parameters:
        pedestal_file: .root file containing a histogram named "peds_2D_hist;1"
        
    Returns:
        pedvars: numpy array with shape (32, 32)
    """

    with uproot.open(pedestal_file) as f:
        h = f["peds_2D_hist;1"]
        image = h.values(flow=False).astype(np.float32)
        
        return image

def read_tel_pedvars(pedvar_file):
    """
    Read the data from a .root file derived from .pcap data of a single telescope and return a numpy array

    Parameters:
        pedvar_file: .root file containing a histogram named "pedvars_2D_hist;1"
        
    Returns:
        pedvars: numpy array with shape (32, 32)
    """

    with uproot.open(pedvar_file) as f:
        h = f["pedvars_2D_hist;1"]
        image = h.values(flow=False).astype(np.float32)
    
    return image