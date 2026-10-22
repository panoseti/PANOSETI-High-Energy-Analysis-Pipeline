"""read_pff

Functions for reading .pff data into numpy arrays
"""

import pypff

def read_pff(filename):
    #Reads in pff file and returns data and metadata
    pff = pypff.io.datapff(filename)
    data, metadata = pff.readpff(metadata=True)
    print("Data read in, number of events: ",len(data))
    return(data,metadata)

def read_pff_hk(filename="hk.pff"):
    #Read in housekeeping pff file and return data, default called hk.pff, gonna need to adjust path
    hkpff = pypff.io.hkpff(filename)
    hk = hkpff.readhk()
    return(hk)
