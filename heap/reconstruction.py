"""reconstruction

Arrival direction reconstruction from the intersection of image axes, weighted by size,
elongation, and angle between axes (Eventdisplay's method).
"""
import numpy as np
import pandas as pd


def reconstruct_direction(event):
    """
    Reconstruct one event's arrival direction in camera coordinates.

    Parameters:
        event: one event's images (rows with MeanX, MeanY, Phi, Size, Length, Width), see
            heap.events.apply_cuts()

    Returns:
        (Xoffset, Yoffset) in the images' units, or None if the image weights are invalid
    """
    NTel = len(event)

    m = []
    x = []
    y = []
    s = []
    l = []

    # get relevant image parameters
    for t in range(NTel):
        tel=event.iloc[t]

        s.append(tel.Size)
        x.append(tel.MeanX)
        y.append(tel.MeanY)
        phi_rad = tel.Phi * np.pi/180
        m.append(np.tan(phi_rad))
        if(tel.Length > 0):
            l.append(tel.Width/tel.Length)
        else:
            l.append(1)

    # don't do anything if angle between image axis is too small (for 2 images only)
    fiangdiff = 0
    if NTel == 2:
        fiangdiff = -1*np.abs(np.arctan(m[0]) - np.arctan(m[1])) * 180/np.pi

    # direction reconstruction
    itotweight = 0.
    iweight = 1.
    ixs = 0.
    iys = 0.
    iangdiff = 0.
    b1 = 0.
    b2 = 0.
    v_xs = []
    v_ys = []
    fmean_iangdiff = 0.
    fmean_iangdiffN = 0.

    for i in range(NTel):
        for j in range(NTel):
            if i >= j:
                continue

            # check minimum angle between image lines; ignore if too small
            iangdiff = np.abs( np.arctan( m[j] ) - np.arctan( m[i] ) )
            if( iangdiff < 0 or np.abs(np.pi - iangdiff ) < 0 ):
                continue

            # mean angle between images
            if( iangdiff < np.pi/2 ):
                fmean_iangdiff += iangdiff * 180/np.pi
            else:
                fmean_iangdiff += ( 180 - (iangdiff * 180/np.pi))
            fmean_iangdiffN += 1

            # weight is sin of angle between image lines
            iangdiff = np.abs( np.sin( np.abs( np.arctan( m[j] ) - np.arctan( m[i] ) ) ) )

            b1 = y[i] - m[i] * x[i]
            b2 = y[j] - m[j] * x[j]

            # line intersection
            if( m[i] != m[j] ):
                xs = ( b2 - b1 )  / ( m[i] - m[j] )
            else:
                xs = 0.

            ys = m[i] * xs + b1

            iweight  = 1. / ( 1. / s[i] + 1. / s[j] ) # weight 1: size of images
            iweight *= ( 1. - l[i] ) * ( 1. - l[j] ) # weight 2: elongation of images (width/length)
            iweight *= iangdiff                      # weight 3: angular differences between the two image axis
            iweight *= iweight                       # use squared value

            ixs += xs * iweight
            iys += ys * iweight
            itotweight += iweight

            v_xs.append( xs )
            v_ys.append( ys )

    # average difference between image pairs
    if( fmean_iangdiffN > 0. ):
        fmean_iangdiff /= fmean_iangdiffN
    else:
        fmean_iangdiff = 0.

    if( NTel > 2 ):
        fiangdiff = fmean_iangdiff

    # check validity of weight
    if( itotweight > 0. ):
        return ixs / itotweight, iys / itotweight
    return None


def reconstruct_directions(array):
    """
    reconstruct_direction() for every event with 2+ images.

    Parameters:
        array: images passing cuts, see heap.events.apply_cuts()

    Returns:
        DataFrame with Event, Date, Run, Xoffset, Yoffset, one row per reconstructed event
    """
    rows = []
    for e, event in array.groupby("Event", sort=True):
        # check that there are enough images
        if len(event) < 2:
            continue
        direction = reconstruct_direction(event)
        if direction is None:
            print("Image weights invalid")
            print(e)
            continue
        rows.append((e, event.iloc[0].Date, event.iloc[0].Run, *direction))

    return pd.DataFrame(rows, columns=["Event", "Date", "Run", "Xoffset", "Yoffset"])
