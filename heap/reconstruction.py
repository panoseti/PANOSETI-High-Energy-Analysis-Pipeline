"""reconstruction

Arrival direction and shower core reconstruction from the intersection of image axes, weighted by
size, elongation, and angle between axes (Hofmann et al. 1999, Astropart. Phys. 12, 135, method 1).

Data and simulations are reconstructed by the same functions. They differ only in the camera's axes
on the sky, which tilted_positions() takes as unit vectors: tilted_frame() for the simulations' alt-az
cameras, heap.events.sky_frame() for the data's cameras fixed to the sky.
"""
import numpy as np
import pandas as pd


def reconstruct_direction(event):
    """
    Reconstruct one event's arrival direction in camera coordinates.

    Parameters:
        event: one event's images (rows with x_c, y_c, phi, size, length, width), see
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

        s.append(tel["size"]) # not tel.size: that is the Series length
        x.append(tel.x_c)
        y.append(tel.y_c)
        phi_rad = tel.phi * np.pi/180
        m.append(np.tan(phi_rad))
        if(tel.length > 0):
            l.append(tel.width/tel.length)
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


def tilted_frame(pointing_alt, pointing_az):
    """
    Unit vectors of the plane perpendicular to the pointing, along the camera's axes, in CORSIKA's frame
    (x north, y west, z up).

    Parameters:
        pointing_alt, pointing_az: [deg], azimuth from north towards east

    Returns:
        (e_x, e_y): towards higher altitude and towards lower azimuth, as sim_telarray's camera +x and +y
        (see heap.process_simulation.camera_position())
    """
    alt, az = np.radians(pointing_alt), np.radians(pointing_az)
    e_x = np.array([-np.sin(alt) * np.cos(az), np.sin(alt) * np.sin(az), np.cos(alt)])
    e_y = np.array([np.sin(az), np.cos(az), 0.])
    return e_x, e_y


def tilted_positions(array, tel_pos, e_x, e_y):
    """
    Each image's telescope position in the plane perpendicular to the pointing.

    Parameters:
        array: images, see heap.events.apply_cuts()
        tel_pos: {telescope: position [m] in CORSIKA's frame (x north, y west, z up)}
        e_x, e_y: the camera's axes in that frame: tilted_frame() (simulations) or
            heap.events.sky_frame() (data)

    Returns:
        array of shape (n images, 2) [m], one row per row of array
    """
    pos = np.array([tel_pos[tel] for tel in array.Telescope], dtype=float).reshape(-1, 3)
    return np.stack([(pos * e_x).sum(axis=1), (pos * e_y).sum(axis=1)], axis=1)


def _least_squares_point(x, y, angle, w):
    """
    The point minimizing the w-weighted squared perpendicular distances to the lines through (x, y) at
    angle [rad]; None if the lines are all parallel.
    """
    n = np.stack([-np.sin(angle), np.cos(angle)], axis=1) # normals
    a = (w[:, None, None] * n[:, :, None] * n[:, None, :]).sum(axis=0)
    b = (w[:, None] * n * (n[:, 0] * x + n[:, 1] * y)[:, None]).sum(axis=0)
    if abs(np.linalg.det(a)) < 1e-9 * np.trace(a)**2:
        return None
    return np.linalg.solve(a, b)


def reconstruct_cores(array, positions):
    """
    Reconstruct each event's shower core in the plane of tilted_positions(): a shower's image points
    from its source position along the line from the telescope to the core, so the core is where the
    lines through the telescopes along their image axes (phi) cross: the point closest to all of them, each
    weighted by (size * (1 - width/length))^2.

    Parameters:
        array: images passing cuts, see heap.events.apply_cuts()
        positions: tilted_positions() of array

    Returns:
        DataFrame with Event, core_x, core_y [m], one row per reconstructed event (2+ images)
    """
    array = array.assign(tel_x=positions[:, 0], tel_y=positions[:, 1], angle=np.radians(array.phi.to_numpy()))
    rows = []
    for e, event in array.groupby("Event", sort=True):
        if len(event) < 2:
            continue
        w = (event["size"] * (1 - event.width / event.length))**2
        core = _least_squares_point(event.tel_x.to_numpy(), event.tel_y.to_numpy(), event.angle.to_numpy(), w.to_numpy())
        if core is not None:
            rows.append((e, *core))
    return pd.DataFrame(rows, columns=["Event", "core_x", "core_y"])


def core_distances(array, positions, cores):
    """
    Each image's telescope distance [m] from its event's reconstructed core, in the plane of
    tilted_positions(): the reconstructed impact distance. NaN for events without a core.

    Parameters:
        array: images, see heap.events.apply_cuts()
        positions: tilted_positions() of array
        cores: reconstruct_cores()

    Returns:
        array of distances, one per row of array
    """
    core = cores.set_index("Event").reindex(array.Event)
    return np.hypot(positions[:, 0] - core.core_x.to_numpy(dtype=float), positions[:, 1] - core.core_y.to_numpy(dtype=float))
