"""significance

On/off region counting, Li & Ma significance, and sky maps for reconstructed arrival directions
(see heap.reconstruction).

Each run is counted in its own camera coordinates (degrees from that run's pointing, see make_wcs()
and heap.events.run_pointings()), so runs with different wobble offsets combine correctly: a sky
position is converted to every run's camera coordinates, and off regions are placed around it there.
Only regions entirely inside the camera are counted (in_camera()).
"""
import numpy as np
import pandas as pd
from astropy import wcs
from scipy.spatial import cKDTree

from heap.parameterize import CAMERA_HALF_WIDTH


def make_wcs(center):
    """
    Tangent-plane WCS centered on the pointing, with 1 deg "pixels" so camera coordinates in
    degrees convert directly: w.wcs_pix2world(Xoffset, Yoffset, 1).

    Camera coordinates are in the preflip orientation (heap.events.load_camera_frame() rotates
    postflip images 180 deg to it), where +x (columns) points east and +y (rows) points south. Verified against
    catalog star positions in pedvar maps on the Sept 2026 nights, all four telescopes.

    Parameters:
        center: SkyCoord of the pointing
    """
    w = wcs.WCS(naxis=2)
    w.wcs.crpix = [0, 0]
    w.wcs.crval = [center.ra.deg, center.dec.deg]
    w.wcs.cdelt = [1, -1] # +x east, +y south
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    return w


def camera_to_sky(df, pointings, x="Xoffset", y="Yoffset"):
    """
    RA/DEC (deg) of camera positions, each row converted with its own Run's pointing.

    Parameters:
        df: rows with Run and camera coordinates x, y (deg)
        pointings: {Run: SkyCoord} each run's pointing, see heap.events.run_pointings()

    Returns:
        (ra, dec) arrays aligned with df
    """
    ra = np.full(len(df), np.nan)
    dec = np.full(len(df), np.nan)
    for run, idx in df.groupby("Run").indices.items():
        ra[idx], dec[idx] = make_wcs(pointings[run]).wcs_pix2world(df[x].to_numpy()[idx], df[y].to_numpy()[idx], 1)
    return ra, dec


class RegionCounter:
    """
    Counts one pointing's events in circular regions, in its camera coordinates. An event counts in a
    region if its reconstructed direction is within theta of the region center, at least one of its
    images passes position_cuts, and every such image's centroid is within max_distance of the
    center.

    Parameters:
        directions: this pointing's reconstructed events, see heap.reconstruction.reconstruct_directions()
        array: images passing cuts, see heap.events.apply_cuts()
        w: this pointing's WCS from make_wcs()
        theta: region radius (deg)
        max_distance: max distance cut (deg); None for no cut
        position_cuts: optional {column: (mode, threshold)} cuts on distance, alpha, and/or miss
            recomputed against the region center, applied in order (see heap.events.apply_cut());
            "nsigma" statistics are over this pointing's images of the events within theta
    """

    def __init__(self, directions, array, w, theta, max_distance, position_cuts=None):
        self.directions = directions
        self.w = w
        self.theta = theta
        self.max_distance = max_distance
        self.position_cuts = position_cuts or {}
        self.x = directions.Xoffset.to_numpy()
        self.y = directions.Yoffset.to_numpy()
        self.tree = cKDTree(np.column_stack([self.x, self.y]))

        # each event's images, padded with nan (Telescope with -1): (n_events, max images per event)
        img = array[array.Event.isin(directions.Event)]
        row = img.Event.map(pd.Series(np.arange(len(directions)), index=directions.Event)).to_numpy()
        col = img.groupby("Event").cumcount().to_numpy()
        shape = (len(directions), col.max() + 1)
        self.img_x = np.full(shape, np.nan)
        self.img_y = np.full(shape, np.nan)
        self.img_phi = np.full(shape, np.nan)
        self.img_tel = np.full(shape, -1)
        self.img_x[row, col] = img.x_c
        self.img_y[row, col] = img.y_c
        self.img_phi[row, col] = img.phi
        self.img_tel[row, col] = pd.factorize(img.Telescope)[0]
        self.img_index = np.full(shape, -1)
        self.img_index[row, col] = np.arange(len(img))
        self.img = img

    def to_camera(self, ra, dec):
        """Sky position(s) (deg) in this pointing's camera coordinates (deg)."""
        return self.w.wcs_world2pix(ra, dec, 1)

    def shower_params(self, center_x, center_y, idx=slice(None)):
        """distance, alpha, miss of events idx's images relative to (center_x, center_y)."""
        dx = self.img_x[idx] - center_x
        dy = self.img_y[idx] - center_y

        distance = np.hypot(dx, dy)
        psi = np.degrees(np.arctan2(dy, dx))
        diff = (self.img_phi[idx] - psi) % 180
        alpha = np.where(diff <= 90, diff, 180 - diff)
        miss = distance * np.sin(np.radians(alpha))

        return {"distance": distance, "alpha": alpha, "miss": miss}

    def passing_images(self, center_x, center_y, idx=slice(None)):
        """Mask of events idx's images passing position_cuts relative to (center_x, center_y), and
        their shower_params()."""
        params = self.shower_params(center_x, center_y, idx)
        tel = self.img_tel[idx]
        keep = tel >= 0
        for column, (mode, threshold) in self.position_cuts.items():
            values = params[column]
            if mode == "nsigma":
                cut = np.full(values.shape, np.nan)
                for t in np.unique(tel[keep]):
                    sel = keep & (tel == t)
                    std = values[sel].std(ddof=1) if sel.sum() > 1 else np.nan # pandas' std
                    cut[tel == t] = values[sel].mean() + threshold*std
            elif mode == "value":
                cut = threshold
            else:
                raise ValueError(f"Unknown cut mode for {column}: {mode!r}")
            keep &= values < cut
        return keep, params

    def max_image_distance(self, center_x, center_y, idx=slice(None)):
        """Per-event distance of the farthest image centroid passing position_cuts from
        (center_x, center_y); nan for events with no such image."""
        keep, params = self.passing_images(center_x, center_y, idx)
        distance = np.where(keep, params["distance"], -np.inf).max(axis=1)
        return np.where(keep.any(axis=1), distance, np.nan)

    def passing_events(self, center_x, center_y, idx=slice(None)):
        """Mask of events idx passing position_cuts and the max distance cut."""
        distance = self.max_image_distance(center_x, center_y, idx)
        if self.max_distance is None:
            return ~np.isnan(distance)
        return distance < self.max_distance

    def in_region(self, center_x, center_y):
        """Row indices into directions of events counted in the region at camera (center_x, center_y)."""
        idx = np.array(self.tree.query_ball_point((center_x, center_y), self.theta), dtype=int)
        angle = np.hypot(self.x[idx] - center_x, self.y[idx] - center_y)
        idx = idx[angle*angle < self.theta**2]
        return idx[self.passing_events(center_x, center_y, idx)]


def make_reflected_regions(x_on, y_on, theta):
    """
    Reflected off regions: the on region rotated about the camera center (the pointing) in equal
    steps, as many as fit without overlapping each other or the on region. Every off region sits
    at the on region's offset from the pointing, so has the same acceptance.

    Parameters:
        x_on, y_on: on region center in camera coordinates (deg)
        theta: region radius (deg)

    Returns:
        list of camera (x, y) region centers; empty if the on region is within theta of the pointing
    """
    r = np.hypot(x_on, y_on)
    if r <= theta:
        return []
    n = int(np.floor(np.pi / np.arcsin(theta / r))) # regions (on included) around the circle of radius r
    angles = np.arctan2(y_on, x_on) + 2*np.pi*np.arange(1, n)/n
    return list(zip(r*np.cos(angles), r*np.sin(angles)))


def make_off_regions(testPosX, testPosY, theta):
    """
    Rings of off regions of radius theta around the test position, restricted to its vicinity.

    Returns:
        list of region centers, in the test position's coordinates
    """
    off_regions=[]
    spacing = 2*theta
    buffer = 4*theta
    limit = 4*theta

    r_start = max(buffer, spacing)  # first ring can't be closer than one circle-width out
    n_rings = int(np.floor((limit - r_start) / spacing)) + 1

    for k in range(n_rings):
        r = r_start + k*spacing
        n_points = max(1, int(np.floor(2*np.pi*r / spacing)))
        for p in range(n_points):
            angle = 2*np.pi*p/n_points
            center_x = testPosX + r*np.cos(angle)
            center_y = testPosY + r*np.sin(angle)
            off_regions.append((center_x,center_y))

    return off_regions


OFF_REGION_METHODS = {"reflected": make_reflected_regions, "ring": make_off_regions}


def in_camera(x, y, theta):
    """True if the region of radius theta centered at camera (x, y) (deg) lies entirely inside the camera."""
    return max(abs(x), abs(y)) + theta <= CAMERA_HALF_WIDTH


def calc_alpha(off_regions, theta):
    """Ratio of on region area to total off region area; 0 if there are no off regions."""
    if len(off_regions) == 0:
        return 0.
    area_B = np.pi*theta**2
    area_A = len(off_regions) * np.pi * theta**2
    alpha=area_B/area_A
    return alpha


def combine_on_off(counts):
    """
    Totals over pointings with different alphas (as gammapy's MapDatasetOnOff.stack()): N_on and N_off summed,
    alpha = sum(alpha_i*N_off_i)/sum(N_off_i). Pointings without off regions (alpha 0) are left out.

    Parameters:
        counts: DataFrame with on, off, alpha, e.g. OnOffCounter.per_run()

    Returns:
        (N_on, N_off, alpha)
    """
    counts = counts[counts.alpha > 0]
    if len(counts) == 0:
        return 0, 0, 0.
    N_on = counts.on.sum()
    N_off = counts.off.sum()
    alpha = (counts.alpha*counts.off).sum()/N_off if N_off > 0 else counts.alpha.mean()
    return N_on, N_off, alpha


class OnOffCounter:
    """
    On/off counts at a sky position over several pointings: each Run's events are counted in its
    own camera coordinates (RegionCounter), with off regions placed around the on region there.

    Parameters:
        directions: reconstructed events with Date and Run, see heap.reconstruction.reconstruct_directions()
        array: images passing cuts, see heap.events.apply_cuts()
        pointings: {Run: SkyCoord} each run's pointing, see heap.events.run_pointings()
        theta: region radius (deg)
        max_distance: max distance cut (deg); None for no cut
        off_method: "reflected" (around each pointing, default) or "ring" (around the on region),
            see make_reflected_regions()/make_off_regions(); only off regions entirely inside the
            camera are counted and make up alpha, see off_regions()
        position_cuts: optional {column: (mode, threshold)} cuts on distance, alpha, and/or miss
            relative to each region center, see RegionCounter
    """

    def __init__(self, directions, array, pointings, theta, max_distance, off_method="reflected", position_cuts=None):
        self.theta = theta
        self.make_off_regions = OFF_REGION_METHODS[off_method]
        self.counters = {
            run: RegionCounter(events, array[array.Date == events.Date.iloc[0]], make_wcs(pointings[run]), theta, max_distance, position_cuts)
            for run, events in directions.groupby("Run")
        }

    def off_regions(self, x_on, y_on):
        """Camera (x, y) centers of off_method's off regions for the on region at camera (x_on, y_on)
        that lie entirely inside the camera (in_camera()); none if the on region doesn't, so that Run
        is left out (alpha 0)."""
        if not in_camera(x_on, y_on, self.theta):
            return []
        return [region for region in self.make_off_regions(x_on, y_on, self.theta) if in_camera(*region, self.theta)]

    def _on_off(self, counter, x_on, y_on):
        off_regions = self.off_regions(x_on, y_on)
        on = counter.in_region(x_on, y_on)
        off = [counter.in_region(*region) for region in off_regions]
        return on, off, calc_alpha(off_regions, self.theta)

    def per_run(self, ra, dec):
        """On/off counts and alpha at (ra, dec) for each Run: DataFrame with Date, Run, on, off, alpha."""
        rows = []
        for run, counter in self.counters.items():
            on, off, alpha = self._on_off(counter, *counter.to_camera(ra, dec))
            rows.append((counter.directions.Date.iloc[0], run, len(on), sum(len(o) for o in off), alpha))
        return pd.DataFrame(rows, columns=["Date", "Run", "on", "off", "alpha"])

    def per_date(self, ra, dec):
        """per_run() combined over each Date's runs (combine_on_off()): DataFrame with Date, on, off, alpha."""
        rows = [(date, *combine_on_off(runs)) for date, runs in self.per_run(ra, dec).groupby("Date")]
        return pd.DataFrame(rows, columns=["Date", "on", "off", "alpha"])

    def count(self, ra, dec):
        """(N_on, N_off, alpha) at (ra, dec), combined over Runs (combine_on_off())."""
        return combine_on_off(self.per_run(ra, dec))

    def events(self, ra, dec):
        """(on, off) events at (ra, dec); an event is counted once per off region it falls in."""
        on_parts, off_parts = [], []
        for counter in self.counters.values():
            on, off, _ = self._on_off(counter, *counter.to_camera(ra, dec))
            on_parts.append(counter.directions.iloc[on])
            off_parts += [counter.directions.iloc[o] for o in off]
        return pd.concat(on_parts), pd.concat(off_parts)

    def off_region_centers(self, ra, dec):
        """Sky (RA, DEC) centers of the off regions for the on region at (ra, dec), each listed once:
        Runs with the same pointing share them."""
        centers = []
        for counter in self.counters.values():
            regions = self.off_regions(*counter.to_camera(ra, dec))
            if regions:
                centers += [tuple(c) for c in counter.w.wcs_pix2world(np.array(regions), 1)]
        return list(dict.fromkeys(centers))

    def theta_square(self, ra, dec):
        """Each event's squared angular distance from (ra, dec) and max image distance (among
        images passing position_cuts, nan if none do), in its own camera coordinates: DataFrame with
        Event, ThetaSquare, Distance. Does not apply the max distance cut."""
        parts = []
        for counter in self.counters.values():
            x_on, y_on = counter.to_camera(ra, dec)
            parts.append(pd.DataFrame({
                "Event": counter.directions.Event.to_numpy(),
                "ThetaSquare": np.hypot(counter.x - x_on, counter.y - y_on)**2,
                "Distance": counter.max_image_distance(x_on, y_on),
            }))
        return pd.concat(parts, ignore_index=True)

    def images(self, ra, dec):
        """Images passing position_cuts relative to (ra, dec), of events also passing the max
        distance cut, with distance, alpha, miss recomputed relative to (ra, dec) in each image's
        camera coordinates. No theta cut."""
        parts = []
        for counter in self.counters.values():
            x_on, y_on = counter.to_camera(ra, dec)
            keep, params = counter.passing_images(x_on, y_on)
            keep &= counter.passing_events(x_on, y_on)[:, None]
            img = counter.img.iloc[counter.img_index[keep]].copy()
            for column, values in params.items():
                img[column] = values[keep]
            parts.append(img)
        return pd.concat(parts, ignore_index=True)


# "borrowed" from Eventdisplay VStatistics.h
def significance(Non,Noff,alpha):
    """Li & Ma significance, Equation 17."""

    if( alpha == 0. ):
        Sig17 = 0.0
        return Sig17

    alphasq = alpha * alpha
    oneplusalpha = 1.0 + alpha
    oneplusalphaoveralpha = oneplusalpha / alpha

    Nsig = Non - alpha * Noff
    Ntot = Non + Noff

    # Sig5
    if( Non + alphasq * Noff > 0. ):
        Sig5 = Nsig / np.sqrt( Non + alphasq* Noff )
    else:
        Sig5 = 0.

    # Sig9
    if( alpha * Ntot > 0. ):
        Sig9 = Nsig / np.sqrt( alpha* Ntot )
    else:
        Sig9 = 0.

    # Sig17
    if( Ntot == 0. ):
        Sig17 = 0.
    elif( Non == 0 and Noff != 0. ):
        Sig17 = np.sqrt( 2.*( Noff* np.log( oneplusalpha * ( Noff / Ntot ) ) ) )
    elif( Non != 0 and Noff == 0. ):
        Sig17 = np.sqrt( 2.*( Non* np.log( oneplusalphaoveralpha * ( Non / Ntot ) ) ) )
    else:
        Sig17 = 2.*( Non* np.log( oneplusalphaoveralpha * ( Non / Ntot ) ) + Noff* np.log( oneplusalpha * ( Noff / Ntot ) ) )
        # value in brackets can be a small negative number
        if Sig17 > 0.:
            Sig17 = np.sqrt( Sig17 )
        else:
            Sig17 = 0.


    if( Nsig < 0 ):
        Sig17 = -Sig17

    # return Sig5
    # return Sig9
    return Sig17


def calc_sigmap(counter, bins_RA, bins_DEC, bin_width):
    """
    On/off counts and significance at every test position of the sky map.

    Parameters:
        counter: OnOffCounter
        bins_RA, bins_DEC: sky map bin edges; test positions are edge + bin_width/2
        bin_width: sky map bin width (deg)

    Returns:
        DataFrame with testPosX, testPosY, significance, on, off, alpha
    """
    rows = []

    for x in bins_RA:
        x=x+bin_width/2 #bin center
        for y in bins_DEC:
            y=y+bin_width/2 #bin center

            N_on, N_off, alpha = counter.count(x, y)

            sigma=significance(N_on,N_off,alpha)
            rows.append((x,y,sigma,N_on,N_off,alpha))

    return pd.DataFrame(rows, columns=['testPosX','testPosY','significance','on','off','alpha'])


def plot_sky_map(ax, sig, weights, label, bins_RA, bins_DEC, on_region=None, vmin=None, vmax=None):
    """
    Plots one calc_sigmap() quantity as a sky map (RA increasing to the left).

    Parameters:
        ax: matplotlib Axes
        sig: DataFrame from calc_sigmap()
        weights: per-test-position values, e.g. sig.significance
        label: colorbar label
        on_region: optional (RA, DEC) to mark with a +
        vmin, vmax: optional colorbar range

    Returns:
        the hist2d return value
    """
    hist = ax.hist2d(sig.testPosX, sig.testPosY, weights=weights, bins=[bins_RA, bins_DEC], cmap='viridis', vmin=vmin, vmax=vmax)

    ax.invert_xaxis()
    ax.set_aspect('equal')
    ax.set_xlabel('RA (degrees)')
    ax.set_ylabel('DEC (degrees)')

    ax.figure.colorbar(hist[3], ax=ax, label=label)
    if on_region is not None:
        ax.scatter(*on_region, marker='+', color='white', s=150)
    return hist
