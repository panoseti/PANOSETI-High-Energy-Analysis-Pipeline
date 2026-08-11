#!/usr/bin/env python3

# pointing.py: Calculate pointing solutions from reference stars.

import numpy as np
import math


def walk_uphill(x0, y0, image, xgrid=None, ygrid=None):
    """
    Walk uphill from an initial point to find a local maximum in an image.

    The search starts at the pixel nearest to ``(x0, y0)`` and then moves to the
    neighboring pixel with the highest image value, continuing until no higher
    neighboring value is found.  The returned coordinates are the grid coordinates
    corresponding to the final pixel.

    Args:
        x0 (float): Starting X coordinate.
        y0 (float): Starting Y coordinate.
        image (array_like): 2D image array.
        xgrid (array_like, optional): X coordinates of the image columns. If not
            provided, a linear grid over the pixel indices is used.
        ygrid (array_like, optional): Y coordinates of the image rows. If not
            provided, a linear grid over the pixel indices is used.

    Returns:
        tuple: (x, y, value, ix, iy) where ``x`` and ``y`` are the coordinates of
        the maximum, ``value`` is the image value at that position, and ``ix``/``iy``
        are the corresponding pixel indices.
    """
    image = np.asarray(image)
    if image.ndim != 2:
        raise ValueError("image must be a 2D array")

    ny, nx = image.shape

    if xgrid is None:
        xgrid = np.linspace(0.0, nx - 1.0, nx)
    else:
        xgrid = np.asarray(xgrid, dtype=float)
        if xgrid.size != nx:
            raise ValueError("xgrid must have one entry per image column")

    if ygrid is None:
        ygrid = np.linspace(0.0, ny - 1.0, ny)
    else:
        ygrid = np.asarray(ygrid, dtype=float)
        if ygrid.size != ny:
            raise ValueError("ygrid must have one entry per image row")

    ix = int(np.argmin(np.abs(xgrid - x0)))
    iy = int(np.argmin(np.abs(ygrid - y0)))

    current_value = image[iy, ix]
    current_ix, current_iy = ix, iy

    while True:
        best_value = current_value
        best_ix, best_iy = current_ix, current_iy

        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                if dx == 0 and dy == 0:
                    continue

                nx_ix = current_ix + dx
                ny_iy = current_iy + dy
                if 0 <= nx_ix < nx and 0 <= ny_iy < ny:
                    candidate_value = image[ny_iy, nx_ix]
                    if candidate_value > best_value:
                        best_value = candidate_value
                        best_ix, best_iy = nx_ix, ny_iy

        if best_ix == current_ix and best_iy == current_iy:
            break

        current_ix, current_iy = best_ix, best_iy
        current_value = best_value

    return xgrid[current_ix], ygrid[current_iy], current_value, current_ix, current_iy


class PointingSolution:
    """
    Represents a pointing solution for translating between sky and image coordinates.
    """
    def __init__(self, ra0, dec0, plate_scale=1.0, theta=0.0, east_on_left=True, pos0=(0, 0), xy_units=None):
        """
        Initializes a pointing solution with explicit parameters.

        Args:
            ra0 (float): RA of the reference point (degrees).
            dec0 (float): Dec of the reference point (degrees).
            plate_scale (float): Plate scale in degrees/pixel (default 1.0).
            theta (float): Rotation angle in degrees (default 0.0).
            east_on_left (bool): If True, East is to the left (standard astronomical view).
            pos0 (tuple): Image coordinates (x, y) for the reference point (default (0, 0)).
        """
        self.ra0 = ra0
        self.dec0 = dec0
        self.plate_scale = plate_scale
        self.theta = math.radians((theta + 180.0) % 360.0 - 180.0)
        self.east_on_left = east_on_left
        self.x1, self.y1 = pos0
        if xy_units is None:
            units = { 1.0: 'deg', 
                      1.0/60.0: 'arcmin',
                      1.0/3600.0: 'arcsec',
                      180.0/np.pi: 'rad',
                      180.0/(1000.0*np.pi): 'mrad' }
            xy_units = units.get(plate_scale, 'pix')
        self.xy_units = xy_units

    @classmethod
    def solve(cls, star1, pos1, star2, pos2, east_on_left=True, plate_scale=None, ref_pos=(0, 0), ref_sky=None):
        """
        Calculates a pointing solution from two stars.

        Args:
            star1 (dict): Reference star 1 with 'ra_deg' and 'dec_deg'.
            pos1 (tuple): Image coordinates (x, y) for star 1.
            star2 (dict): Reference star 2 with 'ra_deg' and 'dec_deg'.
            pos2 (tuple): Image coordinates (x, y) for star 2.
            east_on_left (bool): If True, East is to the left.
            plate_scale (float, optional): If provided, force this plate scale (deg/pix).
                                           The distance between pos1 and pos2 will be
                                           ignored for scale but used for orientation.
            ref_pos (tuple, optional): Image coordinates (x, y) to set as the reference point. Default is (0, 0).
            ref_sky (tuple, optional): Sky coordinates (ra, dec) to set as the reference point.
                                       If provided, overrides ref_pos.
        """
        # We'll use star1 as the provisional reference point
        ps = cls(star1['ra_deg'], star1['dec_deg'], east_on_left=east_on_left, pos0=pos1)

        dx = pos2[0] - pos1[0]
        dy = pos2[1] - pos1[1]
        d_pix = math.sqrt(dx**2 + dy**2)

        if d_pix == 0:
            raise ValueError("Star positions in image must be distinct.")

        f = -1.0 if east_on_left else 1.0
        phi_img = math.atan2(dy, f * dx)

        def compute_params(ra0, dec0):
            ps.ra0 = ra0
            ps.dec0 = dec0
            xi1, eta1 = ps._sky_to_tangent(star1['ra_deg'], star1['dec_deg'])
            xi2, eta2 = ps._sky_to_tangent(star2['ra_deg'], star2['dec_deg'])
            
            dxi = xi2 - xi1
            deta = eta2 - eta1
            d_tp = math.sqrt(dxi**2 + deta**2)
            
            scale = plate_scale if plate_scale is not None else d_tp / d_pix
            
            phi_tp = math.atan2(deta, dxi)
            theta = phi_tp - phi_img
            ps.theta = (theta + math.pi) % (2 * math.pi) - math.pi
            ps.plate_scale = scale
            
            xi1_s = xi1 / scale
            eta1_s = eta1 / scale
            cos_th = math.cos(ps.theta)
            sin_th = math.sin(ps.theta)
            
            x0 = pos1[0] - f * (xi1_s * cos_th + eta1_s * sin_th)
            y0 = pos1[1] - (-xi1_s * sin_th + eta1_s * cos_th)
            
            ps.x1 = x0
            ps.y1 = y0

        # Initial solve with star1 as the tangent point
        compute_params(star1['ra_deg'], star1['dec_deg'])

        if ref_sky is not None:
            compute_params(ref_sky[0], ref_sky[1])
            return ps

        if ref_pos is not None:
            target_x, target_y = ref_pos
            ra0, dec0 = ps.image_to_sky(target_x, target_y)
            for _ in range(20):
                compute_params(ra0, dec0)
                if abs(ps.x1 - target_x) < 1e-10 and abs(ps.y1 - target_y) < 1e-10:
                    break
                ra0, dec0 = ps.image_to_sky(target_x, target_y)

        return ps

    def _sky_to_tangent(self, ra, dec):
        """Convert RA/Dec to tangent plane coordinates (degrees) relative to (ra0, dec0)."""
        ra_rad = math.radians(ra)
        dec_rad = math.radians(dec)
        ra0_rad = math.radians(self.ra0)
        dec0_rad = math.radians(self.dec0)

        cos_d = math.cos(dec_rad)
        sin_d = math.sin(dec_rad)
        cos_d0 = math.cos(dec0_rad)
        sin_d0 = math.sin(dec0_rad)
        cos_dra = math.cos(ra_rad - ra0_rad)
        sin_dra = math.sin(ra_rad - ra0_rad)

        denom = sin_d * sin_d0 + cos_d * cos_d0 * cos_dra
        xi = cos_d * sin_dra / denom
        eta = (sin_d * cos_d0 - cos_d * sin_d0 * cos_dra) / denom

        return math.degrees(xi), math.degrees(eta)

    def _tangent_to_sky(self, xi, eta):
        """Convert tangent plane coordinates (degrees) back to RA/Dec."""
        xi_rad = math.radians(xi)
        eta_rad = math.radians(eta)
        ra0_rad = math.radians(self.ra0)
        dec0_rad = math.radians(self.dec0)

        rho = math.sqrt(xi_rad**2 + eta_rad**2)
        if rho == 0:
            return self.ra0, self.dec0

        c = math.atan(rho)
        cos_c = math.cos(c)
        sin_c = math.sin(c)
        cos_dec0 = math.cos(dec0_rad)
        sin_dec0 = math.sin(dec0_rad)

        dec = math.asin(cos_c * sin_dec0 + (eta_rad * sin_c * cos_dec0) / rho)
        ra = ra0_rad + math.atan2(xi_rad * sin_c, rho * cos_dec0 * cos_c - eta_rad * sin_dec0 * sin_c)

        return math.degrees(ra) % 360, math.degrees(dec)

    def sky_to_image(self, ra, dec):
        """Converts RA/Dec to image coordinates (x, y)."""
        xi, eta = self._sky_to_tangent(ra, dec)

        xi_s = xi / self.plate_scale
        eta_s = eta / self.plate_scale

        cos_th = math.cos(self.theta)
        sin_th = math.sin(self.theta)

        # Rotate back:
        x_p = xi_s * cos_th + eta_s * sin_th
        dy = -xi_s * sin_th + eta_s * cos_th

        f = -1.0 if self.east_on_left else 1.0
        dx = x_p * f

        return self.x1 + dx, self.y1 + dy

    def image_to_sky(self, x, y):
        """Converts image coordinates (x, y) to RA/Dec."""
        dx = x - self.x1
        dy = y - self.y1
        f = -1.0 if self.east_on_left else 1.0

        cos_th = math.cos(self.theta)
        sin_th = math.sin(self.theta)

        xi = self.plate_scale * (f * dx * cos_th - dy * sin_th)
        eta = self.plate_scale * (f * dx * sin_th + dy * cos_th)

        return self._tangent_to_sky(xi, eta)

    def print_stars(self, stars):
        """
        Prints a formatted table of stars including their image coordinates.

        Args:
            stars (list): List of star dictionaries from get_bright_stars().
        """
        from .stars import star_table
        table = star_table(stars)
        if not stars:
            print(table)
            return

        lines = table.splitlines()

        # Extend header
        lines[0] += f" {'X [' + self.xy_units[:6] + ']' if self.xy_units != ''else 'X':>10}"
        lines[0] += f" {'Y [' + self.xy_units[:6] + ']' if self.xy_units != '' else 'Y':>10}"
        lines[1] += "-" * 22

        # Extend data rows
        for i, star in enumerate(stars):
            x, y = self.sky_to_image(star['ra_deg'], star['dec_deg'])
            lines[i+2] += f" {x:10.2f} {y:10.2f}"

        print("\n".join(lines))

    def reset_reference_point(self, pos=None, sky=None):
        """
        Returns a new PointingSolution with the reference point reset to a 
        new image position (pos) or sky coordinates (sky). All translations 
        between sky and image will be the same (to leading order).
        
        Args:
            pos (tuple, optional): (x, y) image coordinates.
            sky (tuple, optional): (ra, dec) sky coordinates in degrees.
        """
        if pos is not None and sky is not None:
            raise ValueError("Provide either pos or sky, not both.")
        
        if pos is not None:
            new_x, new_y = pos
            new_ra, new_dec = self.image_to_sky(new_x, new_y)
        elif sky is not None:
            new_ra, new_dec = sky
            new_x, new_y = self.sky_to_image(new_ra, new_dec)
        else:
            raise ValueError("Must provide either pos or sky.")
            
        pos2 = (new_x, new_y + 1.0)
        ra2, dec2 = self.image_to_sky(*pos2)
        
        new_ps = PointingSolution.solve(
            star1={'ra_deg': new_ra, 'dec_deg': new_dec},
            pos1=(new_x, new_y),
            star2={'ra_deg': ra2, 'dec_deg': dec2},
            pos2=pos2,
            east_on_left=self.east_on_left,
            plate_scale=self.plate_scale,
            ref_pos=None
        )
        new_ps.xy_units = self.xy_units
        return new_ps

    def shift_origin(self, dx, dy):
        """
        Returns a new PointingSolution shifted by (dx, dy).
        
        Args:
            dx (float): Shift in X image coordinates.
            dy (float): Shift in Y image coordinates.
        """
        return PointingSolution(
            self.ra0, self.dec0, self.plate_scale, math.degrees(self.theta),
            self.east_on_left, (self.x1 + dx, self.y1 + dy), self.xy_units
        )

    def rotate(self, angle_deg, origin=(0, 0)):
        """
        Returns a new PointingSolution rotated by angle_deg around origin.
        
        Args:
            angle_deg (float): Rotation angle in degrees (counter-clockwise).
            origin (tuple): Image coordinates (x, y) around which to rotate (default (0, 0)).
        """
        angle_rad = math.radians(angle_deg)
        cos_a = math.cos(angle_rad)
        sin_a = math.sin(angle_rad)
        
        # Rotate the reference point (x1, y1) around the origin
        ox, oy = origin
        dx = self.x1 - ox
        dy = self.y1 - oy
        
        new_x1 = ox + dx * cos_a - dy * sin_a
        new_y1 = oy + dx * sin_a + dy * cos_a
        
        # New rotation angle
        new_theta_deg = math.degrees(self.theta) + angle_deg
        
        return PointingSolution(
            self.ra0, self.dec0, self.plate_scale, new_theta_deg,
            self.east_on_left, (new_x1, new_y1), self.xy_units
        )

    def zoom(self, factor, origin=(0, 0)):
        """
        Returns a new PointingSolution zoomed by factor around origin.
        
        Args:
            factor (float): Zoom factor (> 1 is zoom in, < 1 is zoom out).
                            Image coordinates are multiplied by this factor.
            origin (tuple): Image coordinates (x, y) around which to zoom (default (0, 0)).
        """
        # New plate scale (deg/pix). 
        # Zooming in (factor > 1) means fewer degrees per pixel.
        new_scale = self.plate_scale / factor
        
        # Scale the reference point (x1, y1) relative to origin
        ox, oy = origin
        new_x1 = ox + (self.x1 - ox) * factor
        new_y1 = oy + (self.y1 - oy) * factor
        
        return PointingSolution(
            self.ra0, self.dec0, new_scale, math.degrees(self.theta),
            self.east_on_left, (new_x1, new_y1), self.xy_units
        )

    def __repr__(self):
        return (f"<PointingSolution RA0={self.ra0:.4f} Dec0={self.dec0:.4f} "
                f"scale={self.plate_scale:.6f} deg/pix theta={math.degrees(self.theta):.2f}deg "
                f"east_on_left={self.east_on_left}>")

def hms_to_deg(hms):
    """
    Convert an angle string in hours/minutes/seconds to degrees.

    Accepts separators spaces, colons, or h/m/s notation.
    Examples:
        "19 59 59.8520303688"
        "19:59:59.8520303688"
        "19h59m59.8520303688s"
    """
    if isinstance(hms, (int, float)):
        return float(hms)

    text = str(hms).strip().lower()
    if not text:
        raise ValueError("Empty HMS string")

    text = text.replace("h", " ").replace("m", " ").replace("s", " ")
    text = text.replace(":", " ")
    parts = [p for p in text.split() if p]

    if not parts:
        raise ValueError(f"Invalid HMS value: {hms!r}")

    if len(parts) == 1:
        hours = float(parts[0])
        minutes = 0.0
        seconds = 0.0
    elif len(parts) == 2:
        hours = float(parts[0])
        minutes = float(parts[1])
        seconds = 0.0
    else:
        hours = float(parts[0])
        minutes = float(parts[1])
        seconds = float(parts[2])

    sign = -1.0 if hours < 0 else 1.0
    hours = abs(hours)
    return sign * (hours + minutes / 60.0 + seconds / 3600.0) * 15.0

def dms_to_deg(dms):
    """
    Convert an angle string in degrees/minutes/seconds to decimal degrees.

    Accepts separators spaces, colons, or d/m/s or degree symbols.
    Examples:
        "19 59 59.8520303688"
        "19:59:59.8520303688"
        "19d59m59.8520303688s"
        "-19 59 59.85"
    """
    if isinstance(dms, (int, float)):
        return float(dms)

    text = str(dms).strip().lower()
    if not text:
        raise ValueError("Empty DMS string")

    # Normalize various symbols to spaces
    for ch in ['d', '°', "'", 'm', '"', 's', ':']:
        text = text.replace(ch, ' ')

    parts = [p for p in text.split() if p]
    if not parts:
        raise ValueError(f"Invalid DMS value: {dms!r}")

    if len(parts) == 1:
        deg = float(parts[0])
        minutes = 0.0
        seconds = 0.0
    elif len(parts) == 2:
        deg = float(parts[0])
        minutes = float(parts[1])
        seconds = 0.0
    else:
        deg = float(parts[0])
        minutes = float(parts[1])
        seconds = float(parts[2])

    sign = -1.0 if str(parts[0]).lstrip().startswith('-') or deg < 0 else 1.0
    deg = abs(deg)
    return sign * (deg + minutes / 60.0 + seconds / 3600.0)

OBSERVATION_TARGETS = {
    "crab":            ( 83.63308,  22.01450),
    "mrk421":          (166.11379,  38.20883),
    "mrk501":          (253.46758,  39.76017),
    "galactic_center": (266.41683, -29.00781),
    "hessj1837-069":   (279.41500,  -6.95000),
    "1es1959+650":     (299.99938,  65.14851),
    "mgroj2019+37":    (304.646  ,  36.847),
}

def get_target_coordinates(name):
    """
    Returns the RA and Dec of an observation target in degrees.

    Args:
        name (str): The name/key of the target (e.g., 'crab', 'mrk421').

    Returns:
        tuple: (RA, Dec) in decimal degrees, or None if not found.
    """
    return OBSERVATION_TARGETS.get(name.lower())

