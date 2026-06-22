#!/usr/bin/env python3

# pointing.py: Calculate pointing solutions from reference stars.

import numpy as np
import math

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
        self.theta = math.radians(theta)
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
    def solve(cls, star1, pos1, star2, pos2, east_on_left=True, plate_scale=None):
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
        """
        # Create temporary instance to use _sky_to_tangent (which depends on self.ra0/dec0)
        # We'll use star1 as the reference point
        ps = cls(star1['ra_deg'], star1['dec_deg'], east_on_left=east_on_left, pos0=pos1)

        # Tangent plane coordinates of star2 relative to star1
        xi2, eta2 = ps._sky_to_tangent(star2['ra_deg'], star2['dec_deg'])
        d_tp = math.sqrt(xi2**2 + eta2**2)

        # Image coordinates of star2 relative to star1
        dx = pos2[0] - pos1[0]
        dy = pos2[1] - pos1[1]
        d_pix = math.sqrt(dx**2 + dy**2)

        if d_pix == 0:
            raise ValueError("Star positions in image must be distinct.")

        if plate_scale is not None:
            ps.plate_scale = plate_scale
        else:
            ps.plate_scale = d_tp / d_pix

        # Parity factor
        f = -1.0 if east_on_left else 1.0

        # Calculate rotation angle theta
        phi_img = math.atan2(dy, f * dx)
        phi_tp = math.atan2(eta2, xi2)
        ps.theta = phi_tp - phi_img

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

OBSERVATION_TARGETS = {
    "crab":            (83.63308, 22.01450),
    "mrk421":          (166.11379, 38.20883),
    "mrk501":          (253.46758, 39.76017),
    "galactic_center": (266.41683, -29.00781),
    "hessj1837-069":   (279.41500, -6.95000),
    "mgroj2019+37":    (303.15708, 36.18417),
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
