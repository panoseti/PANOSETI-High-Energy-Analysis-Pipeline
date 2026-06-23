import numpy as np
import scipy.interpolate

class RayBundle:
    def __init__(self, pos, dir, t=None, valid=None, energy_eV=None):
        """
        Initialize a bundle of rays.
        
        Args:
            pos: numpy array of shape (N, 3) representing the positions (x, y, z)
            dir: numpy array of shape (N, 3) representing the directions (dx, dy, dz)
            t: numpy array of shape (N,) representing the time of each ray
            valid: numpy array of shape (N,) boolean valid flags
            energy_eV: numpy array of shape (N,) representing photon energy in eV
        """
        self.pos = np.array(pos, dtype=float)
        self.dir = np.array(dir, dtype=float)
        
        # Normalize directions
        norms = np.linalg.norm(self.dir, axis=1, keepdims=True)
        # Avoid division by zero
        norms[norms == 0] = 1.0
        self.dir = self.dir / norms
        
        N = self.pos.shape[0]
        
        if t is None:
            self.t = np.zeros(N, dtype=float)
        else:
            self.t = np.array(t, dtype=float)
            
        if valid is None:
            self.valid = np.ones(N, dtype=bool)
        else:
            self.valid = np.array(valid, dtype=bool)
            
        if energy_eV is None:
            self.energy_eV = np.ones(N, dtype=float)
        else:
            if np.isscalar(energy_eV):
                self.energy_eV = np.full(N, energy_eV, dtype=float)
            else:
                self.energy_eV = np.array(energy_eV, dtype=float)

    @classmethod
    def create_parallel_beam(cls, num_rays, direction, diameter, distance, energy_eV=None):
        """
        Create a bundle of parallel rays originating from a circular disk.
        
        Args:
            num_rays: Number of rays to generate.
            direction: The common direction vector for all rays (list or array of 3 elements).
            diameter: The diameter of the circular disk where rays originate.
            distance: Distance from the origin to the center of the disk along the direction vector.
                      (Usually negative to have rays approach the origin).
            energy_eV: Photon energy in eV (scalar, array of length num_rays, or a callable/generator function).
        """
        direction = np.array(direction, dtype=float)
        norm = np.linalg.norm(direction)
        if norm == 0:
            raise ValueError("Direction vector cannot be zero.")
        d_hat = direction / norm
        
        # Center of the disk
        center = distance * d_hat
        
        # Create an orthonormal basis (u, v) for the plane orthogonal to d_hat
        if abs(d_hat[0]) < 0.9:
            arbitrary = np.array([1.0, 0.0, 0.0])
        else:
            arbitrary = np.array([0.0, 1.0, 0.0])
            
        u = np.cross(d_hat, arbitrary)
        u = u / np.linalg.norm(u)
        v = np.cross(d_hat, u)
        
        # Generate random points on a 2D disk (uniform areal distribution)
        theta = np.random.uniform(0, 2 * np.pi, num_rays)
        r = (diameter / 2.0) * np.sqrt(np.random.uniform(0, 1.0, num_rays))
        
        # Calculate positions
        pos = center + (r * np.cos(theta))[:, np.newaxis] * u + (r * np.sin(theta))[:, np.newaxis] * v
        
        # All rays have the same direction
        dirs = np.tile(d_hat, (num_rays, 1))
        
        # Handle generator function for energy
        if callable(energy_eV):
            try:
                # Try passing the number of rays (common for numpy random generators)
                energy_val = energy_eV(num_rays)
            except TypeError:
                # Fallback to calling it without arguments N times
                energy_val = np.array([energy_eV() for _ in range(num_rays)])
            except Exception:
                # If there's an issue with arguments, still try 0-arg fallback
                energy_val = np.array([energy_eV() for _ in range(num_rays)])
        else:
            energy_val = energy_eV
            
        return cls(pos, dirs, energy_eV=energy_val)

    def propagate_to_plane_y(self, d, speed=1.0):
        """
        Propagate rays forward in time to a plane defined by y + d = 0 (y = -d).
        
        Args:
            d: The d parameter for the plane y = -d
            speed: The speed of the rays (default 1.0, e.g., c in vacuum)
        """

        # We want to solve for distance s: pos_y + s * dir_y = -d
        # s = (-d - pos_y) / dir_y
        
        with np.errstate(divide='ignore', invalid='ignore'):
            s = (-d - self.pos[:, 1]) / self.dir[:, 1]
            
        # Rays moving parallel to the plane or in the wrong direction won't intersect in the forward direction
        invalid_mask = (self.dir[:, 1] == 0) | (np.isnan(s)) | (s < 0)
        self.valid[invalid_mask] = False
        
        # For valid rays, update position and time
        valid_idx = self.valid
        
        # Update pos: p = p0 + s * dir
        self.pos[valid_idx] += s[valid_idx, np.newaxis] * self.dir[valid_idx]
        
        # Update time: t = t0 + s / speed
        self.t[valid_idx] += s[valid_idx] / speed

    def propagate_to_y_plane(self, y_target, speed=1.0):
        """
        Propagate rays forward in time to a given y-coordinate.
        
        Args:
            y_target: The y-coordinate to propagate to.
            speed: The speed of the rays (default 1.0)
        """
        self.propagate_to_plane_y(-y_target, speed)

    def _refract(self, normal, n1, n2):
        """
        Helper function to perform vector refraction.
        normal: shape (N, 3), pointing towards the medium 1 (opposite to incident ray)
        """
        N = self.pos.shape[0]
        if np.isscalar(n1):
            n1 = np.full(N, n1, dtype=float)
        else:
            n1 = np.asarray(n1, dtype=float)
            if n1.size == 1:
                n1 = np.full(N, n1.item(), dtype=float)
                
        if np.isscalar(n2):
            n2 = np.full(N, n2, dtype=float)
        else:
            n2 = np.asarray(n2, dtype=float)
            if n2.size == 1:
                n2 = np.full(N, n2.item(), dtype=float)
                
        n_ratio = n1 / n2
        
        # dot product of dir and normal
        # normal should point against the incident ray for the standard formula
        # cos_theta_i = - (dir . normal)
        cos_i = -np.sum(self.dir * normal, axis=1)
        
        # If ray hits the back of the normal, flip the normal
        flip_mask = cos_i < 0
        normal[flip_mask] = -normal[flip_mask]
        cos_i[flip_mask] = -cos_i[flip_mask]
        
        sin2_t = (n_ratio ** 2) * (1.0 - cos_i ** 2)
        
        # Total internal reflection mask
        tir_mask = sin2_t > 1.0
        self.valid[tir_mask] = False
        
        valid_idx = self.valid
        
        cos_t = np.sqrt(1.0 - sin2_t[valid_idx])
        
        # Snell's law in vector form:
        # dir_out = n_ratio * dir_in + (n_ratio * cos_i - cos_t) * normal
        n_ratio_v = n_ratio[valid_idx]
        term1 = n_ratio_v[:, np.newaxis] * self.dir[valid_idx]
        term2 = (n_ratio_v * cos_i[valid_idx] - cos_t)[:, np.newaxis] * normal[valid_idx]
        
        self.dir[valid_idx] = term1 + term2
        
        # Re-normalize just in case
        norms = np.linalg.norm(self.dir[valid_idx], axis=1, keepdims=True)
        self.dir[valid_idx] = self.dir[valid_idx] / norms

    def refract_into_lens_y(self, n1, n2):
        """
        Refract rays at their current location into a lens with refractive index n2.
        The normal is aligned with the y-axis.
        """

        N = self.pos.shape[0]
        # Assuming the normal is [0, -1, 0] to oppose a ray traveling in +y direction,
        # but _refract handles normal flipping automatically.
        normal = np.zeros((N, 3))
        normal[:, 1] = 1.0
        
        self._refract(normal, n1, n2)

    def refract_out_of_lens_poly(self, n1, n2, poly_coeffs):
        """
        Refract rays out of the lens at their current position defined by a polynomial
        with a given sag vs rho^2.
        
        The zeroth element of poly_coeffs is ignored as the ray is assumed to be at 
        the outgoing surface.
        
        Args:
            n1: refractive index of current medium
            n2: refractive index of next medium
            poly_coeffs: list or array of coefficients [a0, a1, a2, ...] 
                         for the polynomial a0 + a1*rho^2 + a2*rho^4 + ...
                         The polynomial is assumed to be y = P(rho^2)
        """
        
        # rho^2 = x^2 + z^2
        x = self.pos[:, 0]
        z = self.pos[:, 2]
        rho2 = x**2 + z**2
        
        # The surface is F(x, y, z) = y - P(rho^2) = 0
        # Gradient of F:
        # dF/dx = -P'(rho^2) * 2x
        # dF/dy = 1
        # dF/dz = -P'(rho^2) * 2z
        
        # Calculate P'(rho^2)
        # P(u) = a0 + a1*u + a2*u^2 + ...
        # P'(u) = a1 + 2*a2*u + 3*a3*u^2 + ...
        
        # Ignore zeroth element by slicing from index 1
        # Polyval wants coefficients in highest-to-lowest order, 
        # so we handle evaluation manually for clarity since we have an array of rho2
        
        dp_du = np.zeros_like(rho2)
        
        for i in range(1, len(poly_coeffs)):
            power = i - 1
            coef = poly_coeffs[i]
            dp_du += i * coef * (rho2 ** power)
            
        df_dx = -dp_du * 2 * x
        df_dy = np.ones_like(rho2)
        df_dz = -dp_du * 2 * z
        
        N = self.pos.shape[0]
        normal = np.zeros((N, 3))
        normal[:, 0] = df_dx
        normal[:, 1] = df_dy
        normal[:, 2] = df_dz
        
        # Normalize the normal vectors
        norms = np.linalg.norm(normal, axis=1, keepdims=True)
        # Avoid div by zero (shouldn't happen since df_dy = 1)
        normal = normal / norms
        
        self._refract(normal, n1, n2)


def trace_telescope(ray_bundle, n_func, poly_coeffs, aperture_diameter, focal_plane_y):
    """
    Traces a bundle of rays through the telescope optics.
    
    Args:
        ray_bundle: RayBundle object containing the rays to trace.
        n_func: A callable that takes an array of photon energies (in eV) and 
                returns an array or scalar of refractive indices for the lens material.
        poly_coeffs: Coefficients for the polynomial describing the lens exit surface sag.
        aperture_diameter: Diameter of the entrance aperture.
        focal_plane_y: The y-coordinate of the focal plane.
        
    Returns:
        The updated ray_bundle.
    """
    # 1. Propagate to entrance aperture at y=0
    ray_bundle.propagate_to_y_plane(0.0)
    
    # 2. Entrance aperture cut
    r2 = ray_bundle.pos[:, 0]**2 + ray_bundle.pos[:, 2]**2
    aperture_mask = r2 > (aperture_diameter / 2.0)**2
    ray_bundle.valid[aperture_mask] = False
    
    # 3. Refract into lens at y=0 plane
    # Calculate refractive index for each ray's energy
    n_lens = n_func(ray_bundle.energy_eV)
    # Refract from air (n=1) into lens (n=n_lens)
    ray_bundle.refract_into_lens_y(1.0, n_lens)
    
    # 4. Refract out of lens at polynomial surface
    # Refract from lens (n=n_lens) into air (n=1)
    ray_bundle.refract_out_of_lens_poly(n_lens, 1.0, poly_coeffs)
    
    # 5. Propagate to focal plane
    ray_bundle.propagate_to_y_plane(focal_plane_y)
    
    return ray_bundle

def create_n_interpolator(datapack):
    """
    Returns an interpolator function for refractive index n(energy_eV) 
    using the given datapack.
    """
    ev = datapack["refractive_index"]["ev"]
    n = datapack["refractive_index"]["n"]
    
    # Ensure energies are sorted for interpolation
    sort_idx = np.argsort(ev)
    ev = np.array(ev)[sort_idx]
    n = np.array(n)[sort_idx]
    
    return scipy.interpolate.interp1d(ev, n, bounds_error=False, fill_value=(n[0], n[-1]))

def create_energy_generator(datapack, zn=0):
    """
    Returns a generator function that produces random photon energies (in eV)
    weighted by atmospheric transmission, lens transmission, and SiPM QE.
    
    Args:
        datapack: The optical model data pack dictionary.
        zn: Zenith angle in radians (default 0).
    """
    # Create common energy grid covering the typical range (e.g. 1 to 7 eV based on data)
    ev_grid = np.linspace(1.0, 7.0, 1000)
    
    # Atmospheric transmission
    atm_ev = datapack["atmospheric_absorption"]["ev"]
    atm_tau = datapack["atmospheric_absorption"]["tau"]
    sort_idx = np.argsort(atm_ev)
    atm_ev = np.array(atm_ev)[sort_idx]
    atm_tau = np.array(atm_tau)[sort_idx]
    tau_interp = np.interp(ev_grid, atm_ev, atm_tau, left=0, right=0)
    
    cos_zn = np.cos(zn)
    if cos_zn <= 0:
        cos_zn = 1e-6 # Prevent division by zero or negative
    t_atm = np.exp(-tau_interp / cos_zn)
    
    # PMMA transmission
    pmma_ev = datapack["pmma_transmission"]["ev"]
    pmma_eff = datapack["pmma_transmission"]["eff"]
    sort_idx = np.argsort(pmma_ev)
    pmma_ev = np.array(pmma_ev)[sort_idx]
    pmma_eff = np.array(pmma_eff)[sort_idx]
    t_pmma = np.interp(ev_grid, pmma_ev, pmma_eff, left=0, right=0)
    
    # SiPM QE
    sipm_ev = datapack["sipm_pde"]["ev"]
    sipm_eff = datapack["sipm_pde"]["eff"]
    sort_idx = np.argsort(sipm_ev)
    sipm_ev = np.array(sipm_ev)[sort_idx]
    sipm_eff = np.array(sipm_eff)[sort_idx]
    t_sipm = np.interp(ev_grid, sipm_ev, sipm_eff, left=0, right=0)
    
    # Total probability density function
    prob = t_atm * t_pmma * t_sipm
    
    # Normalize to create a CDF
    prob_sum = np.sum(prob)
    if prob_sum == 0:
        pdf = np.ones_like(prob) / len(prob)
    else:
        pdf = prob / prob_sum
    cdf = np.cumsum(pdf)
    
    def generator(num_rays=1):
        # Inverse transform sampling
        u = np.random.uniform(0, 1, num_rays)
        return np.interp(u, cdf, ev_grid)
        
    return generator

def calc_psf(direction, num_rays, datapack, zn=0):
    """
    Calculates the point spread function (PSF) by generating a bundle of rays
    and tracing them through the telescope.
    
    Args:
        direction: The direction of the incoming parallel rays (e.g. [0, -1, 0] for on-axis).
        num_rays: Number of rays to simulate.
        datapack: The loaded optical model data pack.
        zn: Zenith angle in radians for atmospheric absorption.
        
    Returns:
        The RayBundle object after propagating to the focal plane.
    """
    optics = datapack["thin_optical_model"]
    D = optics["D"]
    F = optics["F"]
    poly_coeffs = optics["p_out"]
    
    # Assuming rays travel roughly in -y direction towards the telescope at y=0.
    # To start them before the telescope, we place the generating disk behind the origin 
    # relative to the ray direction.
    distance = -10.0
    
    bundle = RayBundle.create_parallel_beam(
        num_rays=num_rays,
        direction=direction,
        diameter=D * 1.5, # Slightly larger than aperture to fully illuminate
        distance=distance,
        energy_eV=create_energy_generator(datapack, zn)
    )
    
    n_func = create_n_interpolator(datapack)
    
    # Focal plane is at y = -F (assuming lens is at y=0 and rays travel towards -y)
    traced_bundle = trace_telescope(
        ray_bundle=bundle,
        n_func=n_func,
        poly_coeffs=poly_coeffs,
        aperture_diameter=D,
        focal_plane_y=-F
    )
    
    return traced_bundle
