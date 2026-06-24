import json
import os
import numpy as np
import scipy.interpolate

_DEFAULT_DATAPACK_FILENAME = "optical_model_data_pack.json"

def load_datapack(filename=None):
    """
    Load the optical model data pack from a JSON file.

    Args:
        filename: Path to the JSON file.  If None (the default) the bundled
                  ``resources/optical_model_data_pack.json`` that ships with
                  this package is used.

    Returns:
        The data pack as a Python dictionary.
    """
    if filename is None:
        filename = os.path.join(
            os.path.dirname(__file__), "resources", _DEFAULT_DATAPACK_FILENAME
        )
    with open(filename, "r") as f:
        return json.load(f)


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

    def scatter_directions(self, sigma_theta):
        """
        Scatters the valid ray directions using a 2D Gaussian profile.
        
        Args:
            sigma_theta: The standard deviation of the scattering angle in radians.
        """
        if sigma_theta <= 0:
            return
            
        valid_idx = self.valid
        num_valid = np.sum(valid_idx)
        if num_valid == 0:
            return
            
        # Sample polar angle theta from Rayleigh distribution (equivalent to 2D Gaussian)
        # and azimuthal angle phi uniformly
        u1 = np.random.uniform(0, 1, num_valid)
        u2 = np.random.uniform(0, 1, num_valid)
        
        theta = sigma_theta * np.sqrt(-2.0 * np.log(u1 + 1e-12))
        phi = 2.0 * np.pi * u2
        
        cos_theta = np.cos(theta)[:, np.newaxis]
        sin_theta = np.sin(theta)[:, np.newaxis]
        cos_phi = np.cos(phi)[:, np.newaxis]
        sin_phi = np.sin(phi)[:, np.newaxis]
        
        V = self.dir[valid_idx]
        
        # Create orthogonal basis (U, W) for each V
        A1 = np.array([1.0, 0.0, 0.0])
        A2 = np.array([0.0, 1.0, 0.0])
        
        U1 = np.cross(V, A1)
        norm1 = np.linalg.norm(U1, axis=1)
        U2 = np.cross(V, A2)
        
        mask = norm1 > 0.5
        U = np.where(mask[:, np.newaxis], U1, U2)
        U = U / np.linalg.norm(U, axis=1, keepdims=True)
        
        W = np.cross(V, U)
        
        # Rotate V by theta and phi
        V_new = V * cos_theta + (U * cos_phi + W * sin_phi) * sin_theta
        V_new = V_new / np.linalg.norm(V_new, axis=1, keepdims=True)
        
        self.dir[valid_idx] = V_new

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


def trace_telescope(ray_bundle, n_func, poly_coeffs, aperture_diameter, focal_plane_y, scattering_sigma_theta=0.0):
    """
    Traces a bundle of rays through the telescope optics.
    
    Args:
        ray_bundle: RayBundle object containing the rays to trace.
        n_func: A callable that takes an array of photon energies (in eV) and 
                returns an array or scalar of refractive indices for the lens material.
        poly_coeffs: Coefficients for the polynomial describing the lens exit surface sag.
        aperture_diameter: Diameter of the entrance aperture.
        focal_plane_y: The y-coordinate of the focal plane.
        scattering_sigma_theta: Roughness scattering angle std dev in radians.
        
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
    
    # Add surface scattering due to roughness
    if scattering_sigma_theta > 0:
        ray_bundle.scatter_directions(scattering_sigma_theta)
    
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

def trace_parallel_ray_bundle(direction, num_rays, datapack, zn=0, focal_offset=0.0, energy_eV=None):
    """
    Generate a bundle of rays and trace them through the telescope.
    
    Args:
        direction: The direction of the incoming parallel rays (e.g. [0, -1, 0] for on-axis).
        num_rays: Number of rays to simulate.
        datapack: The loaded optical model data pack.
        zn: Zenith angle in radians for atmospheric absorption.
        focal_offset: Offset added to the nominal focal distance.
        
    Returns:
        The RayBundle object after propagating to the focal plane.
    """
    optics = datapack["thin_optical_model"]
    D = optics["D"]
    F = optics["F"]
    poly_coeffs = optics["p_out"]
    roughness = optics.get("roughness", 0.0)
    sigma_theta = roughness / F if F > 0 else 0.0
    
    # Assuming rays travel roughly in -y direction towards the telescope at y=0.
    # To start them before the telescope, we place the generating disk behind the origin 
    # relative to the ray direction.
    distance = -D
    
    energy_generator = None
    if energy_eV is None:
        energy_generator = create_energy_generator(datapack, zn)
    else:
        def generator(num_rays=1):
            return np.full(num_rays, energy_eV)
        energy_generator = generator

    bundle = RayBundle.create_parallel_beam(
        num_rays=num_rays,
        direction=direction,
        diameter=D * 1.01, # Slightly larger than aperture to fully illuminate
        distance=distance,
        energy_eV=energy_generator
    )
    
    n_func = create_n_interpolator(datapack)
    
    # Focal plane is at y = -(F + focal_offset)
    traced_bundle = trace_telescope(
        ray_bundle=bundle,
        n_func=n_func,
        poly_coeffs=poly_coeffs,
        aperture_diameter=D,
        focal_plane_y=-(F + focal_offset),
        scattering_sigma_theta=sigma_theta
    )
    
    return traced_bundle


def generate_psf_image(x, y, num_rays, datapack, npixel=None, pixel_spacing=None, zn=0, focal_offset=0.0, energy_eV=None):
    """
    Generate a PSF image by tracing a parallel ray bundle through the telescope
    and histogramming the ray positions on the focal plane.

    The ray bundle direction is chosen so that the prime (on-axis) ray of the
    bundle arrives at position (x_phys, z_phys) = (x * pixel_spacing, y * pixel_spacing)
    on the focal plane, where x and y are pixel-grid coordinates measured from
    the centre of the array (fractional pixel values are allowed).

    In the telescope coordinate frame the focal plane is the XZ plane (at the
    appropriate Y coordinate).  The user-facing ``y`` argument therefore maps
    to the Z axis of the telescope frame.

    Args:
        x: Target position along the X axis, in pixels from the centre of the
           image (can be fractional).
        y: Target position along the Z axis (telescope frame), in pixels from
           the centre of the image (can be fractional).
        num_rays: Number of rays to trace.
        datapack: The loaded optical model data pack.
        npixel: Number of pixels on each side of the square output image.  If
              None, ``datapack["thin_optical_model"]["npixel"]`` is used.
        pixel_spacing: Physical size of one pixel (same units as the datapack
                  geometry, typically metres).  If None,
                  ``datapack["thin_optical_model"]["pixel_spacing"]`` is used.
        zn: Zenith angle in radians for atmospheric absorption (default 0).
        focal_offset: Offset added to the nominal focal distance (default 0).
        energy_eV: If given, a fixed photon energy (eV) for all rays.  If
                   None the energy is drawn from the combined atmospheric /
                   lens / SiPM efficiency spectrum.

    Returns:
        image: numpy array of shape (npixel, npixel) containing the number of rays
               that landed in each pixel.  The first axis corresponds to X and
               the second axis to Z (telescope frame) / Y (user frame).
    """
    optics = datapack["thin_optical_model"]
    F = optics["F"]
    focal_length = F + focal_offset

    if npixel is None:
        npixel = optics["npixel"]
    if pixel_spacing is None:
        pixel_spacing = optics["pixel_spacing"]

    # Physical target position on the focal plane (metres)
    x_target = x * pixel_spacing
    z_target = y * pixel_spacing   # user's y → telescope Z

    # Direction of the incoming parallel beam.
    # A ray travelling along direction (dx, dy, dz) with dy < 0 arrives at the
    # focal plane (y = -focal_length) at approximately:
    #   x_fp ≈ (dx / |dy|) * focal_length
    #   z_fp ≈ (dz / |dy|) * focal_length
    # So to target (x_target, z_target) we need:
    #   dx/|dy| = x_target / focal_length
    #   dz/|dy| = z_target / focal_length
    # Choosing |dy| = 1 (before normalisation):
    dx = x_target / focal_length
    dy = -1.0
    dz = z_target / focal_length
    direction = np.array([dx, dy, dz])

    # Trace the bundle
    bundle = trace_parallel_ray_bundle(
        direction=direction,
        num_rays=num_rays,
        datapack=datapack,
        zn=zn,
        focal_offset=focal_offset,
        energy_eV=energy_eV,
    )

    # Select only valid rays
    valid = bundle.valid
    if not np.any(valid):
        return np.zeros((npixel, npixel), dtype=int)

    # Focal-plane positions of valid rays (X and Z in telescope frame)
    x_fp = bundle.pos[valid, 0]
    z_fp = bundle.pos[valid, 2]

    # Convert physical positions to pixel indices.
    # Pixel index 0 corresponds to the range [-npixel/2 * pixel_spacing, (-npixel/2 + 1) * pixel_spacing).
    # Centre of the image is between pixels npixel//2 - 1 and npixel//2 for even npixel,
    # or at pixel npixel//2 for odd npixel.
    half = npixel / 2.0
    ix = np.floor(x_fp / pixel_spacing + half).astype(int)
    iz = np.floor(z_fp / pixel_spacing + half).astype(int)

    # Keep only rays that fall within the image boundaries
    in_bounds = (ix >= 0) & (ix < npixel) & (iz >= 0) & (iz < npixel)
    ix = ix[in_bounds]
    iz = iz[in_bounds]

    # Histogram into the image array
    image = np.zeros((npixel, npixel), dtype=int)
    np.add.at(image, (ix, iz), 1)

    return image
