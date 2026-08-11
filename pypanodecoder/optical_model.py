import json
import os
import numpy as np
import scipy.interpolate
import scipy.optimize

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

    def _select_values(self, values, key):
        """Return a subset of a ray attribute using NumPy-style indexing."""
        selected = np.asarray(values[key])
        if values.ndim == 2 and selected.ndim == 1:
            selected = selected.reshape(1, -1)
        elif values.ndim == 1 and selected.ndim == 0:
            selected = selected.reshape(1)
        return selected

    def __getitem__(self, key):
        """Return a new RayBundle containing the rays selected by key."""
        return RayBundle(
            self._select_values(self.pos, key),
            self._select_values(self.dir, key),
            t=self._select_values(self.t, key),
            valid=self._select_values(self.valid, key),
            energy_eV=self._select_values(self.energy_eV, key),
        )

    def __setitem__(self, key, value):
        """Assign another RayBundle or a subset of one into the selected rays."""
        if not isinstance(value, RayBundle):
            raise TypeError("RayBundle assignment requires another RayBundle")

        target_pos = self.pos[key]
        source_pos = value.pos
        if np.shape(target_pos) != np.shape(source_pos):
            raise ValueError("Assigned RayBundle has incompatible shape")

        self.pos[key] = source_pos
        self.dir[key] = value.dir
        self.t[key] = value.t
        self.valid[key] = value.valid
        self.energy_eV[key] = value.energy_eV

    def copy(self):
        """Return a copy of the ray bundle with independent arrays."""
        return RayBundle(
            self.pos.copy(),
            self.dir.copy(),
            t=self.t.copy(),
            valid=self.valid.copy(),
            energy_eV=self.energy_eV.copy(),
        )

    def _ray_parameter(self, value, name):
        """
        Return a scalar or per-ray parameter as an array of shape (N,).
        """
        N = self.pos.shape[0]
        value = np.asarray(value, dtype=float)
        if value.ndim == 0 or value.size == 1:
            return np.full(N, value.item(), dtype=float)
        if value.size != N:
            raise ValueError(f"{name} must be a scalar or have one value per ray")
        return value.reshape(N)

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

    @classmethod
    def create_grid_parallel_beam(cls, num_rays, direction, diameter, distance, energy_eV=None):
        """
        Create a bundle of parallel rays originating from a regular square grid on a circular disk.
        
        Args:
            num_rays: Approximate number of rays to generate. The actual number may slightly differ due to grid mapping.
            direction: The common direction vector for all rays (list or array of 3 elements).
            diameter: The diameter of the circular disk where rays originate.
            distance: Distance from the origin to the center of the disk along the direction vector.
                      (Usually negative to have rays approach the origin).
            energy_eV: Photon energy in eV (scalar, array of length actual_num_rays, or a callable/generator function).
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
        
        # Calculate grid spacing
        radius = diameter / 2.0
        area = np.pi * radius**2
        
        if num_rays <= 0:
            raise ValueError("num_rays must be greater than 0.")
            
        spacing = np.sqrt(area / num_rays)
        
        # Generate grid points (with one ray at 0,0)
        max_idx = int(np.floor(radius / spacing))
        indices = np.arange(-max_idx, max_idx + 1)
        i_grid, j_grid = np.meshgrid(indices, indices)
        
        # Filter points within the circle
        mask = (i_grid * spacing)**2 + (j_grid * spacing)**2 <= radius**2
        i_valid = i_grid[mask]
        j_valid = j_grid[mask]
        
        actual_num_rays = len(i_valid)
        
        # Calculate positions
        pos = center + (i_valid * spacing)[:, np.newaxis] * u + (j_valid * spacing)[:, np.newaxis] * v
        
        # All rays have the same direction
        dirs = np.tile(d_hat, (actual_num_rays, 1))
        
        # Handle generator function for energy
        if callable(energy_eV):
            try:
                # Try passing the actual number of rays (common for numpy random generators)
                energy_val = energy_eV(actual_num_rays)
            except TypeError:
                # Fallback to calling it without arguments N times
                energy_val = np.array([energy_eV() for _ in range(actual_num_rays)])
            except Exception:
                # If there's an issue with arguments, still try 0-arg fallback
                energy_val = np.array([energy_eV() for _ in range(actual_num_rays)])
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

    def propagate_to_cone_y(self, m, d, speed=1.0):
        """
        Propagate rays forward in time to a cone aligned with the y-axis.

        The cone surface is defined by

            y + d = m * sqrt(x^2 + z^2)

        so m=0 is the plane y + d = 0, matching propagate_to_plane_y(d).

        Args:
            m: The cone slope in y per radial distance. Can be a scalar or
               one value per ray.
            d: The cone vertex offset; the vertex is at (0, -d, 0). Can be a
               scalar or one value per ray.
            speed: The speed of the rays (default 1.0, e.g., c in vacuum)
        """
        m = self._ray_parameter(m, "m")
        d = self._ray_parameter(d, "d")

        if np.all(m == 0):
            self.propagate_to_plane_y(d, speed)
            return

        N = self.pos.shape[0]
        x0 = self.pos[:, 0]
        y0 = self.pos[:, 1] + d
        z0 = self.pos[:, 2]
        vx = self.dir[:, 0]
        vy = self.dir[:, 1]
        vz = self.dir[:, 2]
        m2 = m * m

        # Solve (y0 + s*vy)^2 = m^2*((x0 + s*vx)^2 + (z0 + s*vz)^2).
        a = vy**2 - m2 * (vx**2 + vz**2)
        b = 2.0 * (y0 * vy - m2 * (x0 * vx + z0 * vz))
        c = y0**2 - m2 * (x0**2 + z0**2)

        eps = 1e-14
        tol = 1e-10
        candidates = np.full((N, 3), np.inf)
        active = self.valid

        plane_idx = active & (m == 0.0)
        if np.any(plane_idx):
            with np.errstate(divide='ignore', invalid='ignore'):
                s_plane = -y0[plane_idx] / vy[plane_idx]
            valid_plane = (vy[plane_idx] != 0) & np.isfinite(s_plane) & (s_plane >= 0.0)
            plane_rows = np.flatnonzero(plane_idx)
            candidates[plane_rows[valid_plane], 0] = s_plane[valid_plane]

        disc = b**2 - 4.0 * a * c
        cone_idx = active & (m != 0.0)
        quad_idx = cone_idx & (np.abs(a) > eps) & (disc >= -tol)
        if np.any(quad_idx):
            sqrt_disc = np.sqrt(np.maximum(disc[quad_idx], 0.0))
            candidates[quad_idx, 0] = (-b[quad_idx] - sqrt_disc) / (2.0 * a[quad_idx])
            candidates[quad_idx, 1] = (-b[quad_idx] + sqrt_disc) / (2.0 * a[quad_idx])

        linear_idx = cone_idx & (np.abs(a) <= eps) & (np.abs(b) > eps)
        if np.any(linear_idx):
            candidates[linear_idx, 0] = -c[linear_idx] / b[linear_idx]

        on_surface_idx = (
            cone_idx
            & (np.abs(a) <= eps)
            & (np.abs(b) <= eps)
            & (np.abs(c) <= tol)
            & (np.abs(y0 - m * np.sqrt(x0**2 + z0**2)) <= tol)
        )
        candidates[on_surface_idx, 0] = 0.0

        s = candidates
        finite_candidates = np.isfinite(s)
        x = np.full_like(s, np.inf)
        y = np.full_like(s, np.inf)
        z = np.full_like(s, np.inf)
        finite_rows, finite_cols = np.where(finite_candidates)
        finite_s = s[finite_rows, finite_cols]
        x[finite_rows, finite_cols] = x0[finite_rows] + finite_s * vx[finite_rows]
        y[finite_rows, finite_cols] = y0[finite_rows] + finite_s * vy[finite_rows]
        z[finite_rows, finite_cols] = z0[finite_rows] + finite_s * vz[finite_rows]
        rho = np.full_like(s, np.inf)
        residual = np.full_like(s, np.inf)
        residual_tol = np.full_like(s, np.inf)
        rho[finite_rows, finite_cols] = np.sqrt(
            x[finite_rows, finite_cols]**2 + z[finite_rows, finite_cols]**2
        )
        residual[finite_rows, finite_cols] = (
            y[finite_rows, finite_cols]
            - m[finite_rows] * rho[finite_rows, finite_cols]
        )
        residual_tol[finite_rows, finite_cols] = tol * (
            1.0
            + np.abs(y[finite_rows, finite_cols])
            + np.abs(m[finite_rows] * rho[finite_rows, finite_cols])
        )

        valid_candidates = finite_candidates & (s >= 0.0) & (np.abs(residual) <= residual_tol)
        candidates = np.where(valid_candidates, candidates, np.inf)
        s_hit = np.min(candidates, axis=1)

        hit_idx = active & np.isfinite(s_hit)
        self.valid[active & ~hit_idx] = False

        self.pos[hit_idx] += s_hit[hit_idx, np.newaxis] * self.dir[hit_idx]
        self.t[hit_idx] += s_hit[hit_idx] / speed

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

    def _normalize_normal(self, normal):
        """Normalize a surface-normal array of shape (N, 3)."""
        normal = np.asarray(normal, dtype=float)
        if normal.ndim == 1:
            normal = normal.reshape(1, -1)

        norms = np.linalg.norm(normal, axis=1, keepdims=True)
        norms[norms == 0] = 1.0
        return normal / norms

    def _normal_poly(self, poly_coeffs):
        """Return the outward unit normal for a polynomial surface y = P(rho^2)."""
        poly_coeffs = np.asarray(poly_coeffs, dtype=float)

        x = self.pos[:, 0]
        z = self.pos[:, 2]
        rho2 = x**2 + z**2

        dp_du = np.zeros_like(rho2)
        for i in range(1, len(poly_coeffs)):
            power = i - 1
            coef = poly_coeffs[i]
            dp_du += i * coef * (rho2**power)

        normal = np.zeros((self.pos.shape[0], 3), dtype=float)
        normal[:, 0] = -2.0 * dp_du * x
        normal[:, 1] = 1.0
        normal[:, 2] = -2.0 * dp_du * z
        return self._normalize_normal(normal)

    def _normal_cone(self, m):
        """Return the outward unit normal for a cone y + d = m * sqrt(x^2 + z^2)."""
        m = self._ray_parameter(m, "m")

        x = self.pos[:, 0]
        z = self.pos[:, 2]
        rho = np.sqrt(x**2 + z**2)

        normal = np.zeros((self.pos.shape[0], 3), dtype=float)
        normal[:, 1] = 1.0

        cone_idx = m != 0.0
        non_vertex_idx = cone_idx & (rho > 0.0)
        normal[non_vertex_idx, 0] = -m[non_vertex_idx] * x[non_vertex_idx] / rho[non_vertex_idx]
        normal[non_vertex_idx, 2] = -m[non_vertex_idx] * z[non_vertex_idx] / rho[non_vertex_idx]
        return self._normalize_normal(normal)

    def _reflect(self, normal):
        """Reflect the ray directions across the given surface normal vectors."""
        normal = self._normalize_normal(normal)
        dot = np.sum(self.dir * normal, axis=1)
        self.dir = self.dir - 2.0 * dot[:, np.newaxis] * normal

    def reflect_in_y_plane(self):
        """Reflect rays in a plane parallel to the x-z plane."""
        normal = np.zeros((self.pos.shape[0], 3), dtype=float)
        normal[:, 1] = 1.0
        self._reflect(normal)

    def reflect_in_poly(self, poly_coeffs):
        """Reflect rays in a polynomial surface y = P(rho^2)."""
        self._reflect(self._normal_poly(poly_coeffs))

    def reflect_in_cone(self, m):
        """Reflect rays in a cone defined by y + d = m * sqrt(x^2 + z^2)."""
        self._reflect(self._normal_cone(m))

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

    def refract_out_of_lens_cone(self, n1, n2, m):
        """
        Refract rays out of the lens at their current position on a cone.

        The cone surface is assumed to be y + d = m * sqrt(x^2 + z^2).
        The offset d is not needed because it does not affect the surface
        normal.

        Args:
            n1: refractive index of current medium
            n2: refractive index of next medium
            m: The cone slope in y per radial distance. Can be a scalar or
               one value per ray.
        """
        normal = self._normal_cone(m)
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
        
        normal = self._normal_poly(poly_coeffs)
        self._refract(normal, n1, n2)


def trace_telescope_thin(ray_bundle, n_func, optics, focal_offset=0.0):
    """
    Traces a bundle of rays through the telescope optics with the thin lens approximation.
    
    Args:
        ray_bundle: RayBundle object containing the rays to trace.
        n_func: A callable that takes an array of photon energies (in eV) and 
                returns an array or scalar of refractive indices for the lens material.
        optics: An object containing the optical properties of the lens.
        focal_offset: The offset of the focal plane from the default position.
        
    Returns:
        The updated ray_bundle.
    """

    # 0. Extract optical parameters
    aperture_diameter       = optics["D"]
    focal_length            = optics["F"]
    poly_coeffs             = optics["p_out"]
    roughness               = optics.get("roughness", 0.0)
    reflection_probability  = max(optics.get("reflection_probability", 0.0), 0.0)
    focal_plane_y = -(focal_length + focal_offset)
    scattering_sigma_theta = roughness / focal_length if focal_length > 0 else 0.0

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
    
    # 4. Reflect some rays due to surface reflection if configured
    if reflection_probability > 0:
        eligable_rays = ray_bundle.valid.copy()
        while True:
            eligable_rays &= np.random.uniform(0.0, 1.0, len(eligable_rays)) < reflection_probability
            if not np.any(eligable_rays):
                break
            reflect_bundle = ray_bundle[eligable_rays]
            reflect_bundle.reflect_in_poly(poly_coeffs)
            reflect_bundle.reflect_in_y_plane()
            ray_bundle[eligable_rays] = reflect_bundle

    # 5. Refract out of lens at polynomial surface
    # Refract from lens (n=n_lens) into air (n=1)
    ray_bundle.refract_out_of_lens_poly(n_lens, 1.0, poly_coeffs)
    
    # 6. Add surface scattering due to roughness
    if scattering_sigma_theta > 0:
        ray_bundle.scatter_directions(scattering_sigma_theta)
    
    # 7. Propagate to focal plane
    ray_bundle.propagate_to_y_plane(focal_plane_y)
    
    return ray_bundle

def _groove_configuration(aperture_diameter, groove_width, draft_angle, poly_coeffs):
    """
    Compute the groove configuration for the telescope optics.
    """
    poly_derivative_coeffs = np.polyder(poly_coeffs)
    ngroove = int(np.ceil(aperture_diameter/2/groove_width))
    groove_rho_inner = np.arange(ngroove) * groove_width
    groove_rho_mid = groove_rho_inner + groove_width/2
    groove_m = np.polyval(poly_derivative_coeffs, groove_rho_mid**2) * 2 * groove_rho_mid
    groove_d = groove_m * groove_rho_inner
    groove_rho_outer = groove_rho_inner + groove_width/(1 + groove_m*np.tan(np.deg2rad(draft_angle)))
    return ngroove, groove_rho_inner, groove_rho_outer,  groove_m, groove_d

def trace_telescope_thick(ray_bundle, n_func, optics, focal_offset=0.0):
    """
    Traces a bundle of rays through the telescope optics with the thick lens approximation.
    
    Args:
        ray_bundle: RayBundle object containing the rays to trace.
        n_func: A callable that takes an array of photon energies (in eV) and 
                returns an array or scalar of refractive indices for the lens material.
        optics: An object containing the optical properties of the lens.
        focal_offset: The offset of the focal plane from the default position.
        
    Returns:
        The updated ray_bundle.
    """

    # 0. Extract optical parameters and compute grooves
    aperture_diameter   = optics["D"]
    focal_length        = optics["F"]
    roughness           = optics.get("roughness", 0.0)
    thickness           = optics["thickness"]
    groove_width        = optics["groove_width"]
    draft_angle         = max(optics.get("draft_angle", 0.0), 0.0)

    focal_plane_y = -(focal_length + focal_offset)
    poly_coeffs = np.flipud(np.array(optics["p_out"]))
    scattering_sigma_theta = roughness / focal_length

    ngroove, groove_rho_inner, groove_rho_outer,  groove_m, groove_d = _groove_configuration(aperture_diameter, groove_width, draft_angle, poly_coeffs)

    # 1. Propagate to entrance aperture at y=0
    ray_bundle.propagate_to_y_plane(thickness)
    
    # 2. Entrance aperture cut
    r2 = ray_bundle.pos[:, 0]**2 + ray_bundle.pos[:, 2]**2
    aperture_mask = r2 > (aperture_diameter / 2.0)**2
    ray_bundle.valid[aperture_mask] = False
    
    # 3. Refract into lens
    # Calculate refractive index for each ray's energy
    n_lens = n_func(ray_bundle.energy_eV)
    # Refract from air (n=1) into lens (n=n_lens)
    ray_bundle.refract_into_lens_y(1.0, n_lens)
    
    # 4. Propagate through the lens thickness to find groove
    test_ray_bundle = ray_bundle.copy()
    test_ray_bundle.propagate_to_y_plane(0.0)
    rhoexit = np.sqrt(test_ray_bundle.pos[:, 0]**2 + test_ray_bundle.pos[:, 2]**2)
    igroove = np.minimum(np.floor(rhoexit / groove_width).astype(int), ngroove - 1)

    # 5. Propagate to the nominal conical groove surface
    test_ray_bundle = ray_bundle.copy()
    test_ray_bundle.propagate_to_cone_y(groove_m[igroove], groove_d[igroove])
    rhogroove = np.sqrt(test_ray_bundle.pos[:, 0]**2 + test_ray_bundle.pos[:, 2]**2)
    test_ray_bundle.valid &= (rhogroove >= groove_rho_inner[igroove]) & (rhogroove <= groove_rho_outer[igroove])

    # 6. Test propagation to the previous conical groove surface
    igroove2 = np.maximum(np.minimum(np.floor(rhoexit / groove_width).astype(int) - 1, ngroove - 1), 0)
    test_ray_bundle2 = ray_bundle.copy()
    test_ray_bundle2.propagate_to_cone_y(groove_m[igroove2], groove_d[igroove2])
    rhogroove2 = np.sqrt(test_ray_bundle2.pos[:, 0]**2 + test_ray_bundle2.pos[:, 2]**2)
    test_ray_bundle2.valid &= (rhogroove2 >= groove_rho_inner[igroove2]) & (rhogroove2 <= groove_rho_outer[igroove2])

    # 7. Merge the two groove tests
    ray_bundle = test_ray_bundle
    ray_bundle[test_ray_bundle2.valid] = test_ray_bundle2[test_ray_bundle2.valid]
    igroove[test_ray_bundle2.valid] = igroove2[test_ray_bundle2.valid]

    # 8. Refract out of lens at conical surface
    # Refract from lens (n=n_lens) into air (n=1)
    ray_bundle.refract_out_of_lens_cone(n_lens, 1.0, groove_m[igroove])
    
    # Add surface scattering due to roughness
    if scattering_sigma_theta > 0:
        ray_bundle.scatter_directions(scattering_sigma_theta)
    
    # 9. Propagate to focal plane
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
    Source photon flux spectrum is assumed to be flat (Cherenkov) - we 
    should also suuport simple thermal spectra with a given B-V.
    
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

def trace_parallel_ray_bundle(direction, num_rays, datapack, thick_lens=False, zn=0, focal_offset=0.0, energy_eV=None, use_grid=False, fixed_n=None):
    """
    Generate a bundle of rays and trace them through the telescope.
    
    Args:
        direction: The direction of the incoming parallel rays (e.g. [0, -1, 0] for on-axis).
        num_rays: Number of rays to simulate.
        datapack: The loaded optical model data pack.
        thick_lens: If True, use the thick lens model; otherwise, use the thin lens model.
        zn: Zenith angle in radians for atmospheric absorption.
        focal_offset: Offset added to the nominal focal distance.
        energy_eV: Photon energy in eV (default None, drawn from standard spectrum).
        use_grid: If True, generate rays on a regular square grid instead of uniformly randomly.
        fixed_n: Fixed refractive index to use for all rays, overriding the datapack interpolator.
        
    Returns:
        The RayBundle object after propagating to the focal plane.
    """
    optics = datapack["optical_model"]
    D = optics["D"]
    
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

    if use_grid:
        bundle = RayBundle.create_grid_parallel_beam(
            num_rays=num_rays,
            direction=direction,
            diameter=D * 1.01, # Slightly larger than aperture to fully illuminate
            distance=distance,
            energy_eV=energy_generator
        )
    else:
        bundle = RayBundle.create_parallel_beam(
            num_rays=num_rays,
            direction=direction,
            diameter=D * 1.01, # Slightly larger than aperture to fully illuminate
            distance=distance,
            energy_eV=energy_generator
        )
    
    if fixed_n is None:
        n_func = create_n_interpolator(datapack)
    else:
        def n_func(energy):
            if np.isscalar(energy):
                return fixed_n
            return np.full_like(energy, fixed_n, dtype=float)
    
    # Focal plane is at y = -(F + focal_offset)
    if thick_lens:
        traced_bundle = trace_telescope_thick(
            ray_bundle=bundle,
            n_func=n_func,
            optics=optics,
            focal_offset=focal_offset
        )
    else:
        traced_bundle = trace_telescope_thin(
            ray_bundle=bundle,
            n_func=n_func,
            optics=optics,
            focal_offset=focal_offset
        )
    
    return traced_bundle


def generate_psf_image(x, y, num_rays, datapack, thick_lens=False, npixel=None, pixel_spacing=None, zn=0, focal_offset=0.0, energy_eV=None, calc_diameter=False, diameter_quantile=0.80):
    """
    Generate a PSF image by tracing a parallel ray bundle through the telescope
    and histogramming the ray positions on the focal plane.

    The ray bundle direction is chosen so that the principal ray of the
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
        thick_lens: If True, use the thick lens model; otherwise, use the thin lens model.
        npixel: Number of pixels on each side of the square output image.  If
              None, ``datapack["optical_model"]["npixel"]`` is used.
        pixel_spacing: Physical size of one pixel (same units as the datapack
                  geometry, typically cm).  If None,
                  ``datapack["optical_model"]["pixel_spacing"]`` is used.
        zn: Zenith angle in radians for atmospheric absorption (default 0).
        focal_offset: Offset added to the nominal focal distance (default 0).
        energy_eV: If given, a fixed photon energy (eV) for all rays.  If
                   None the energy is drawn from the combined atmospheric /
                   lens / SiPM efficiency spectrum.
        calc_diameter: If True, also compute and return the centre and diameter
                       of the smallest circle enclosing ``diameter_quantile``
                       of the valid rays.  Default is False.
        diameter_quantile: Fraction of rays to enclose when computing the
                           diameter (e.g. 0.80 for d80, 0.50 for d50).
                           Only used when ``calc_diameter`` is True.
                           Default is 0.80.

    Returns:
        image: numpy array of shape (npixel, npixel) containing the number of
               rays that landed in each pixel.  The first axis corresponds to X
               and the second axis to Z (telescope frame) / Y (user frame).
        If calc_diameter is True, two additional values are returned:
        center: (cx, cz) tuple giving the centre of the enclosing circle in
                pixels, measured from the centre of the image (same convention
                as the x and y input arguments).
        diameter: diameter of the smallest circle enclosing
                  ``diameter_quantile`` of all valid rays on the focal plane
                  (including those outside the image bounds), in pixels.
    """
    optics = datapack["optical_model"]
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
        thick_lens=thick_lens,
        zn=zn,
        focal_offset=focal_offset,
        energy_eV=energy_eV,
    )

    # Select only valid rays
    valid = bundle.valid
    if not np.any(valid):
        image = np.zeros((npixel, npixel), dtype=int)
        if calc_diameter:
            return image, (np.nan, np.nan), np.nan
        return image

    # Focal-plane positions of valid rays (X and Z in telescope frame)
    x_fp = bundle.pos[valid, 0]
    z_fp = bundle.pos[valid, 2]
    nvalid = len(x_fp)

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
    np.add.at(image, (iz, ix), 1)

    # Return median and men positions in pixel units measured from the image
    # centre, for consistency with the x and y input arguments.
    x_pix = x_fp / pixel_spacing   # pixels from image centre
    z_pix = z_fp / pixel_spacing
    r_pix = np.sqrt(x_pix**2 + z_pix**2)

    center_median = ( np.median(x_pix), np.median(z_pix) )
    center_mean = ( np.mean(x_pix), np.mean(z_pix) )
    r_median = np.median(r_pix)
    r_mean = np.median(r_pix)

    if not calc_diameter:
        return nvalid, image, center_mean, r_mean, center_median, r_median

    # --- Smallest circle enclosing 100*Q% of all valid rays (dQ) ---
    # All valid rays on the focal plane are included, even those outside the
    # image bounds.  

    # For a candidate centre the radius needed to enclose 100*Q% of the rays
    # is the 100*Qth percentile of the distances.  Minimising over all candidate
    # centres gives the smallest such circle.
    def _radius_for_fraction(centre):
        dx = x_pix - centre[0]
        dz = z_pix - centre[1]
        return np.percentile(np.sqrt(dx**2 + dz**2), diameter_quantile * 100.0)

    x0 = np.array(center_median)
    result = scipy.optimize.minimize(
        _radius_for_fraction,
        x0,
        method="Nelder-Mead",
        options={"xatol": 1e-3, "fatol": 1e-3},  # tolerances in pixels
    )
    center_Q = (result.x[0], result.x[1])
    diameter_Q = 2.0 * result.fun

    return nvalid, image, center_mean, r_mean, center_median, r_median, center_Q, diameter_Q
