import numpy as np

from .transforms import transform_points, transform_dir, lab_frame_points, lab_frame_dir
from .constants import SURV_CONST, MAX_INTERSECTION_ITERATIONS, INTERSECTION_CONVERGENCE_TOLERANCE, MAX_REFRACTION_ITERATIONS, REFRACTION_CONVERGENCE_TOLERANCE
from .exceptions import NormalizationError, NotOnSurfaceError
from .geometry import geometry

from typing import Union, List

class RayGroup:
    """Class for a group of rays and their vectorized propagation through surfaces.

    Attributes
    ----------
    P : np.ndarray of shape (N,3)
        Positions of rays in the lab frame.
    D : np.ndarray of shape (N,3)
        Direction cosines for rays in the lab frame.
    N : np.ndarray of shape (N,)
        Index of refraction for each ray.
    wvl : np.ndarray of shape (N,)
        Wavelength of each ray in microns.
    active : np.ndarray of bool, shape (N,)
        Boolean array indicating if each ray is still active.
    P_hist : list of np.ndarray
        History of positions for each propagation step.
    D_hist : list of np.ndarray
        History of directions for each propagation step.
    normal : np.ndarray of shape (N,3)
        Surface normals at intersection points (updated per surface).
    """

    def __init__(self, P: Union[List[List[float]], np.ndarray],
                       D: Union[List[List[float]], np.ndarray],
                       wvl: Union[float, List[float]] = 0.55,
                       N_0: float = 1.0):
        self.P = np.array(P)  # shape (N,3)
        self.D = np.array(D)  # shape (N,3)
        N_rays = self.P.shape[0]
        # If wvl is a scalar, replicate for all rays.
        if np.isscalar(wvl):
            self.wvl = np.full(N_rays, wvl, dtype=np.float64)
        else:
            self.wvl = np.asarray(wvl, dtype=np.float64)  # type: ignore
        self.N = np.full(N_rays, N_0)
        self.active = np.ones(N_rays, dtype=bool)
        # Check that all direction vectors are normalized.
        norms = np.linalg.norm(self.D, axis=1)
        if not np.allclose(norms, 1.0, atol=0.01):
            raise NormalizationError()
        self.P_hist = [self.P.copy()]
        self.D_hist = [self.D.copy()]
        self.normal = np.zeros_like(self.P)

    def transform(self, surface: geometry) -> None:
        """Transforms positions and directions into the surface's reference frame."""
        self.P = transform_points(surface.R, surface, self.P)
        self.D = transform_dir(surface.R, surface, self.D)

    def find_intersections(self, surface: geometry) -> None:
        """Vectorized computation of intersections of rays with a surface."""
        active_idx = np.where(self.active)[0]
        if active_idx.size == 0:
            return
        # For active rays, compute initial guess for s (parameter along ray)
        P_active = self.P[active_idx]
        D_active = self.D[active_idx]
        # Initial guess: s0 = -P_z / D_z
        s0 = -P_active[:, 2] / D_active[:, 2]
        # Initialize s with s0 instead of zeros
        s = s0.copy()
        error = np.full(s0.shape, 1.0)
        n_iter = 0
        max_iter = int(MAX_INTERSECTION_ITERATIONS)
        tol = INTERSECTION_CONVERGENCE_TOLERANCE
        # Newton-Raphson iteration for all active rays
        while np.any(error > tol) and n_iter < max_iter:
            idx = error > tol
            if not np.any(idx):
                break
            # Update guess points for these rays
            P_temp = P_active[idx] + D_active[idx] * s[idx, None]
            try:
                func, deriv = surface.get_surface_vector(P_temp)
            except Exception:
                # If evaluation fails, mark these rays as inactive
                self.active[active_idx[idx]] = False
                error[idx] = 0.0
                continue
            # Compute derivative along ray direction: dot(deriv, D)
            deriv_dot = np.sum(deriv * D_active[idx], axis=1)
            # Newton update: s_new = s - func / deriv_dot
            delta_s = func / deriv_dot
            s[idx] = s[idx] - delta_s
            error[idx] = np.abs(func)
            n_iter += 1
        # Update positions for active rays based on final s
        P_new = P_active + D_active * s[:, None]
        self.P[active_idx] = P_new
        # Update surface normals for active rays
        try:
            _, norm = surface.get_surface_vector(P_new)
            self.normal[active_idx] = norm
        except Exception:
            self.active[active_idx] = False

    def interact(self, surface: geometry) -> None:
        """Updates ray directions based on interaction type (reflection or refraction).
        
        Uses vectorized operations for active rays.
        """
        active_idx = np.where(self.active)[0]
        if active_idx.size == 0:
            return
        # Determine refractive index ratio mu for active rays.
        if hasattr(surface, 'glass'):
            mu = self.N[active_idx] / surface.glass(self.wvl[active_idx])
        else:
            mu = self.N[active_idx] / surface.N
        # Compute constant a = mu * dot(D, normal) / ||normal||^2
        D_active = self.D[active_idx]
        norm_active = self.normal[active_idx]
        dot_dn = np.sum(D_active * norm_active, axis=1)
        norm_sq = np.sum(norm_active ** 2, axis=1)
        a = mu * dot_dn / norm_sq
        # Compute b = (mu**2 - 1) / norm_sq
        b = (mu**2 - 1) / norm_sq
        if surface.action == 'stop':
            # No change to direction.
            return
        elif surface.action == 'reflection':
            # Reflection: update D = D - 2*a * normal
            self.D[active_idx] = D_active - 2 * a[:, None] * norm_active
        elif surface.action == 'refraction':
            # Identify rays undergoing total internal reflection.
            tir = b > a**2
            if np.any(tir):
                idx_tir = active_idx[tir]
                self.D[idx_tir] = self.D[idx_tir] - 2 * a[tir, None] * norm_active[tir]
            # For rays that are meant to refract:
            ref_idx = active_idx[~tir]
            if ref_idx.size > 0:
                D_ref = self.D[ref_idx]
                norm_ref = self.normal[ref_idx]
                # a_ref and b_ref correspond to the rays that refract.
                a_ref = a[~tir]
                b_ref = b[~tir]
                mu_ref = mu[~tir]  # refractive index ratios for these rays
                # Initial guess for G: -b/(2*a)
                G = -b_ref / (2 * a_ref)
                error = np.full(G.shape, 1.0)
                n_iter = 0
                max_iter = int(MAX_REFRACTION_ITERATIONS)
                tol = REFRACTION_CONVERGENCE_TOLERANCE
                # Newton–Raphson iteration (updating only rays that haven't converged)
                while np.any(error > tol) and n_iter < max_iter:
                    mask = error > tol  # update only nonconverged rays
                    G[mask] = (G[mask]**2 - b_ref[mask]) / (2 * (G[mask] + a_ref[mask]))
                    error[mask] = np.abs(G[mask]**2 + 2 * a_ref[mask] * G[mask] + b_ref[mask])
                    n_iter += 1
                # Identify which rays did not converge
                not_converged = error > tol
                if np.any(not_converged):
                    self.active[ref_idx[not_converged]] = False
                # For rays that converged, update the ray directions and indices.
                converged = ~not_converged
                if np.any(converged):
                    self.D[ref_idx[converged]] = (
                        mu_ref[converged, None] * D_ref[converged]
                        + G[converged, None] * norm_ref[converged]
                    )
                    if hasattr(surface, 'glass'):
                        self.N[ref_idx[converged]] = surface.glass(self.wvl[ref_idx[converged]])
                    else:
                        self.N[ref_idx[converged]] = surface.N

    def ray_lab_frame(self, surface: geometry) -> None:
        """Transforms positions and directions back to the lab frame."""
        self.P = lab_frame_points(surface.R, surface, self.P)
        self.D = lab_frame_dir(surface.R, surface, self.D)

    def update_history(self) -> None:
        """Appends current positions and directions to history arrays."""
        self.P_hist.append(self.P.copy())
        self.D_hist.append(self.D.copy())

    def to_ray_list(self):
        dummy_list = []
        # self.P_hist is a list of arrays [P0, P1, P2, ...]
        # each entry is shape (N, 3) for N rays
        for i in range(self.P.shape[0]):
            dummy = type("DummyRay", (), {})()
            
            # Build the entire per‐surface history for ray i
            dummy.P_hist = [p[i] for p in self.P_hist]
            dummy.D_hist = [d[i] for d in self.D_hist]
            
            # Keep any other attributes needed
            dummy.active = self.active[i]
            dummy.P = self.P[i]
            dummy.D = self.D[i]
            dummy.wvl = self.wvl[i]
            dummy.N = self.N[i]
            dummy.normal = self.normal[i]
            
            dummy_list.append(dummy)
        return dummy_list

    def propagate(self, surfaces: List[geometry]) -> None:
        """Propagates the ray group through a list of surfaces.

        Parameters
        ----------
        surfaces : list of geometry objects
            Surfaces to propagate through in order.
        """
        for surface in surfaces:
            self.transform(surface)
            self.find_intersections(surface)
            if not np.any(self.active):
                break
            self.interact(surface)
            if not np.any(self.active):
                break
            self.ray_lab_frame(surface)
            self.update_history()