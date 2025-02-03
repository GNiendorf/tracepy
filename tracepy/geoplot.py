import numpy as np
import matplotlib.pyplot as plt

from .ray import RayGroup
from .geometry import geometry
from .transforms import lab_frame_points

from typing import List, Dict, Optional, Any

def _gen_points(surfaces: List[geometry]) -> np.ndarray:
    """Generates the mesh points for each surface in the obj frame.

    Parameters
    ----------
    surfaces : list of geometry objects
        Surfaces in propagation order to plot.

    Returns
    -------
    surfpoints : 2d np.array
        Points (X, Y, Z) in the surface's reference frame.
    """
    surfpoints = []
    for surface in surfaces:
        bound = surface.Diam / 2.
        # General mesh points.
        linspace = np.linspace(-bound, bound, 200)
        x_mesh, y_mesh = np.meshgrid(linspace, -linspace)
        # Used for cross sections.
        x_points = np.append(x_mesh, [linspace, np.zeros(200)])
        y_points = np.append(y_mesh, [np.zeros(200), linspace])
        meshpoints_2d = np.vstack((x_points.ravel(), y_points.ravel())).T
        # Function values in obj. frame.
        z_points = surface.get_surface_plot(meshpoints_2d)
        meshpoints = np.vstack((x_points.ravel(), y_points.ravel(), z_points)).T
        # Append to surfpoints which holds all surfaces' points.
        surfpoints.append(meshpoints)
    return np.array(surfpoints)

def _plot_rays(rays: np.ndarray,
               axes: List[int],
               pltparams: Dict) -> None:
    """Plots 2d ray history points. Takes list axes to specify axes (0, 1, 2) to plot.

    Parameters
    ----------
    rays : np.array
        Array of dummy ray objects to be plotted.
    axes : list of length 2 with integers from range [0,2]
        Axes (X, Y, Z) to plot from ray points.
    pltparams : dictionary
        Plot characteristics of rays such as colors and alpha.
    """
    for ray in rays:
        if len(ray.P_hist) < 2:
            continue
        for idx in range(len(ray.P_hist) - 1):
            F, G = ray.P_hist[idx][axes]
            F_p, G_p = ray.P_hist[idx + 1][axes]
            #Alpha, beta, or gamma rotations.
            H_p, I_p = ray.D_hist[idx + 1][axes]
            plt.plot([G, G_p], [F, F_p], **pltparams)
        # Only plot final direction if the loop ran at least once.
        if len(ray.P_hist) >= 2:
            plt.plot([G_p, G_p + I_p], [F_p, F_p + H_p], **pltparams)

def _clip_lens(surfaces: List[geometry], surfpoints: np.ndarray, idx: int) -> np.ndarray:
    """Clips points outside of a lens intersection point.

    Parameters
    ----------
    surfaces : list of geometry objects
        Surfaces whose reference frames to transform from.
    surfpoints : 2d np.array
        Points from surface for each row that will be clipped.
    idx : int
        Index of surface in propagation order.

    Returns
    -------
    surfpoints : 2d np.array
        Points after clipping.
    """
    surf1, surf2 = surfaces[idx], surfaces[idx + 1]
    d = np.sqrt(np.sum(np.square(surf1.P - surf2.P)))
    points1, points2 = np.nan_to_num(surfpoints[idx]), np.nan_to_num(surfpoints[idx + 1])
    points2[:, 2] += d
    clipped_idx = (points2[:, 2] - points1[:, 2]) <= 0.
    surfpoints[idx][:, 2][clipped_idx] = np.nan
    surfpoints[idx + 1][:, 2][clipped_idx] = np.nan
    return surfpoints

def _plot_surfaces(geo_params: List[Dict], axes: List[int]) -> None:
    """Plots 2d surface cross sections. Takes list axes to specify axes (0, 1, 2) to plot surfaces in.

    Note
    ----
    Directly plots the surfaces rather than returning something.

    Parameters
    ----------
    geo_params : list of dictionaries
        Surfaces in propagation order to plot.
    axes : list of length 2 with integers from range [0,2]
        Axes (X, Y, Z) to plot surfaces in.
    """
    surfaces = [geometry(surf) for surf in geo_params]
    surfpoints = _gen_points(surfaces)
    lens_check = 0
    start = None
    for idx, surf in enumerate(surfaces):
        lens_condition = (idx + 1 < len(surfaces) and
                          surfaces[idx].action == surfaces[idx + 1].action == 'refraction')
        if lens_condition:
            surfpoints = _clip_lens(surfaces, surfpoints, idx)
        with np.errstate(invalid='ignore'):
            if np.any(np.mod(surf.D / np.pi, 1) != 0) and surf.c == 0 and surf.diam == 0:
                # Find cross section points.
                cross_idx = abs(surfpoints[idx][:, axes[1]]) == 0
            else:
                # Find cross section points.
                cross_idx = abs(surfpoints[idx][:, 1 - axes[0]]) == 0
        cross_points = surfpoints[idx][cross_idx]
        # Transform to lab frame.
        points = lab_frame_points(surf.R, surf, cross_points)
        F, G = points[:, axes[0]], points[:, axes[1]]
        # Connect the surfaces in a lens
        if surfaces[idx].action == surfaces[idx - 1].action == 'refraction' and start is not None:
            lens_check = 1 - lens_check
            if lens_check == 1:
                start_2 = np.array([F[0], G[0]])
                end_2 = np.array([F[-1], G[-1]])
                dis1 = np.sqrt(np.sum(np.square(start - start_2)))
                dis2 = np.sqrt(np.sum(np.square(start - end_2)))
                if dis1 <= dis2:
                    idxs = [0, -1]
                else:
                    idxs = [-1, 0]
                F = np.insert(F, idxs, [start[0], end[0]])
                G = np.insert(G, idxs, [start[1], end[1]])
        if lens_condition:
            start = np.array([F[0], G[0]])
            end = np.array([F[-1], G[-1]])
        plt.plot(G, F, 'k')

def plotxz(geo_params: List[Dict],
           ray_list: Any,
           pltparams: Dict = {'c': 'red', 'alpha': 0.3},
           both: Optional[bool] = None) -> None:
    """Plots the xz coordinates of all rays and surface cross sections.

    Parameters
    ----------
    geo_params : list of dictionaries
        Surfaces in propagation order to plot.
    ray_list : list or RayGroup object
        Dummy ray objects or RayGroup that are going to be plotted.
    pltparams : dictionary
        Plot characteristics of rays such as colors and alpha.
    both (optional) : bool
        Flag for overriding a 111 subplot for self.plot2d.
    """
    if hasattr(ray_list, "to_ray_list"):
        ray_list = ray_list.to_ray_list()
    rays = np.array([rayiter for rayiter in ray_list if rayiter.active])
    if both is None:
        plt.subplot(1, 1, 1, aspect='equal')
    _plot_rays(rays, [0, 2], pltparams)
    _plot_surfaces(geo_params, axes=[0, 2])
    plt.xlabel("Z")
    plt.ylabel("X")

def plotyz(geo_params: List[Dict],
           ray_list: Any,
           pltparams: Dict = {'c': 'red', 'alpha': 0.3},
           both: Optional[bool] = None) -> None:
    """Plots the yz coordinates of all rays and surface cross sections.

    Parameters
    ----------
    geo_params : list of dictionaries
        Surfaces in propagation order to plot.
    ray_list : list or RayGroup object
        Dummy ray objects or RayGroup that are going to be plotted.
    pltparams : dictionary
        Plot characteristics of rays such as colors and alpha.
    both (optional) : bool
        Flag for overriding a 111 subplot for self.plot2d.
    """
    if hasattr(ray_list, "to_ray_list"):
        ray_list = ray_list.to_ray_list()
    rays = np.array([rayiter for rayiter in ray_list if rayiter.active])
    if both is None:
        plt.subplot(1, 1, 1, aspect='equal')
    _plot_rays(rays, [1, 2], pltparams)
    _plot_surfaces(geo_params, axes=[1, 2])
    plt.xlabel("Z")
    plt.ylabel("Y")

def plot2d(geo_params: List[Dict],
           rays: Any,
           pltparams: Dict = {'c': 'red', 'alpha': 0.3}) -> None:
    """Plots both xz and yz side-by-side.

    Parameters
    ----------
    geo_params : list of dictionaries
        Surfaces in propagation order to plot.
    rays : list or RayGroup object
        Dummy ray objects or RayGroup that are going to be plotted.
    pltparams : dictionary
        Plot characteristics of rays such as colors and alpha.
    """
    if hasattr(rays, "to_ray_list"):
        rays = rays.to_ray_list()
    plt.subplot(2, 1, 1, aspect='equal')
    plotxz(geo_params, rays, pltparams, both=True)
    plt.subplot(2, 1, 2, aspect='equal')
    plotyz(geo_params, rays, pltparams, both=True)
    plt.tight_layout()