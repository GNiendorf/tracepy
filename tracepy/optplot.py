import numpy as np
import matplotlib.pyplot as plt

from .geometry import geometry
from .transforms import transform_points
from .constants import PLOT_ROUNDING_ACC
from .exceptions import TraceError

from typing import List, Dict, Tuple, Any

def _gen_object_points(surface: geometry,
                       surface_idx: int,
                       rays: List[Any]) -> np.ndarray:
    """Transform intersection points into a surface's reference frame.

    Parameters
    ----------
    surface : geometry object
        Surface whose reference frame the points will be transformed into.
    surface_idx : int
        Integer corresponding to where the surface is in the propagation
        order. For example, 0 means the surface is the first surface rays
        are propagated through.
    rays : list of ray objects
        Rays that were propagated through the geometry list.

    Returns
    -------
    points_obj: 2d np.array
        X, Y pair points in 2d array for easy RMS calculation.
    """
    # Gather the history points for all active rays.
    points = np.array([rayiter.P_hist[surface_idx] for rayiter in rays if rayiter.active])
    if points.size == 0:
        # Instead of raising an error, return a single point at (0,0).
        return np.array([[0.0, 0.0]])
    # Transform to the object (surface) frame.
    points_obj = transform_points(surface.R, surface, points)
    points_obj = np.around(points_obj, PLOT_ROUNDING_ACC)
    # Filter out any rows that contain NaN.
    valid = ~np.isnan(points_obj).any(axis=1)
    if np.sum(valid) == 0:
        return np.array([[0.0, 0.0]])
    return points_obj[valid]

def calculate_rms(points: np.ndarray) -> float:
    """Calculates the RMS of the given points.

    Parameters
    ----------
    points : np.ndarray
        Array of points to calculate RMS for.

    Returns
    -------
    float
        The calculated RMS value.
    """
    # If only one valid point exists, the RMS is zero.
    if points.shape[0] < 2:
        return 0.0
    return np.std(points[:, [0, 1]] - points[:, [0, 1]].mean(axis=0))

def spot_rms(geo_params: List[Dict], rays: Any) -> float:
    """Calculates the RMS of the spot diagram points.

    Parameters
    ----------
    geo_params : list of dictionaries
        Surfaces in order of propagation.
    rays : list of ray objects or RayGroup object
        Rays that were propagated through the geometry list. If a RayGroup
        object is provided, it will be converted automatically.

    Returns
    -------
    float
        The calculated RMS value.
    """
    if hasattr(rays, "to_ray_list"):
        rays = rays.to_ray_list()
    stop = geometry(geo_params[-1])
    points_obj = _gen_object_points(stop, -1, rays)
    return calculate_rms(points_obj)

def spotdiagram(geo_params: List[Dict],
                rays: Any,
                pltparams: Dict = {'color': 'red'}) -> None:
    """Plots the transformed intersection points of rays on the stop surface.

    Parameters
    ----------
    geo_params : list of dictionaries
        Surfaces in order of propagation.
    rays : list of ray objects or RayGroup object
        Rays that were propagated through the geometry list. If a RayGroup
        object is provided, it will be converted automatically.
    pltparams : dictionary
        Plotting attributes of the spot diagram.
    """
    if hasattr(rays, "to_ray_list"):
        rays = rays.to_ray_list()
    stop = geometry(geo_params[-1])
    points_obj = _gen_object_points(stop, -1, rays)
    X, Y = points_obj[:, 0], points_obj[:, 1]
    rms = calculate_rms(points_obj)

    plt.subplot(1, 1, 1, aspect='equal')
    plt.locator_params(axis='x', nbins=8)
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0, 0))
    plt.plot(X, Y, '+', **pltparams)
    plt.text(0.65, 1.015, 'RMS=%.2E' % rms, transform=plt.gca().transAxes)

def plotobject(geo_params: List[Dict],
               rays: Any,
               pltparams: Dict = {'color': 'blue'}) -> None:
    """Plots the initial ray locations in the initial object's frame.

    Note
    ----
    Directly plots the object spot diagram rather than returning something.

    Parameters
    ----------
    geo_params : list of dictionaries
        Surfaces in order of propagation.
    rays : list of ray objects or RayGroup object
        Rays that were propagated through the geometry list. If a RayGroup
        object is provided, it will be converted automatically.
    pltparams : dictionary
        Plotting attributes of the object spot diagram.
    """
    if hasattr(rays, "to_ray_list"):
        rays = rays.to_ray_list()
    start = geometry(geo_params[0])
    points_obj = _gen_object_points(start, 0, rays)
    X, Y = points_obj[:, 0], points_obj[:, 1]
    rms = np.std(points_obj[:, [0, 1]] - points_obj[:, [0, 1]].mean(axis=0))
    plt.subplot(1, 1, 1, aspect='equal')
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0, 0))
    plt.text(0.65, 1.015, 'RMS=%.2E' % rms, transform=plt.gca().transAxes)
    plt.plot(X, Y, '+', **pltparams)

def rayaberration(geo_params: List[Dict],
                  rays: Any,
                  pltparams: Dict = {'color': 'red', 'linewidth': 2}) -> None:
    """Plots the ray aberration diagrams for an optical system.

    Note
    ----
    Directly plots the ray aberration curve rather than returning something.

    Parameters
    ----------
    geo_params : list of dictionaries
        Surfaces in order of propagation.
    rays : list of ray objects or RayGroup object
        Rays that were propagated through the geometry list. If a RayGroup
        object is provided, it will be converted automatically.
    pltparams : dictionary
        Plotting attributes of the curve.
    """
    if hasattr(rays, "to_ray_list"):
        rays = rays.to_ray_list()
    start = geometry(geo_params[0])
    start_points = _gen_object_points(start, 0, rays)
    X_start, Y_start = start_points[:, 0], start_points[:, 1]
    stop = geometry(geo_params[-1])
    stop_points = _gen_object_points(stop, -1, rays)
    X_stop, Y_stop = stop_points[:, 0], stop_points[:, 1]
    # Chief ray information
    distances = np.sqrt(np.sum((start_points[:, [0, 1]] - start_points[:, [0, 1]].mean(axis=0))**2, axis=1))
    chief_ray_idx = np.argmin(distances)
    chief_ray = start_points[:, [0, 1]][chief_ray_idx]
    ray_subset_idxs = np.where(X_start == chief_ray[0])
    rel_stops = Y_stop[ray_subset_idxs] - Y_stop[chief_ray_idx]
    # Normalized pupil coordinates
    norm_coords = Y_start[ray_subset_idxs] - chief_ray[1]
    norm_coords[np.where(norm_coords < 0)] /= abs(np.min(norm_coords))
    norm_coords[np.where(norm_coords > 0)] /= abs(np.max(norm_coords))
    # Used to show axes when there is no aberration
    if np.all(rel_stops == 0.):
        plt.plot([0., 0.], [-0.04, 0.04], 'k')
    else:
        plt.plot([0] * len(rel_stops), rel_stops, 'k')
    plt.plot(norm_coords, [0] * len(norm_coords), 'k')
    plt.plot(norm_coords, rel_stops, **pltparams)