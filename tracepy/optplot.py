# Authors: Gavin Niendorf <gavinniendorf@gmail.com>
#
# Functions for plotting useful features of an optical system.
#
# License: MIT

import numpy as np

import matplotlib.pyplot as plt

from .geometry import geometry
from .transforms import transform
from .exceptions import TraceError

def _gen_object_points(surface, surface_idx, rays):
    """Transform intersection points into a surfaces' reference frame.

    Parameters
    ----------
    surface : geometry object
        Surface whos reference frame the points will be transformed into.
    surface_idx : int
        Integer corresponding to where the surface is in the propagation
        order. For example, 0 means the surface is the first surface rays
        are propagated through.
    rays : list of ray objects
        Rays that were propagated through the geometry list.

    Returns
    -------
    X, Y : np.array of len(rays)
        X, Y pair in the surface's reference frame.
    points_obj: 2d np.array
        X, Y pair points in 2d array for easy rms calculation.

    """

    points = np.array([rayiter.P_hist[surface_idx] for rayiter in rays if rayiter.P is not None])
    if points.size == 0:
        #No rays survived
        raise TraceError()
    #Get X,Y points in obj. reference frame.
    points_obj = transform(surface.R, surface, points)
    #Round arrays to upper bound on accuracy.
    points_obj = np.around(points_obj, 14)
    if points_obj.ndim == 2:
        X, Y = points_obj[:,0], points_obj[:,1]
    elif points_obj.ndim == 1:
        X, Y = points_obj[0], points_obj[1]
        points_obj = np.array([points_obj])
    return X, Y, points_obj

def spotdiagram(geo_params, rays, pltparams = {'color': 'red'}, optimizer=False):
    """Plots the transformed intersection points of rays on the stop surface.

    Note
    ----
    Directly plots the spot diagram rather than returning something.

    Parameters
    ----------
    geo_params : list of dictionaries
        Surfaces in order of propagation.
    rays : list of ray objects
        Rays that were propagated through the geometry list.
    pltparams : dictionary
        Plotting attributes of the spot diagram.
    optimizer : bool
        Flag for whether the optimizer is calling the function to get
        the rms of the spot diagram.

    """

    stop = geometry(geo_params[-1])
    X, Y, points_obj = _gen_object_points(stop, -1, rays)
    rms = np.std(points_obj[:,[0,1]] - points_obj[:,[0,1]].mean(axis=0))
    if optimizer:
        return rms
    plt.subplot(1,1,1, aspect='equal')
    plt.locator_params(axis='x', nbins=8)
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    plt.plot(X, Y, '+', **pltparams)
    plt.text(0.65, 1.015, 'RMS=%.2E' % rms, transform=plt.gca().transAxes)

def plotobject(geo_params, rays, pltparams = {'color': 'blue'}):
    """Plots the initial ray locations in the initial objects' frame.

    Note
    ----
    Directly plots the object spot diagram rather than returning something.

    Parameters
    ----------
    geo_params : list of dictionaries
        Surfaces in order of propagation.
    rays : list of ray objects
        Rays that were propagated through the geometry list.
    pltparams : dictionary
        Plotting attributes of the object spot diagram.

    """

    start = geometry(geo_params[0])
    X, Y, points_obj = _gen_object_points(start, 0, rays)
    rms = np.std(points_obj[:,[0,1]] - points_obj[:,[0,1]].mean(axis=0))
    plt.subplot(1,1,1, aspect='equal')
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    plt.text(0.65, 1.015, 'RMS=%.2E' % rms, transform=plt.gca().transAxes)
    plt.plot(X, Y, '+', **pltparams)

def rayaberration(geo_params, rays, pltparams= {'color': 'red', 'linewidth': 2}):
    """Plots the ray aberration diagrams for an optical system.

    Note
    ----
    Directly plots the ray aberration curve rather than returning something.

    Parameters
    ----------
    geo_params : list of dictionaries
        Surfaces in order of propagation.
    rays : list of ray objects
        Rays that were propagated through the geometry list.
    pltparams : dictionary
        Plotting attributes of the curve.

    """

    start = geometry(geo_params[0])
    stop = geometry(geo_params[-1])
    X_start, Y_start, start_points = _gen_object_points(start, 0, rays)
    X_stop, Y_stop, _ = _gen_object_points(stop, -1, rays)
    #Chief ray information
    distances = np.sqrt(np.sum((start_points[:,[0,1]] - start_points[:,[0,1]].mean(axis=0))**2, axis=1))
    chief_ray_idx = np.argmin(distances)
    chief_ray = start_points[:,[0,1]][chief_ray_idx]
    #Rays along principal plane
    ray_subset_idxs = np.where(X_start == chief_ray[0])
    rel_stops = Y_stop[ray_subset_idxs] - Y_stop[chief_ray_idx]
    #Normalized pupil coordinates
    norm_coords = Y_start[ray_subset_idxs] - chief_ray[1]
    norm_coords[np.where(norm_coords < 0)] /= abs(np.min(norm_coords))
    norm_coords[np.where(norm_coords > 0)] /= abs(np.max(norm_coords))
    #Used to show axes when there is no aberration
    if np.all(rel_stops == 0.):
        plt.plot([0., 0.], [-.04, .04], 'k')
    else:
        plt.plot([0]*len(rel_stops), rel_stops, 'k')
    plt.plot(norm_coords, [0]*len(norm_coords), 'k')
    plt.plot(norm_coords, rel_stops, **pltparams)
