import numpy as np

from .ray import ray
from .geometry import geometry

from typing import Union, List, Dict

def ray_plane(geo_params: List[Dict],
              pos: Union[List[float], float, int],
              radius: Union[float, int],
              d: List[float],
              nrays: int = 100,
              wvl: Union[float, int] = 0.55) -> List[ray]:
    """Creates a plane of rays and propagates them through geometry.

    Note
    ----
    This is done by creating a square grid with side length 2 x radius
    and then chopping off all points outside of a radius. But this means
    the number of rays is significantly less than the parameter nrays in
    some cases. This needs to be fixed eventually.

    Parameters
    ----------
    geo_params : list of dictionaries
        Dictionaries correspond to surfaces.
    pos : list of length 3 or float/int
        X, Y, and Z values of the plane or just the Z value.
    d : list of length 3
        Direction cosines for each ray.
    radius : float/int
        Radius of ray plane.
    nrays : int
        Number of rays in the ray plane.
    wvl: float/int
        Wavelength of the ray in microns (default 550nm --> 0.55).

    Returns
    -------
    list of ray objects
        Returns rays after propagating them through geometry list.

    """

    x_mesh = np.linspace(-radius, radius, int(4./np.pi*np.sqrt(nrays)))
    y_mesh = np.linspace(-radius, radius, int(4./np.pi*np.sqrt(nrays)))
    x_points, y_points = np.meshgrid(x_mesh, y_mesh)
    xs, ys = x_points.ravel(), y_points.ravel()
    dis = np.sqrt(pow(xs,2) + pow(ys, 2))
    if isinstance(pos, list):
        x_circ, y_circ = xs[dis < radius] + pos[0], ys[dis < radius] + pos[1]
        #To-do: Should transform for d
        z_circ = [pos[2]]*len(x_circ)
    else:
        x_circ, y_circ = xs[dis < radius], ys[dis < radius]
        #To-do: Should transform for d
        z_circ = [pos]*len(x_circ)
    P_arr = np.vstack((x_circ, y_circ, z_circ)).T
    D_arr = np.array([d] *len(x_circ))
    #Initialize rays.
    rays = [ray(params={'D':D, 'P':P, 'wvl':wvl}) for D,P in zip(D_arr, P_arr)]
    geo = [geometry(surf) for surf in geo_params]
    #Propagate rays through geometry.
    for rayiter in rays:
        rayiter.propagate(geo)
    return rays
