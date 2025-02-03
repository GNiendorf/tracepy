import numpy as np

from .ray import RayGroup
from .geometry import geometry

from typing import Union, List, Dict

def ray_plane(geo_params: List[Dict],
              pos: Union[List[float], float, int],
              radius: Union[float, int],
              d: List[float],
              nrays: int = 100,
              wvl: Union[float, int] = 0.55) -> RayGroup:
    """Creates a plane of rays and propagates them through geometry.

    Note
    ----
    This is done by creating a square grid with side length 2*radius
    and then selecting only those points within the circle of given radius.

    Parameters
    ----------
    geo_params : list of dictionaries
        Dictionaries corresponding to surfaces.
    pos : list of length 3 or float/int
        X, Y, and Z values of the plane or just the Z value.
    d : list of length 3
        Direction cosines for each ray.
    radius : float/int
        Radius of ray plane.
    nrays : int
        Approximate number of rays in the ray plane.
    wvl: float/int
        Wavelength of the rays in microns.

    Returns
    -------
    RayGroup
        RayGroup object after propagation through the geometry.
    """
    # Create a square grid and select points within a circle of given radius.
    grid_size = int(4./np.pi*np.sqrt(nrays))  
    x_vals = np.linspace(-radius, radius, grid_size)
    y_vals = np.linspace(-radius, radius, grid_size)
    x_mesh, y_mesh = np.meshgrid(x_vals, y_vals)
    xs = x_mesh.ravel()
    ys = y_mesh.ravel()
    dis = np.sqrt(xs**2 + ys**2)
    mask = dis < radius
    xs = xs[mask]
    ys = ys[mask]
    if isinstance(pos, list):
        x_circ = xs + pos[0]
        y_circ = ys + pos[1]
        z_circ = np.full(xs.shape, pos[2])
    else:
        x_circ = xs
        y_circ = ys
        z_circ = np.full(xs.shape, pos)
    P_arr = np.vstack((x_circ, y_circ, z_circ)).T
    D_arr = np.tile(d, (P_arr.shape[0], 1))
    # Initialize ray group.
    raygroup = RayGroup(P_arr, D_arr, wvl=wvl)
    # Convert geometry dictionaries to geometry objects.
    geo_list = [geometry(surf) for surf in geo_params]
    # Propagate rays through geometry.
    raygroup.propagate(geo_list)
    return raygroup