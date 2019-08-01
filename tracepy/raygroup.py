import numpy as np
from numpy import pi, sqrt

from .ray import ray
from .geometry import geometry

def ray_plane(geo_params, pos, radius, d, ang=None, nrays=1000):
    """ Creates a plane of rays and propagates them through geometry. """
    x_mesh = np.linspace(-radius, radius, int(4./pi*sqrt(nrays)))
    y_mesh = np.linspace(-radius, radius, int(4./pi*sqrt(nrays)))
    x_points, y_points = np.meshgrid(x_mesh, y_mesh)
    xs, ys = x_points.ravel(), y_points.ravel()
    dis = sqrt(pow(xs,2) + pow(ys, 2))
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
    rays = [ray(params={'D':D, 'P':P}) for D,P in zip(D_arr, P_arr)]
    geo = [geometry(surf) for surf in geo_params]
    #Propagate rays through geometry.
    for rayiter in rays:
        rayiter.propagate(geo)
    return rays
