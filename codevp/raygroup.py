import numpy as np
from numpy import cos, sin, pi, sqrt

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
        z_circ = [pos[2]]*len(x_circ) #To-do: Should transform for d
    else:
        x_circ, y_circ = xs[dis < radius], ys[dis < radius]
        z_circ = [pos]*len(x_circ) #To-do: Should transform for d
    P_arr = np.vstack((x_circ, y_circ, z_circ)).T
    D_arr = np.array([d] *len(x_circ)) 
    rays = [ray(params={'D':D, 'P':P}) for D,P in zip(D_arr, P_arr)] #Initialize rays.
    geo = [geometry(surf) for surf in geo_params]
    for rayiter in rays: #Propagate rays through geometry.
        rayiter.propagate(geo)
    return rays

def ray_fan(geo_params, pos, radius, d, ang=None, nrays=1000):
    """ Creates a fan of rays and propagates them through geometry. """
    #To-do: Need to account for d and angle of fan
    D_arr = 10*np.random.rand(nrays, 3)
    D_arr[:,1]*=np.random.uniform(-1,1,nrays)
    D_arr[:,0]*=np.random.uniform(-1,1,nrays)
    if isinstance(pos, list):
        P_arr = [pos] * (2*nrays)
    else:
        P_arr = np.zeros((2*nrays,3))
        P_arr[:,2] = pos
    rays = [ray(params={'D':D, 'P':P}) for D,P in zip(D_arr, P_arr)] #Initialize rays.
    geo = [geometry(surf) for surf in geo_params]
    for rayiter in rays: #Propagate rays through geometry.
        rayiter.propagate(geo)
    return rays
