import numpy as np
from scipy.optimize import least_squares

from .optplot import optplot
from .geometry import geometry
from .raygroup import ray_plane

def update_geometry(inputs, geoparams, vary_dicts):
    """ Return the geometry requested by the optimization algorithm. """
    for dict_ in vary_dicts:
        name = dict_["name"]
        vary_list = dict_["vary"]
        for surface in geoparams:
            try:
                if surface["name"] == name:
                    for idx, item in enumerate(vary_list):
                        surface[item] = inputs[idx]
            except:
                continue
    return geoparams

def get_rms(inputs, geoparams, raygroup, vary_dicts):
    """ Return the rms of an updated geometry. """
    params_iter = update_geometry(inputs, geoparams, vary_dicts)
    raygroup_iter = ray_plane(params_iter, [0., 0., 0.], 1.1, d=[0.,0.,1.], nrays=100)
    try:
        oplt = optplot(params_iter, raygroup_iter)
        rms = oplt.spotdiagram(optimizer=True)
    except:
        rms = 999.
    return rms

def optimize(geoparams, raygroup, vary_dicts, typeof='least_squares', max_iter=None):
    """ Optimize a given geometry for a given raygroup and varylist and return the new geometry. """
    initial_guess = []
    for dict_ in vary_dicts:
        name = dict_["name"]
        vary_list = dict_["vary"]
        for surface in geoparams:
            try:
                if surface["name"] == name:
                    for idx, item in enumerate(vary_list):
                        attr = surface[item]
                        if item == 'P' and (isinstance(item, float) or isinstance(item, int)):
                            attr = attr[2]
                        initial_guess.append(attr)
            except:
                continue
    if typeof == 'least_squares':
        res = least_squares(get_rms, initial_guess, args=(geoparams, raygroup, vary_dicts), max_nfev=max_iter)
    new_params = update_geometry(res.x, geoparams, vary_dicts)
    return new_params