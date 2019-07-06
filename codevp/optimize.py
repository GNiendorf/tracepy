import numpy as np
from scipy.optimize import least_squares

from .optplot import optplot
from .geometry import geometry
from .raygroup import ray_plane

def update_geometry(inputs, geoparams, vary_dicts):
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
    params_iter = update_geometry(inputs, geoparams, vary_dicts)
    raygroup_iter = ray_plane(params_iter, [0., 0., 0.], 1.1, d=[0.,0.,1.], nrays=100)
    try:
        oplt = optplot(params_iter, raygroup_iter)
        rms = oplt.spotdiagram(optimizer=True)
    except:
        rms = 999.
    return rms

def optimize(geoparams, raygroup, vary_dicts, typeof='least_squares'):
    geo_list = [geometry(geo) for geo in geoparams]
    initial_guess = []
    for dict_ in vary_dicts:
        name = dict_["name"]
        vary_list = dict_["vary"]
        for surface in geo_list:
            if surface["name"] == name:
                for idx, item in enumerate(vary_list):
                    attr = surface[item]
                    if item == 'P':
                        attr = attr[2]
                    initial_guess.append(attr)

    if typeof == 'least_squares':
        res = least_squares(get_rms, initial_guess, args=(geoparams, raygroup, vary_dicts))
    new_params = update_geometry(res.x, geoparams, vary_dicts)
    return new_params