import numpy as np
from scipy.optimize import least_squares

from .optplot import spotdiagram
from .geometry import geometry
from .raygroup import ray_plane

def update_geometry(inputs, geoparams, vary_dicts):
    """
        Return the geometry requested by the optimization algorithm.
    """
    vary_idxs = 0
    for dict_ in vary_dicts:
        name = dict_["name"]
        vary_list = dict_["vary"]
        for surface in geoparams:
            if "name" in surface and surface["name"] == name:
                for item in vary_list:
                    surface[item] = inputs[vary_idxs]
                    vary_idxs += 1
    return geoparams

def get_rms(inputs, geoparams, vary_dicts):
    """ 
        Return the rms of an updated geometry.
    """
    params_iter = update_geometry(inputs, geoparams, vary_dicts)
    raygroup_iter = ray_plane(params_iter, [0., 0., 0.], 1.1, d=[0.,0.,1.], nrays=50)
    try:
        rms = spotdiagram(params_iter, raygroup_iter, optimizer=True)
    except:
        rms = 999.
    return rms

def optimize(geoparams, vary_dicts, typeof='least_squares', max_iter=None):
    """ 
        Optimize a given geometry for a given varylist and return the new geometry. 
    """
    initial_guess = []
    param_lb = []
    param_ub = []
    for dict_ in vary_dicts:
        name = dict_["name"]
        vary_list = dict_["vary"]
        for surf_idx, surface in enumerate(geoparams):
            if "name" in surface and surface["name"] == name:
                for idx, item in enumerate(vary_list):
                    attr = surface[item]
                    if item == 'P' and (isinstance(item, float) or isinstance(item, int)):
                         attr = attr[2]
                    if item == 'c' or item == 'kappa':
                        lower_bound = 0 if attr > 0 else -np.inf
                        upper_bound = np.inf if attr >= 0 else 0.
                        param_lb.append(lower_bound)
                        param_ub.append(upper_bound)
                    else:
                        param_lb.append(-np.inf)
                        param_ub.append(np.inf)
                    initial_guess.append(attr)
    param_bounds = [tuple(param_lb), tuple(param_ub)]
    if typeof == 'least_squares':
        res = least_squares(get_rms, initial_guess, args=(geoparams, vary_dicts), max_nfev=max_iter,\
                            bounds = param_bounds) 
    new_params = update_geometry(np.around(res.x, 3), geoparams, vary_dicts)
    return new_params