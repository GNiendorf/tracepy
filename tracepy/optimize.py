# Authors: Gavin Niendorf <gavinniendorf@gmail.com>
#
# Functions for optimizing an optical system.
#
# License: MIT

import numpy as np
from scipy.optimize import least_squares

from .optplot import spotdiagram
from .raygroup import ray_plane
from .exceptions import TraceError

def update_geometry(inputs, geoparams, vary_dicts):
    """Return the geometry requested by the optimization algorithm.

    Parameters
    ----------
    inputs : list of floats/ints
        Values for varying parameters that the optimizer selected.
    geoparams : list of dictionaries
        Surface parameters in dictionary form.
    vary_dicts : list of dictionaries
        Lists the parameters that the user wants to optimize.

    Returns
    -------
    geoparams : list of dictionaries
        Returns the list of paramaters after updating it with inputs.

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
    """Return the rms of an updated geometry.

    Note
    ----
    If no rays survive then a large rms value is chosen
    so that the optimizer is discouraged from looking
    there.

    Parameters
    ----------
    inputs : list of floats/ints
        Values for varying parameters that the optimizer selected.
    geoparams : list of dictionaries
        Surface parameters in dictionary form.
    vary_dicts : list of dictionaries
        Lists the parameters that the user wants to optimize.

    Returns
    -------
    rms : float
        RMS of the spotdiagram.

    """

    params_iter = update_geometry(inputs, geoparams, vary_dicts)
    raygroup_iter = ray_plane(params_iter, [0., 0., 0.], 1.1, d=[0.,0.,1.], nrays=50)
    ratio_surv = np.sum([1 for ray in raygroup_iter if ray.P is not None])/len(raygroup_iter)
    try:
        rms = spotdiagram(params_iter, raygroup_iter, optimizer=True)
    except TraceError:
        rms = 999.
    #Weight of failed propagation.
    surv_const = 100
    return rms + (1-ratio_surv)*surv_const

def optimize(geoparams, vary_dicts, typeof='least_squares', max_iter=None):
    """Optimize a given geometry for a given varylist and return the new geometry.

    Parameters
    ----------
    geoparams : list of dictionaries
        Surface parameters in dictionary form.
    vary_dicts : list of dictionaries
        Lists the parameters that the user wants to optimize.
    typeof : str
        Type of optimization to use.
    max_iter : int
        Max number of function calls until the optimizer stops.

    Returns
    -------
    new_params : list of dictionaries
        List of otimized dictionaries that describe surfaces.

    """

    initial_guess = []
    param_lb = []
    param_ub = []
    for dict_ in vary_dicts:
        name = dict_["name"]
        vary_list = dict_["vary"]
        for surface in geoparams:
            if "name" in surface and surface["name"] == name:
                for item in vary_list:
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
