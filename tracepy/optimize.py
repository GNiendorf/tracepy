import numpy as np
from scipy.optimize import least_squares

from .optplot import spot_rms
from .raygroup import ray_plane
from .constants import SURV_CONST, MAX_RMS
from .exceptions import TraceError

from typing import List, Union, Dict, Optional

def update_geometry(inputs: List[Union[float, int]],
                    geoparams: List[Dict],
                    vary_dicts: List[Dict]) -> List[Dict]:
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

def get_rms(inputs: List[Union[float, int]],
            geoparams: List[Dict],
            vary_dicts: List[Dict]) -> float:
    """Return the RMS of an updated geometry.

    Note
    ----
    If no rays survive then a large RMS value is chosen so that the optimizer is discouraged from looking
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
    raygroup_iter = ray_plane(params_iter, [0., 0., 0.], 1.1, d=[0., 0., 1.], nrays=50)
    rays_list = raygroup_iter.to_ray_list()
    ratio_surv = np.sum([1 for ray in rays_list if ray.active]) / len(rays_list)
    try:
        rms = spot_rms(params_iter, rays_list)
    except TraceError:
        rms = MAX_RMS
    return rms + (1 - ratio_surv) * SURV_CONST

def optimize(geoparams:List[Dict],
             vary_dicts: List[Dict],
             typeof: str = 'least_squares',
             max_iter: Optional[int] = None) -> List[Dict]:
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
