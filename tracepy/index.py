import os
import json
import yaml

import numpy as np
from scipy.optimize import curve_fit

def cauchy_two_term(x: float, B: float, C: float) -> float:
    """
    This is a simple two-term Cauchy equation for the index of refraction as a function of wavelength.
    https://en.wikipedia.org/wiki/Cauchy%27s_equation
    It is used to create a function for refractive index database entries when only tabulated data is available.
    
    Parameters
    ----------
    x : float
        Wavelength (in micrometers).
    B : float
        First term of the Cauchy equation.
    C : float
        Second term of the Cauchy equation.
        
    Returns
    -------
    float
        Refractive index at the given wavelength.
    """
    return B + (C / (x ** 2))

# TODO: This function, written by @MikeMork, needs to be broken up into several smaller functions with explicit typing.
def glass_index(glass):
    '''
    Given a glass name this function will return a function that takes wavelength in microns and returns index of refraction.
    - All values were taken from refractiveindexinfo's database of different optical glasses.

    https://refractiveindex.info/

    A constant value can also be input to this function to allow directly setting index

    Some glass names like F2 are shared by multiple manufacturers.
    Example input specifying the manufacturer: "F2 Schott"
    Example input without manufacturer: "F2"

    If manufacturer is not specified the following priority are used to choose the glass:
    1)schott
    2)ohara
    3)sumita
    4)cdgm
    5)hoya
    6)hikari
    7)vitron
    8)ami
    9)barberini
    10)corning
    11)lightpath
    12)lzos
    13)misc
    14)nsg

    This list isn't really a comprehensive list of manufacturers but rather the directory structure of the refractiveindexinfo yaml database.

    The functions that this function returns take a wavelength in microns and return the index of refractive index.
    '''
    # Load glass dictionary
    with open(os.path.join(os.path.dirname(__file__), "glass", "glass_dict.json"), "r") as f:
        glass_dict = json.load(f)
    glass_catalog = ["schott", "ohara", "sumita", "cdgm", "hoya", "hikari",
                     "vitron", "ami", "barberini", "lightpath", "lzos", "misc", "nsg"]

    # Process input and return glass file:
    if isinstance(glass, (float, int)):
        constant_index = float(glass)
        return lambda x=0.55: constant_index  # scalar; broadcasting works fine
    elif isinstance(glass, str):
        input_list = glass.split()
        glass_name = input_list[0]
        glass_manufacturer = input_list[1].lower() if len(input_list) > 1 else None

        if glass_name in glass_dict:
            if glass_manufacturer:
                if glass_manufacturer in glass_dict[glass_name]:
                    glass_file = glass_dict[glass_name][glass_manufacturer]
                else:
                    raise Exception(f"{glass_name} is not in the glass catalog for {glass_manufacturer}.")
            else:
                available = list(glass_dict[glass_name].keys())
                # choose based on the defined priority order
                priority = sorted(available, key=lambda x: glass_catalog.index(x) if x in glass_catalog else np.inf)
                glass_file = glass_dict[glass_name][priority[0]]
        else:
            raise Exception(f"{glass_name} is not in the glass catalog.")

        # Process glass file and return the correct function (wavelength in microns)
        split_path = glass_file.split("\\")
        if len(split_path) == 4:
            glass_file_path = os.path.join(split_path[1], split_path[2], split_path[3])
        else:
            glass_file_path = os.path.join(split_path[1], split_path[2], split_path[3], split_path[4])
        f = yaml.load(open(os.path.join(os.path.dirname(__file__), glass_file_path)), Loader=yaml.BaseLoader)
        type_index_function = f["DATA"][0]["type"]
        
        # Use np.asarray(x) so that the operations are vectorized
        if type_index_function == 'formula 1':
            coeffs = [float(x) for x in f["DATA"][0]['coefficients'].split()]
            return lambda x=0.55: np.sqrt(
                1 + coeffs[0] +
                coeffs[1] / (1 - (coeffs[2] / np.asarray(x))**2) +
                coeffs[3] / (1 - (coeffs[4] / np.asarray(x))**2)
            )
        elif type_index_function == 'formula 2':
            coeffs = [float(x) for x in f["DATA"][0]['coefficients'].split()]
            if len(coeffs) >= 7:
                return lambda x=0.55: np.sqrt(
                    1 +
                    coeffs[1] / (1 - coeffs[2] / np.asarray(x)**2) +
                    coeffs[3] / (1 - coeffs[4] / np.asarray(x)**2) +
                    coeffs[5] / (1 - coeffs[6] / np.asarray(x)**2)
                )
            else:
                return lambda x=0.55: np.sqrt(
                    1 +
                    coeffs[1] / (1 - coeffs[2] / np.asarray(x)**2) +
                    coeffs[3] / (1 - coeffs[4] / np.asarray(x)**2)
                )
        elif type_index_function == 'formula 3':
            coeffs = [float(x) for x in f["DATA"][0]['coefficients'].split()]
            return lambda x=0.55: np.sqrt(
                coeffs[0] -
                coeffs[1] * np.asarray(x)**coeffs[2] +
                coeffs[3] * np.asarray(x)**coeffs[4] +
                coeffs[5] * np.asarray(x)**coeffs[6] +
                coeffs[7] * np.asarray(x)**coeffs[8] +
                coeffs[9] * np.asarray(x)**coeffs[10]
            )
        elif type_index_function == 'formula 5':
            coeffs = [float(x) for x in f["DATA"][0]['coefficients'].split()]
            return lambda x=0.55: (coeffs[0] -
                                     coeffs[1] * np.asarray(x)**coeffs[2] +
                                     coeffs[3] * np.asarray(x)**(-coeffs[4]))
        elif type_index_function in ['tabulated n', 'tabulated nk']:
            raw_data = f["DATA"][0]['data'].strip().split('\n')
            wavelength = np.array([float(line.split()[0]) for line in raw_data])
            index = np.array([float(line.split()[1]) for line in raw_data])
            popt, _ = curve_fit(cauchy_two_term, wavelength, index)
            return lambda x=0.55: popt[0] + popt[1] / (np.asarray(x)**2)
    else:
        raise TypeError("Glass must be either a float/int or a string.")