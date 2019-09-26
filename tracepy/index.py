import os

from scipy.optimize import curve_fit

import json
import yaml
import numpy as np

def cauchy_two_term(x,B,C):
    '''
    This is a simple two term cauchy for index of refraction as function of wavelength.
    https://en.wikipedia.org/wiki/Cauchy%27s_equation
    It is used to create a function for refractive index database entries when only tabulated data is available.
    '''
    return B + (C/(x**2))

def glass_index (glass):
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
    #Load glass dictionary
    with open(os.path.join(os.path.dirname(__file__),os.path.join("glass","glass_dict.json")),"r") as f:
        glass_dict = eval(json.dumps(json.load(f))) #Solution to dict keys formatting
    glass_catalog = ["schott", "ohara", "sumita", "cdgm", "hoya", "hikari", "vitron", "ami", "barberini", "lightpath", "lzos", "misc", "nsg"]
    #Process input and return glass file:
    if isinstance(glass, float): #Allow index to still be set by constant value
        constant_index = float(glass)
        return lambda x=0.55: constant_index
    elif isinstance(glass, str): #Input is a glass name
        input_list = glass.split(" ")
        glass_name = input_list[0]
        if len(input_list) > 1:
            glass_manufacturer = input_list[1].lower()
        else:
            glass_manufacturer = None
        if glass_name in glass_dict:
            if glass_manufacturer:
                if glass_manufacturer in glass_dict[glass_name].keys():
                    glass_file = glass_dict[glass_name][glass_manufacturer]
                else:
                    raise Exception("{} is not in the glass catalog for {}.".format(glass_name,glass_manufacturer))
            else:
                glass_file = glass_dict[glass_name][glass_catalog[min([glass_catalog.index(x) for x in glass_dict[glass_name].keys()])]]
        else:
            raise Exception("{} is not in the glass catalog.".format(glass_name))

        #Process glass file and return the correct function for the glass (wavelength in microns)
        split_path = glass_file.split("\\") #Issue when not working on windows
        if len(split_path) == 4:
            glass_file_path = os.path.join(os.path.join(split_path[1],split_path[2]),split_path[3])
        else:
            glass_file_path = os.path.join(os.path.join(os.path.join(split_path[1],split_path[2]),split_path[3]),split_path[4])
        f = yaml.load(open(os.path.join(os.path.dirname(__file__),glass_file_path)), Loader=yaml.BaseLoader)
        type_index_function = f["DATA"][0]["type"]
        if type_index_function == 'formula 1':
            coefficents = f["DATA"][0]['coefficients'].split(" ")
            coefficents = [float(x) for x in coefficents]
            return lambda x=0.55: (1+coefficents[0]+coefficents[1]/(1-(coefficents[2]/x)**2)+coefficents[3]/(1-(coefficents[4]/x)**2))**.5
        if type_index_function == 'formula 2':
            coefficents = f["DATA"][0]['coefficients'].split(" ")
            coefficents = [float(x) for x in coefficents]
            if len(coefficents) >=7:
                return lambda x=0.55: (1+coefficents[1]/(1-coefficents[2]/x**2)+coefficents[3]/(1-coefficents[4]/x**2)+coefficents[5]/(1-coefficents[6]/x**2))**.5
            else:
                return lambda x=0.55: (1+coefficents[1]/(1-coefficents[2]/x**2)+coefficents[3]/(1-coefficents[4]/x**2))**.5
        if type_index_function == 'formula 3':
            coefficents = f["DATA"][0]['coefficients'].split(" ")
            coefficents = [float(x) for x in coefficents]
            return lambda x=0.55: (coefficents[0]-coefficents[1]*x**coefficents[2]+coefficents[3]*x**coefficents[4]+coefficents[5]*x**coefficents[6]+coefficents[7]*x**coefficents[8]+coefficents[9]*x**coefficents[10])**.5
        if type_index_function == 'formula 5':
            coefficents = f["DATA"][0]['coefficients'].split(" ")
            coefficents = [float(x) for x in coefficents]
            return lambda x=0.55: coefficents[0]-coefficents[1]*x**coefficents[2]+coefficents[3]*x**-coefficents[4]
        if type_index_function == 'tabulated n':
            raw_data = f["DATA"][0]['data'].split('\n')[:-1]
            wavelength = np.array([float(x.split(' ')[0]) for x in raw_data])
            index = np.array([float(x.split(' ')[1]) for x in raw_data])
            popt, _ = curve_fit(cauchy_two_term, wavelength, index) #Curve fit cauchy two term function
            return lambda x=0.55: popt[0] + (popt[1]/(x**2))
        if type_index_function == 'tabulated nk':
            raw_data = f["DATA"][0]['data'].split('\n')[:-1]
            wavelength = np.array([float(x.split(' ')[0]) for x in raw_data])
            index = np.array([float(x.split(' ')[1]) for x in raw_data])
            popt, _ = curve_fit(cauchy_two_term, wavelength, index) #Curve fit cauchy two term function
            return lambda x=0.55: popt[0] + (popt[1]/(x**2)) #Many of these materials are IR glass so the default0.55 is way out of the tabulated range...