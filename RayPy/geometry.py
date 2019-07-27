import numpy as np
from numpy import cos, sin, pi, sqrt

from .transforms import *

class geometry:
    """ Class for the different surfaces in a optical system. """
    def __init__(self, params):
        self.P = params['P'] #Position of surface
        self.D = np.array(params.get('D', [0., 0., 0.])) #Rotation angles
        self.action = params['action'] #Interaction with surface
        self.Diam = params['Diam'] #Outer diameter
        self.N = params.get('N', 1.) #Index of refraction     
        self.kappa = params.get('kappa', None) #Specifies type of conic
        self.diam = params.get('diam', 0.) #Inner diameter
        self.c = params.get('c', 0.) #Vertex curvature. Default to plane.
        self.name = params.get('name', None) #Name of surface (optional)
        self.R = gen_rot(self.D) #Rotation matrix for surface
        self.check_params()

    def __getitem__(self, item):
        """ Return attribute of geometry. """
        return getattr(self, item)

    def __setitem__(self, item, value):
        """ Set attribute of geometry. """
        return setattr(self, item, value)

    def check_params(self):
        """ Check that required parameters are given and update needed parameters. """
        if isinstance(self.P, float) or isinstance(self.P, int): #Allow on axis integer for P.
            self.P = np.array([0., 0., self.P])
        else:
            self.P = np.array(self.P)
        if self.c != 0:
            if self.kappa is None:
                raise Exception("Specify a kappa for this conic.")
            elif self.kappa>0:
                print("Warning: Specified c value is not used when kappa>0")  
                self.c = sqrt(1/(self.kappa*pow(self.Diam/2.,2))) 
        elif self.c == 0 and self.kappa is None:
            self.kappa = 1. #Used for planes, does not affect calculations.       

    def get_surface(self, point):
        """ Returns the function and derivitive of a surface for a point. """
        return self.conics(point)

    def get_surface_plot(self, points):
        """ Returns the function value for an array of points. """
        return self.conics_plot(points)
    
    def conics(self, point):
        """ Returns function and derivitive for conics and sphere surfaces. """
        X,Y,Z = point
        rho = sqrt(pow(X,2) + pow(Y, 2))
        if rho > self.Diam/2. or rho < self.diam/2.: #Not on surface
            raise Exception()
        function = Z - self.c*pow(rho, 2)/(1 + pow((1-self.kappa*pow(self.c, 2)*pow(rho,2)), 0.5)) #Conic equation.
        E = self.c / pow((1-self.kappa*pow(self.c, 2)*pow(rho,2)), 0.5) #See Spencer, Murty section on rotational surfaces.
        derivitive = [-X*E, -Y*E, 1.]
        return function, derivitive

    def conics_plot(self, point):
        """ Returns Z value for an array of points for plotting conics. """
        X, Y = point[:,0], point[:,1]
        rho = sqrt(pow(X,2) + pow(Y,2))
        function = np.zeros(len(point))
        nan_idx = (rho > self.Diam/2.) + (rho < self.diam/2.)
        rho = sqrt(pow(X[~nan_idx],2) + pow(Y[~nan_idx], 2))
        function[nan_idx] = np.nan
        function[~nan_idx] = self.c*pow(rho, 2)/(1 + pow((1-self.kappa*pow(self.c, 2)*pow(rho,2)), 0.5))
        return function 
