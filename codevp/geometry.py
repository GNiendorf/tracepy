import numpy as np
from numpy import cos, sin, pi, sqrt

from .transforms import *

class geometry:
    """ Class for the different surfaces in a optical system. """
    def __init__(self, params):
        self.P = np.array(params['P']) #Position of surface
        self.D = np.array(params['D']) #Rotation angles of surface?
        self.shape = params['shape'] #Shape of surface e.g. conic, plane
        self.inter = params['inter'] #Interaction with surface
        self.N = params.get('N', 1.) #Index of refraction     
        self.kappa = params.get('kappa', None) #Specifies type of conic
        self.Diam = params.get('Diam', None) #Outer radius
        self.diam = params.get('diam', 0.) #Inner radius
        self.c = params.get('c', None) #Vertex curvature
        self.name = params.get('name', None) #Name of surface (optional)
        self.R = gen_rot(self.D) #Rotation matrix for surface
        if self.shape=='conic' and self.kappa>0:
            self.c = sqrt(1/(self.kappa*pow(self.Diam/2.,2)))
        if self.shape=='conic' and self.kappa>0 and self.c is not None:
            print("Warning: c value is not used when kappa>0")
        if self.shape=='conic' and self.c is None:
            raise Exception("Specify a vertex curvature c.")
        if self.shape=='conic' and self.kappa is None:
            raise Exception("Specify a kappa for this conic.")
        if self.shape=='plane' and self.Diam is None:
            raise Exception("Specify a diameter for this plane.")
        if self.shape=='conic' and self.Diam is None:
            raise Exception("Specify a diameter for this conic.")

    def __getitem__(self, item):
        """ Return attribute of geometry. """
        return getattr(self, item)

    def get_surface(self, point, plot=None):
        """ Returns the function and derivitive of a surface. """
        if self.shape == 'plane':
            if plot:
                return self.plane_plot(point)
            return self.plane(point)
        elif self.shape == 'conic':
            if plot:
                return self.conics_plot(point)
            return self.conics(point)
        else:
            raise Exception("The surface chosen is not supported!")  
    
    def plane(self, point):
        """ Returns function and derivitive for plane surfaces. """
        X,Y,Z = point
        rho = sqrt(pow(X,2) + pow(Y,2))
        function = Z
        derivitive = [0., 0., 1.]
        if rho > self.Diam/2. or rho < self.diam/2.: #Not on plane, success=False.
            return 0., derivitive, False
        return function, derivitive, True
    
    def conics(self, point):
        """ Returns function and derivitive for conics and sphere surfaces. """
        X,Y,Z = point
        rho = sqrt(pow(X,2) + pow(Y, 2))
        if rho > self.Diam/2. or rho < self.diam/2.: #Not on surface, success=False.
            return 0., np.ones(3), False
        function = Z - self.c*pow(rho, 2)/(1 + pow((1-self.kappa*pow(self.c, 2)*pow(rho,2)), 0.5)) #Conic equation.
        E = self.c / pow((1-self.kappa*pow(self.c, 2)*pow(rho,2)), 0.5) #See Spencer, Murty section on rotational surfaces.
        derivitive = [-X*E, -Y*E, 1.]
        return function, derivitive, True

    def plane_plot(self, point):
        """ Returns Z=0 value for an array of points for plotting planes. """
        X, Y = point[:,0], point[:,1]
        function = np.zeros(len(point)) #All points in plane obj frame are Z=0.
        rho = sqrt(pow(X,2) + pow(Y,2))
        function[np.array(rho > self.Diam/2.) + np.array(rho < self.diam/2.)] = np.nan #All points outside of plane are nan.
        return function      
    
    def conics_plot(self, point):
        """ Returns Z value for an array of points for plotting conics. """
        X, Y = point[:,0], point[:,1]
        rho = sqrt(pow(X,2) + pow(Y,2))
        function = np.zeros(len(point))
        nan_idx = (rho > self.Diam/2.) + (rho < self.diam/2.) #For kappa<=0 surfaces
        rho = sqrt(pow(X[~nan_idx],2) + pow(Y[~nan_idx], 2))
        function[nan_idx] = np.nan
        function[~nan_idx] = self.c*pow(rho, 2)/(1 + pow((1-self.kappa*pow(self.c, 2)*pow(rho,2)), 0.5))
        return function 
