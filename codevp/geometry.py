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
        self.kappa = params.get('kappa', 0.) #Specifies type of conic
        self.r = params.get('r', 0.) #Outer radius
        self.r2 = params.get('r2', 0.) #Inner radius
        self.c = params.get('c', 0.) #Vertex curvature
        self.name = params.get('name', None) #Name of surface (optional)
        self.R = gen_rot(self.D) #Rotation matrix for surface
        if self.shape!='plane' and self.kappa>0:
            self.c = sqrt(1/(self.kappa*pow(self.r,2)))
        elif self.shape!='plane' and self.kappa <= 0 and self.c == 0:
            raise Exception(" Specify a vertex curvature c.")
        elif self.shape=='plane' and self.r==0.:
            raise Exception("Specify a radius for your plane.")
        elif self.shape=='rot_sym' and self.r==0:
            raise Exception("Specify radius for rotational surface.")
        
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
        if rho > self.r or rho < self.r2: #Not on plane, success=False.
            return 0., derivitive, False
        return function, derivitive, True
    
    def conics(self, point):
        """ Returns function and derivitive for conics and sphere surfaces. """
        X,Y,Z = point
        rho = sqrt(pow(X,2) + pow(Y, 2))
        failure = self.kappa*pow(self.c, 2)*pow(rho,2) > 1 or (self.kappa <= 0 and rho > self.r)
        if failure: #Not on surface, success=False.
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
        function[np.array(rho > self.r) + np.array(rho < self.r2)] = np.nan #All points outside of plane are nan.
        return function      
    
    def conics_plot(self, point):
        """ Returns Z value for an array of points for plotting conics. """
        X, Y = point[:,0], point[:,1]
        function = np.zeros(len(point))
        if self.kappa > 0:
            nan_idx = self.kappa*pow(self.c, 2)*(pow(X,2)+pow(Y, 2)) > 1 #Points not on surface.
        else:
            nan_idx = sqrt(pow(X,2)+pow(Y, 2)) > self.r #For kappa<=0 surfaces
        rho = sqrt(pow(X[~nan_idx],2) + pow(Y[~nan_idx], 2))
        function[nan_idx] = np.nan
        function[~nan_idx] = self.c*pow(rho, 2)/(1 + pow((1-self.kappa*pow(self.c, 2)*pow(rho,2)), 0.5))
        return function 
