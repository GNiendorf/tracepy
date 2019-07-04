import matplotlib.pyplot as plt
import numpy as np
from numpy import cos, sin, pi, sqrt

from .geometry import geometry
from .transforms import transform
from .raygroup import ray_plane

class optplot:
    """ Plots useful features of an optical system such as spot diagrams. """
    def __init__(self, geo_params, rays=None):
        self.rays = rays
        self.geo_params = geo_params
        self.geo = [geometry(surf) for surf in self.geo_params]
        if self.rays is None: #Needs to be fixed for arbitrary initial direction.
            self.rays = ray_plane(self.geo_params, [0., 0., 0.], self.geo[0]["Diam"]/2., d=[0.,0.,1.], nrays=5000)
        self.startpoints = np.array([rayiter.P_hist[0] for rayiter in self.rays])
        self.endpoints = np.array([rayiter.P_hist[-1] for rayiter in self.rays if rayiter.P is not None])
        #Check if any rays survived.
        if self.endpoints.size == 0: #No rays propagated through entire system.
            raise Exception("No rays survived.")
        
    def spotdiagram(self, pltparams = {'color': 'red'}, optimizer=False):
        """ Plots the transformed intersection points of rays on the stop surface. """
        stop = self.geo[-1]
        points_obj = transform(stop.R, stop, self.endpoints) # Get X,Y points in obj. reference frame.
        points_obj = np.around(points_obj, 14) #Round arrays to upper bound on accuracy.            
        X, Y = points_obj[:,0], points_obj[:,1]    
        rms = np.std(points_obj[:,[0,1]] - points_obj[:,[0,1]].mean(axis=0))
        if optimizer:  
            return rms
        plt.subplot(1,1,1, aspect='equal')
        plt.locator_params(axis='x', nbins=8)
        plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
        plt.plot(X, Y, '+', **pltparams)
        plt.text(0.65, 1.015, 'RMS=%.2E' % rms, transform=plt.gca().transAxes)
        
    def plotobject(self, pltparams = {'color': 'blue'}): 
        """ Plots the initial ray locations in the initial objects' frame. """
        start = self.geo[0]
        points_obj = transform(start.R, start, self.startpoints) # Get X,Y points in obj. reference frame.
        X, Y = points_obj[:,0], points_obj[:,1]
        plt.subplot(1,1,1, aspect='equal') 
        plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
        rms = np.std(self.startpoints[:,[0,1]] - self.startpoints[:,[0,1]].mean(axis=0)) 
        plt.text(0.65, 1.015, 'RMS=%.2E' % rms, transform=plt.gca().transAxes) 
        plt.plot(X, Y, '+', **pltparams)
    
    def rayaberration(self):
        """ Plots the ray aberration diagrams for an optical system. """
        pass 
