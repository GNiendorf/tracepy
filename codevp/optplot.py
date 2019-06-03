import matplotlib.pyplot as plt
import numpy as np
from numpy import cos, sin, pi, sqrt

from .geometry import geometry
from .transforms import transform

class optplot:
    """ Plots useful features of an optical system such as spot diagrams. """
    def __init__(self, geo_params, rays=None):
        self.rays = rays
        self.geo_params = geo_params
        self.geo = [geometry(surf) for surf in self.geo_params]
        
    def spotdiagram(self, pltparams = {'color': 'red'}, return_ = None):
        """ Plots the transformed intersection points of rays on the stop surface. """
       # if self.rays is None: 
       #     self.rays = ray_group('plane', self.geo, inters, [0., 0., 0.], 1.8, d=[0.,0.,1.], nrays=100)
        points = np.array([rayiter.P for rayiter in self.rays if rayiter.P is not None])
        if points.size == 0: #No rays propagated through entire system.
            raise Exception("No rays survived.")
        stop = self.geo[-1]
        points_obj = transform(stop.R, stop, points) # Get X,Y points in obj. reference frame. - 
                                                     # Used for stops orientated in arbitrary direction.
        points_obj = np.around(points_obj, 14) #Round arrays to upper bound on accuracy.            
        X, Y = points_obj[:,0], points_obj[:,1]      
        rms = np.std(points_obj[:,[0,1]] - points_obj[:,[0,1]].mean(axis=0))
        plt.subplot(1,1,1, aspect='equal')
        plt.locator_params(axis='x', nbins=8)
        plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
        plt.plot(X, Y, '+', **pltparams)
        plt.text(0.65, 1.015, 'RMS=%.2E' % rms, transform=plt.gca().transAxes)
        if return_ is not None: #Return points if asked.
            return X, Y
        
    def rayaberation(self, return_ = None):
        """ Plots the ray aberation curves for an optical system. """
        pass
