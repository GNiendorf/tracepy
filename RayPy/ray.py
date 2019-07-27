import numpy as np
from numpy import cos, sin, pi, sqrt

from .transforms import *

class ray:
    """ Class for rays and their propagation through surfaces. """
    def __init__(self, params, N_0=1):
        self.P = np.array(params['P'])
        self.D = np.array(params['D'])
        self.P_hist = [self.P]
        self.D_hist = [self.D]
        self.N = N_0
        if abs(np.linalg.norm(self.D)-1.) > .01:
            raise Exception("Ray direction cosines are not normalized.")

    def transform(self, surface):
        """ Updates position and direction of a ray to obj coordinate system. """
        self.P, self.D = transform(surface.R, surface, np.array([self.P]), np.array([self.D]))

    def find_intersection(self, surface):
        """ Finds the intersection point of a ray with a surface. """
        s_0 = -self.P[2]/self.D[2] #Initial guesses, see Spencer, Murty.
        X_1 = self.P[0]+self.D[0]*s_0
        Y_1 = self.P[1]+self.D[1]*s_0
        s_j = [0., 0.]
        error = 1. #Initial error.
        n_iter = 0
        n_max = 1e4 #Max iterations allowed.
        while error > 1e-6 and n_iter < n_max:
            X, Y, Z = [X_1, Y_1, 0.]+np.dot(self.D, s_j[0])
            try:
                func, normal= surface.get_surface([X, Y, Z]) #'normal' is the surface direction numbers.
                deriv = np.dot(normal, self.D)
                s_j = s_j[1], s_j[1]-func/deriv #Newton-raphson method
            except:
                self.P = None
                return None
            error = abs(func) #Error is how far f(X, Y, Z) is from 0.
            n_iter += 1 
        if n_iter == n_max or s_0+s_j[0] < 0 or np.dot(([X, Y, Z]-self.P), self.D) < 0.:
            self.P = None
        else:
            self.normal = normal
            self.P = np.array([X, Y, Z])

    def interact(self, surface, typeof):
        """ Updates new direction of a ray for an interaction type. """
        mu = self.N/surface.N
        a = mu*np.dot(self.D, self.normal)/pow(np.linalg.norm(self.normal), 2)
        b = (pow(mu,2)-1)/pow(np.linalg.norm(self.normal), 2)
        if typeof == 'stop':
            pass
        elif b > pow(a, 2) or typeof == 'reflection': 
            self.reflection(surface, a/mu)
        elif typeof == 'refraction':
            self.refraction(surface, mu, a, b)
            
    def reflection(self, surface, a):
        """ Reflects the ray off a surface and updates the ray's direction. """
        k, l, m = self.D
        K, L, M = self.normal
        self.D = np.array([k-2.*a*K, l-2.*a*L, m-2.*a*M])        
        
    def refraction(self, surface, mu, a, b):
        """ Simulates refraction of a ray into a surface and updates the ray's direction. """
        k, l, m = self.D
        K, L, M = self.normal        
        G = [-b/(2*a), -b/(2*a)]
        error = 1. #Initial error.
        niter = 0 
        nmax = 1e5 #Max iterations allowed.
        while error > 1e-15 and niter < nmax: 
            G = G[1], (pow(G[1],2)-b)/(2*(G[1]+a)) #Newton-raphson method
            error = abs(pow(G[1],2)+2*a*G[1]+b)
            niter += 1
        if niter==nmax:
            self.P = None
            return 0.
        #Update direction and index of refraction of the current material.
        self.D = np.array([mu*k+G[1]*K,mu*l+G[1]*L,mu*m+G[1]*M])
        self.N = surface.N
    
    def ray_lab_frame(self, surface):
        """ Updates position and direction of a ray in the lab frame. """
        self.P, self.D = lab_frame(surface.R, surface, np.array([self.P]), np.array([self.D]))
    
    def update(self):
        """ Updates the P_hist and D_hist arrays from current P and D arrays. """
        self.P_hist.append(self.P)
        self.D_hist.append(self.D)
        
    def propagate(self, surfaces):
        """ Propagates a ray through a given surface with a given interaction. """
        for surface in surfaces:
            self.transform(surface)
            self.find_intersection(surface)
            if self.P is None: #Results from failure to converge.
                break
            self.interact(surface, surface.action)
            if self.P is None: #Results from too many iterations.
                break
            self.ray_lab_frame(surface)
            self.update() #Update current to history arrays.
