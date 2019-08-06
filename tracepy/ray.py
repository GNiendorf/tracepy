# Authors: Gavin Niendorf <gavinniendorf@gmail.com>
#
# Classes and methods for defining rays and their propagation rules.
#
# License: MIT

import numpy as np

from .transforms import *
from .exceptions import NormalizationError, NotOnSurfaceError

class ray:
    """Class for rays and their propagation through surfaces.

    Note
    ----
    Also checks whether the direction cosines are normalized.

    Attributes
    ----------
    P : np.array of 3 floats/ints
        Position of ray in the lab frame.
    D : np.array of 3 floats/ints
        Direction cosines for the ray in the lab frame.
    P_hist : list of P np.arrays
        Previous P np.arrays in a list.
    D_hist : list of D np.arrays
        Previous D np.arrays in a list.
    N : float/int
        Index of refraction of current material.
    wvl: float/int
        Wavelength of the ray in microns 550nm --> 0.55.
    """

    def __init__(self, params, N_0=1):
        self.P = np.array(params['P'])
        self.D = np.array(params['D'])
        self.P_hist = [self.P]
        self.D_hist = [self.D]
        self.N = N_0
        self.wvl = params.get('wvl',0.55) #Added default wavelength 550nm
        if abs(np.linalg.norm(self.D)-1.) > .01:
            #Ray direction cosines are not normalized.
            raise NormalizationError()

    def transform(self, surface):
        """ Updates position and direction of a ray to obj coordinate system. """
        self.P, self.D = transform(surface.R, surface, np.array([self.P]), np.array([self.D]))

    def find_intersection(self, surface):
        """Finds the intersection point of a ray with a surface.

        Note
        ----
        Directly changes the self.P (position) attribute of the ray
        that corresponds to the intersection point. Also be aware
        that my error definition is different from Spencer's paper.
        I found that the more direct error equation of abs(F) allows
        me to tune my max error values to get better accuracy.

        Parameters
        ----------
        surface : geometry object
            Surface to find intersection of ray with.

        """

        #Initial guesses, see Spencer, Murty for explanation.
        s_0 = -self.P[2]/self.D[2]
        X_1 = self.P[0]+self.D[0]*s_0
        Y_1 = self.P[1]+self.D[1]*s_0
        s_j = [0., 0.]
        #Initial error.
        error = 1.
        n_iter = 0
        #Max iterations allowed.
        n_max = 1e4
        while error > 1e-6 and n_iter < n_max:
            X, Y, Z = [X_1, Y_1, 0.]+np.dot(self.D, s_j[0])
            try:
                #'normal' is the surface direction numbers.
                func, normal= surface.get_surface([X, Y, Z])
                deriv = np.dot(normal, self.D)
                #Newton-raphson method
                s_j = s_j[1], s_j[1]-func/deriv
            except NotOnSurfaceError:
                self.P = None
                return None
            #Error is how far f(X, Y, Z) is from 0.
            error = abs(func)
            n_iter += 1
        if n_iter == n_max or s_0+s_j[0] < 0 or np.dot(([X, Y, Z]-self.P), self.D) < 0.:
            self.P = None
        else:
            self.normal = normal
            self.P = np.array([X, Y, Z])

    def interact(self, surface, typeof):
        """Updates new direction of a ray for a given interaction type.

        Note
        ----
        High level method that calls the appropriate method for a given
        interaction.

        Parameters
        ----------
        surface : geometry object
            Surface to find intersection of ray with.
        typeof : str
            Type of interaction
            reflection -> Reflect the ray off the surface.
            refraction -> Refract the ray into the surface.
            stop -> Don't change ray direction.

        """
        if hasattr(surface,'glass'):
            mu = self.N / surface.glass(self.wvl)
        else:
            mu = self.N / surface.N
        a = mu*np.dot(self.D, self.normal)/pow(np.linalg.norm(self.normal), 2)
        b = (pow(mu,2)-1)/pow(np.linalg.norm(self.normal), 2)
        if typeof == 'stop':
            pass
        #Needed for total internal reflection even if typeof is refraction.
        elif b > pow(a, 2) or typeof == 'reflection':
            self.reflection(surface, a/mu)
        elif typeof == 'refraction':
            self.refraction(surface, mu, a, b)

    def reflection(self, surface, a):
        """Reflects the ray off a surface and updates the ray's direction.

        Note
        ----
        This method computes D exactly rather than numerically like in the
        refraction method.

        Parameters
        ----------
        surface : geometry object
            Surface to reflect from.
        a : float/int
            Constant defined in the interact method.

        """

        k, l, m = self.D
        K, L, M = self.normal
        self.D = np.array([k-2.*a*K, l-2.*a*L, m-2.*a*M])

    def refraction(self, surface, mu, a, b):
        """Simulates refraction of a ray into a surface and updates the ray's direction.

        Note
        ----
        My error definition is not in Spencer and Murty's paper but is inspired by my
        unique intersection error definition. We are solving for roots of a quadratic and
        I am defining my error by how far the quadtratic is from 0. See Spencer, Murty for
        derivation of the quadratic.

        Parameters
        ----------
        surface : geometry object
            Surface to refract into.
        mu, a, b : float/int
            Constants defined in the interact method.

        Returns
        -------
        0
            Returns 0 if the number of iterations exceeds the max allowed to converge.

        """

        k, l, m = self.D
        K, L, M = self.normal
        G = [-b/(2*a), -b/(2*a)]
        #Initial error.
        error = 1.
        niter = 0
        #Max iterations allowed.
        nmax = 1e5
        while error > 1e-15 and niter < nmax:
            #Newton-raphson method
            G = G[1], (pow(G[1],2)-b)/(2*(G[1]+a))
            #See Spencer, Murty for where this is inspired by.
            error = abs(pow(G[1],2)+2*a*G[1]+b)
            niter += 1
        if niter==nmax:
            self.P = None
            return 0.
        #Update direction and index of refraction of the current material.
        self.D = np.array([mu*k+G[1]*K,mu*l+G[1]*L,mu*m+G[1]*M])
        if hasattr(surface,'glass'):
            self.N = surface.glass(self.wvl)
        else:
            self.N = surface.N

    def ray_lab_frame(self, surface):
        """ Updates position and direction of a ray in the lab frame. """
        self.P, self.D = lab_frame(surface.R, surface, np.array([self.P]), np.array([self.D]))

    def update(self):
        """ Updates the P_hist and D_hist arrays from current P and D arrays. """
        self.P_hist.append(self.P)
        self.D_hist.append(self.D)

    def propagate(self, surfaces):
        """Propagates a ray through a given surfaces list.

        Note
        ----
        If self.P is None then the ray failed to converge or
        took too many iterations to meet the required accuracy.
        Note that this is used (self.P is None) as a flag in
        many other functions in TracePy.

        Parameters
        ----------
        surfaces : list of geometry objects
            Surfaces to propagate through in order of propagation.

        """

        for surface in surfaces:
            self.transform(surface)
            self.find_intersection(surface)
            #Results from failure to converge.
            if self.P is None:
                break
            self.interact(surface, surface.action)
            #Results from too many iterations.
            if self.P is None:
                break
            self.ray_lab_frame(surface)
            #Update current to history arrays.
            self.update()
