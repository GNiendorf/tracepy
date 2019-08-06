# Authors: Gavin Niendorf <gavinniendorf@gmail.com>
#
# Classes and methods for defining surfaces.
#
# License: MIT

import numpy as np

from .transforms import *
from .exceptions import NotOnSurfaceError
from .index import *

class geometry:
    """Class for the different surfaces in an optical system.

    Attributes
    ----------
    P : list of 3 floats/ints or just float/int
        Position of surface in the lab frame. If only one float or int is given
        then it is considered the position in the Z direction.
    D (optional): list of 3 floats/ints
        Rotation angles alpha, beta, gamma for a surface.
        Defaults to no rotation -> [0., 0., 0.]
    action : str
        Interaction of rays with the specified surface. Either reflection,
        refraction, or stop (do nothing).
    Diam : float/int
        Outer diameter of the specified surface.
    diam : float/int
        Inner diameter of the specified surface.
    N (optional): float/int
        Index of refraction of the material that the rays propagate INTO
        (meaning on the image side). If N is not given then it defaults
        to a value of 1 for air. A constant value can be specified such as 1.72.
    kappa (optional): float/int
        Defines the type of conic that the surface models.
        If kappa is None then the surface is planar.
        kappa < 0 -> hyperboloid
        kappa = 0 -> paraboloid
        0 < kappa < 1 -> hemelipsoid of revolution about major axis
        kappa = 1 -> hemisphere
        kappa > 1 -> hemelipsoid of revolution about minor axis
    c (optional): float/int
        Vertex curvature of the surface.
        If c is 0 then the surface is planar.
    name (optional): str
        Name of the surface, used for optimization
    R (generated): np.matrix((3,3))
        Rotation matrix for the surface from rotation angles D.

    """

    def __init__(self, params):
        self.P = params['P']
        self.D = np.array(params.get('D', [0., 0., 0.]))
        self.action = params['action']
        self.Diam = params['Diam']
        self.N = params.get('N', 1.)
        self.kappa = params.get('kappa', None)
        self.diam = params.get('diam', 0.)
        self.c = params.get('c', 0.)
        self.name = params.get('name', None)
        self.R = gen_rot(self.D)
        if params.get('glass'):
            self.glass = glass_index(params.get('glass'))
        self.check_params()

    def __getitem__(self, item):
        """ Return attribute of geometry. """
        return getattr(self, item)

    def __setitem__(self, item, value):
        """ Set attribute of geometry. """
        return setattr(self, item, value)

    def check_params(self):
        """Check that required parameters are given and update needed parameters.

        Summary
        -------
            If P is given as a float/int then it is converted to a np array
            with that float/int in the Z direction. If c != 0 (in the case
            of a conic) then kappa must be specified, and if kappa is greater
            than 0 then the value of c is redundant by boundary conditions of
            the conic equation. Lastly, if c == 0 in the case of a planar
            surface the None value of kappa needs to be set to a dummy value
            to avoid exceptions in calculating the conic equation. Note that
            this does not affect the calculation since c is 0.

        """

        if isinstance(self.P, float) or isinstance(self.P, int):
            #Allow on axis integer for P.
            self.P = np.array([0., 0., self.P])
        else:
            self.P = np.array(self.P)
        if self.c != 0:
            if self.kappa is None:
                raise Exception("Specify a kappa for this conic.")
            elif self.kappa>0:
                print("Warning: Specified c value is not used when kappa>0")
                self.c = np.sqrt(1/(self.kappa*pow(self.Diam/2.,2)))
        elif self.c == 0 and self.kappa is None:
            #Used for planes, does not affect calculations.
            self.kappa = 1.

    def get_surface(self, point):
        """ Returns the function and derivitive of a surface for a point. """
        return self.conics(point)

    def get_surface_plot(self, points):
        """ Returns the function value for an array of points. """
        return self.conics_plot(points)

    def conics(self, point):
        """Returns function value and derivitive list for conics and sphere surfaces.

        Note
        ----
        This differs from the conics_plot method because we are returning the function
        value from the conic at the point X,Y,Z rather than the Z value for a given X, Y
        pair. A surface is defined by the equation F(X, Y, Z) = 0, so if the ray
        intersects with the surface this function will return a function value of
        0.

        Parameters
        ----------
        point : np.array of length 3
            X, Y, and Z values for the point.

        Returns
        -------
        function : float
            Function value of the surface. A value of 0 corresponds to the ray intersecting
            the surface.
        derivitive : list of length 3
            Function derivitive of the surface at the point given. Used to determine which
            direction the ray needs to travel in and the step size to intersect the surface.

        """

        X,Y,Z = point
        rho = np.sqrt(pow(X,2) + pow(Y, 2))
        if rho > self.Diam/2. or rho < self.diam/2.:
            raise NotOnSurfaceError()
        #Conic equation.
        function = Z - self.c*pow(rho, 2)/(1 + pow((1-self.kappa*pow(self.c, 2)*pow(rho,2)), 0.5))
        #See Spencer, Murty section on rotational surfaces for definition of E.
        E = self.c / pow((1-self.kappa*pow(self.c, 2)*pow(rho,2)), 0.5)
        derivitive = [-X*E, -Y*E, 1.]
        return function, derivitive

    def conics_plot(self, point):
        """Returns Z values for an array of points for plotting conics.

        Parameters
        ----------
        point : 2d np.array with shape (n, 3)
            X, Y, and Z values for each point (for each row).

        Returns
        -------
        function : 1d np.array
            np.array of Z values for each X,Y pair.

        """

        X, Y = point[:,0], point[:,1]
        rho = np.sqrt(pow(X,2) + pow(Y,2))
        #Initialize Z value array
        function = np.zeros(len(point))
        nan_idx = (rho > self.Diam/2.) + (rho < self.diam/2.)
        rho = np.sqrt(pow(X[~nan_idx],2) + pow(Y[~nan_idx], 2))
        function[nan_idx] = np.nan
        function[~nan_idx] = self.c*pow(rho, 2)/(1 + pow((1-self.kappa*pow(self.c, 2)*pow(rho,2)), 0.5))
        return function