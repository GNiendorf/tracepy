"""
Ray tracing and optical design module for Python
================================================

TracePy is a sequential ray tracing package written
in Python for designing optical systems in the geometric
optics regime. It features lens optimization from Scipy.
TracePy is currently in active development and any collaborators
would be welcome.

"""

from .ray import ray
from .geometry import geometry
from .optimize import optimize
from .geoplot import plotxz, plotyz, plot2d
from .optplot import spotdiagram, plotobject, rayaberration
from .iotables import *
from .transforms import *
from .raygroup import *
from.index import *