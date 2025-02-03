"""
Ray tracing and optical design module for Python
================================================

TracePy is a sequential ray tracing package written
in Python for designing optical systems in the geometric
optics regime. It features lens optimization from Scipy.
TracePy is currently in active development and any collaborators
would be welcome.
"""

from .ray import RayGroup
from .geometry import geometry
from .optimize import optimize
from .geoplot import plotxz, plotyz, plot2d
from .optplot import spotdiagram, plotobject, rayaberration, spot_rms
from .iotables import save_optics
from .raygroup import ray_plane
from .index import cauchy_two_term, glass_index
from .utils import gen_rot