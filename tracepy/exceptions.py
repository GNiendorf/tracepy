# Authors: Gavin Niendorf <gavinniendorf@gmail.com>
#
# Custom warnings and error classes used in TracePy.
#
# License: MIT

class NormalizationError(Exception):
    """ Custom exception for unnormalized input. """

class NotOnSurfaceError(Exception):
    """ Error for rays that do not intersect with a surface. """

class TraceError(Exception):
    """ Custom error for lens systems where no rays survive being traced. """
