class NormalizationError(Exception):
    """ Custom exception for unnormalized input. """

class NotOnSurfaceError(Exception):
    """ Error for rays that do not intersect with a surface. """

class TraceError(Exception):
    """ Custom error for lens systems where no rays survive being traced. """

class InvalidGeometry(Exception):
    """ Invalid parameters were given to define a geometry object. """
