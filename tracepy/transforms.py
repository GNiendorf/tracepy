# Authors: Gavin Niendorf <gavinniendorf@gmail.com>
#
# Functions for transforming between reference frames.
#
# License: MIT

import numpy as np
from numpy import cos, sin

def gen_rot(ang):
    """Returns a rotation matrix from 3 rotation angles.

    Note
    ----
    np.matrix is deprecated and needs to be changed to a
    normal np.array but I remember difficulties with the
    regular np.array. Needs to be fixed eventually.

    Parameters
    ----------
    ang : np.array of length 3
        Euler angles alpha, beta, gamma in the lab frame.

    Returns
    -------
    np.matrix((3,3))
        Returns the rotation matrix.

    """

    alpha, beta, gamma = ang
    R_11 = cos(alpha)*cos(gamma)+sin(alpha)*sin(beta)*sin(gamma)
    R_12 = -cos(beta)*sin(gamma)
    R_13 = -sin(alpha)*cos(gamma)+cos(alpha)*sin(beta)*sin(gamma)
    R_21 = cos(alpha)*sin(gamma)-sin(alpha)*sin(beta)*cos(gamma)
    R_22 = cos(beta)*cos(gamma)
    R_23 = -sin(alpha)*sin(gamma)-cos(alpha)*sin(beta)*cos(gamma)
    R_31 = sin(alpha)*cos(beta)
    R_32 = sin(beta)
    R_33 = cos(alpha)*cos(beta)
    R = np.matrix([[R_11, R_12, R_13],\
                   [R_21, R_22, R_23],\
                   [R_31, R_32, R_33]])
    return R

def transform(R, surface, points, D=None):
    """Transforms points into the reference frame of surface.

    Note
    ----
    Arrays are flattened before return if they only have one row.

    Parameters
    ----------
    R : np.matrix((3,3))
        Rotation matrix for surface.
    surface : geometry object
        Surface whos reference frame to transform into.
    points : 2d np.array
        Point for each row that will be transformed.
    D (optional): 2d np.array
        Direction for each row that will be transformed.

    Returns
    -------
    tran_points : 2d np.array
        Points in the transformed reference frame.
    tran_D (optional): 2d np.array
        Directions in the transformed reference frame.

    """

    tran_points = np.array(np.dot(R, (points-surface.P).T).T)
    if D is not None:
        tran_D = np.array(np.dot(R, D.T).T)
        if len(tran_points) == 1:
            tran_points = tran_points.flatten()
            tran_D = tran_D.flatten()
        return tran_points, tran_D
    if len(tran_points) == 1:
        tran_points = tran_points.flatten()
    return tran_points

def lab_frame(R, surface, points, D=None):
    """Transforms points into the lab frame.

    Note
    ----
    Arrays are flattened before return if they only have one row.

    Parameters
    ----------
    R : np.matrix((3,3))
        Rotation matrix for surface.
    surface : geometry object
        Surface whos reference frame to transform from.
    points : 2d np.array
        Point for each row that will be transformed.
    D (optional): 2d np.array
        Direction for each row that will be transformed.

    Returns
    -------
    lab_points : 2d np.array
        Points in the lab reference frame.
    lab_D (optional): 2d np.array
        Directions in the lab reference frame.

    """

    lab_points = np.array(np.dot(R.T, points.T).T)+surface.P
    if D is not None:
        lab_D = np.array(np.dot(R.T, D.T).T)
        if len(lab_points) == 1:
            lab_points = lab_points.flatten()
            lab_D = lab_D.flatten()
        return lab_points, lab_D
    if len(lab_points) == 1:
        lab_points = lab_points.flatten()
    return lab_points
