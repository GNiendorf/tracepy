import numpy as np
from numpy import cos, sin, pi, sqrt

def gen_rot(ang):
    """ Returns rotation matrix from 3 rotation angles. """
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
    """ 
        Inputs of points and D must both be 2d arrays 
        Returns both points and D in transformed frame.
        Arrays are flattened if they have one row. 
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
    """ 
        Inputs of points and D must both be 2d arrays 
        Returns both points and D in lab frame.
        Arrays are flattened if they have one row. 
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
