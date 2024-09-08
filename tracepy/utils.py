import numpy as np

def gen_rot(ang: np.ndarray) -> np.ndarray:
    """Returns a rotation matrix from 3 rotation angles.

    Parameters
    ----------
    ang : np.array of length 3
        Euler angles alpha, beta, gamma in the lab frame.

    Returns
    -------
    np.array((3,3))
        Returns the rotation matrix.

    """

    alpha, beta, gamma = ang
    R_11 = np.cos(alpha)*np.cos(gamma)+np.sin(alpha)*np.sin(beta)*np.sin(gamma)
    R_12 = -np.cos(beta)*np.sin(gamma)
    R_13 = -np.sin(alpha)*np.cos(gamma)+np.cos(alpha)*np.sin(beta)*np.sin(gamma)
    R_21 = np.cos(alpha)*np.sin(gamma)-np.sin(alpha)*np.sin(beta)*np.cos(gamma)
    R_22 = np.cos(beta)*np.cos(gamma)
    R_23 = -np.sin(alpha)*np.sin(gamma)-np.cos(alpha)*np.sin(beta)*np.cos(gamma)
    R_31 = np.sin(alpha)*np.cos(beta)
    R_32 = np.sin(beta)
    R_33 = np.cos(alpha)*np.cos(beta)
    R = np.array([[R_11, R_12, R_13],\
                   [R_21, R_22, R_23],\
                   [R_31, R_32, R_33]])
    return R
