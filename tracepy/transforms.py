import numpy as np

from .geometry import geometry

from typing import Optional, Tuple, Union

def transform(R: np.ndarray,
              surface: geometry,
              points: np.ndarray,
              D: Optional[np.ndarray] = None) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
    """Transforms points into the reference frame of surface.

    Note
    ----
    Arrays are flattened before return if they only have one row.

    Parameters
    ----------
    R : np.array((3,3))
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

def lab_frame(R: np.ndarray,
              surface: geometry,
              points: np.ndarray,
              D: Optional[np.ndarray] = None) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
    """Transforms points into the lab frame.

    Note
    ----
    Arrays are flattened before return if they only have one row.

    Parameters
    ----------
    R : np.array((3,3))
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
