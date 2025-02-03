import numpy as np

from .geometry import geometry

from typing import Optional, Tuple, Union

def transform_points(R: np.ndarray,
                     surface: geometry,
                     points: np.ndarray) -> np.ndarray:
    """Transforms points into the reference frame of surface.

    Note
    ----
    Arrays are flattened before return if they only have one row.

    Parameters
    ----------
    R : np.array((3,3))
        Rotation matrix for surface.
    surface : geometry object
        Surface whose reference frame to transform into.
    points : 2d np.array
        Points for each row that will be transformed.

    Returns
    -------
    tran_points : 2d np.array
        Points in the transformed reference frame.
    """
    points = np.atleast_2d(points)
    tran_points = np.dot(R, (points - surface.P).T).T
    if tran_points.shape[0] == 1:
        tran_points = tran_points.flatten()
    return tran_points

def transform_dir(R: np.ndarray,
                  surface: geometry,
                  D: np.ndarray) -> np.ndarray:
    """Transforms directions into the reference frame of surface.

    Note
    ----
    Arrays are flattened before return if they only have one row.

    Parameters
    ----------
    R : np.array((3,3))
        Rotation matrix for surface.
    surface : geometry object
        Surface whose reference frame to transform into.
    D : 2d np.array
        Directions for each row that will be transformed.

    Returns
    -------
    tran_D : 2d np.array
        Directions in the transformed reference frame.
    """
    D = np.atleast_2d(D)
    tran_D = np.dot(R, D.T).T
    if tran_D.shape[0] == 1:
        tran_D = tran_D.flatten()
    return tran_D

def lab_frame_points(R: np.ndarray,
                     surface: geometry,
                     points: np.ndarray) -> np.ndarray:
    """Transforms points into the lab frame.

    Note
    ----
    Arrays are flattened before return if they only have one row.

    Parameters
    ----------
    R : np.array((3,3))
        Rotation matrix for surface.
    surface : geometry object
        Surface whose reference frame to transform from.
    points : 2d np.array
        Points for each row that will be transformed.

    Returns
    -------
    lab_points : 2d np.array
        Points in the lab reference frame.
    """
    points = np.atleast_2d(points)
    lab_points = np.dot(R.T, points.T).T + surface.P
    if lab_points.shape[0] == 1:
        lab_points = lab_points.flatten()
    return lab_points

def lab_frame_dir(R: np.ndarray,
                  surface: geometry,
                  D: np.ndarray) -> np.ndarray:
    """Transforms directions into the lab frame.

    Note
    ----
    Arrays are flattened before return if they only have one row.

    Parameters
    ----------
    R : np.array((3,3))
        Rotation matrix for surface.
    surface : geometry object
        Surface whose reference frame to transform from.
    D : 2d np.array
        Directions for each row that will be transformed.

    Returns
    -------
    lab_D : 2d np.array
        Directions in the lab reference frame.
    """
    D = np.atleast_2d(D)
    lab_D = np.dot(R.T, D.T).T
    if lab_D.shape[0] == 1:
        lab_D = lab_D.flatten()
    return lab_D