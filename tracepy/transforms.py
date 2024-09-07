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
        Surface whos reference frame to transform into.
    points : 2d np.array
        Point for each row that will be transformed.

    Returns
    -------
    tran_points : 2d np.array
        Points in the transformed reference frame.

    """

    tran_points = np.array(np.dot(R, (points-surface.P).T).T)
    if len(tran_points) == 1:
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
        Surface whos reference frame to transform into.
    D : 2d np.array
        Direction for each row that will be transformed.

    Returns
    -------
    tran_D: 2d np.array
        Directions in the transformed reference frame.

    """

    tran_D = np.array(np.dot(R, D.T).T)
    if len(tran_D) == 1:
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
        Surface whos reference frame to transform from.
    points : 2d np.array
        Point for each row that will be transformed.

    Returns
    -------
    lab_points : 2d np.array
        Points in the lab reference frame.

    """

    lab_points = np.array(np.dot(R.T, points.T).T)+surface.P
    if len(lab_points) == 1:
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
        Surface whos reference frame to transform from.
    D : 2d np.array
        Direction for each row that will be transformed.

    Returns
    -------
    lab_D : 2d np.array
        Directions in the lab reference frame.

    """

    lab_D = np.array(np.dot(R.T, D.T).T)
    if len(lab_D) == 1:
        lab_D = lab_D.flatten()
    return lab_D
