import numpy as np
from pytest import approx

from tracepy import gen_rot, transform, lab_frame

def test_rot():
    """ Simple rotation matrix cases. """
    assert np.all(gen_rot([0., 0., 0.]) == np.identity(3))