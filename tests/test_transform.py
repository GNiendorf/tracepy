import numpy as np

from tracepy import gen_rot

def test_rot():
    """ Simple rotation matrix cases. """
    assert np.all(gen_rot([0., 0., 0.]) == np.identity(3))