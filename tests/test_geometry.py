import numpy as np
from pytest import approx

from tracepy import geometry

surface1 = {
    'P': 2.,
    'D': [0., 0., 2.],
    'action': 'reflection',
    'Diam': 2.,
    'diam': 1.,
    'N': 2.,
    'kappa': -2.,
    'c': 3.,
    'name': 'hello',
}

surface2 = {
    'P': [2., 3., 4.],
    'action': 'stop',
    'Diam': 2.,
    'glass': 'F2 schott',
    'c': -3,
    'kappa': 0.,
}

geo = [surface1, surface2]
geo_obj = [geometry(surf) for surf in geo]

def test_params():
    """ Read parameters to make sure nothing was changed. """
    print(geo_obj[0].P)
    assert np.all(geo_obj[0].P == [0., 0., 2.])
    assert np.all(geo_obj[0].D == [0., 0., 2.])
    assert geo_obj[0].action == 'reflection'
    assert geo_obj[0].Diam == 2.
    assert geo_obj[0].diam == 1.
    assert geo_obj[0].N == 2.
    assert geo_obj[0].kappa == -2.
    assert geo_obj[0].c == 3.
    assert geo_obj[0].name == 'hello'

    assert np.all(geo_obj[1].P == [2., 3., 4.])
    assert np.all(geo_obj[1].D == [0., 0., 0.])
    assert geo_obj[1].action == 'stop'
    assert geo_obj[1].Diam == 2.
    assert geo_obj[1].kappa == 0.
    assert geo_obj[1].c == -3.

def test_conics():
    """ Check simple function cases. """
    assert geo_obj[1].conics(np.array([0., 0., 0.]))[0] == 0.0
    assert np.all(geo_obj[1].conics(np.array([0., 0., 0.]))[1] == [0., 0., 1.])

    assert approx(geo_obj[0].conics(np.array([1., 0., 0.]))[0], 0.00002) == -0.55982
    assert approx(geo_obj[0].conics(np.array([1., 0., 0.]))[1][0], 0.00002) == -0.68825
    assert geo_obj[0].conics(np.array([1., 0., 0.]))[1][1] == -0.0
    assert geo_obj[0].conics(np.array([1., 0., 0.]))[1][2] == 1.0

def test_conics_plot():
    """ Check simple Z value returns. """
    assert geo_obj[1].conics_plot(np.array([[0., 0.]]))[0] == 0.
