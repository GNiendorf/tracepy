import pytest

import numpy as np
import matplotlib.pyplot as plt

import tracepy as tp

lens = {
    'action': 'refraction',
    'P': 2.,
    'kappa': -1.25,
    'c': -2,
    'Diam': 2.2
}

back_lens = {
    'action': 'refraction',
    'P': 1.,
    'N': 1.5,
    'Diam': 2.2
}

stop = {
    'action': 'stop',
    'P': 3.,
    'Diam': 2.2
}

def test_rms():
    geo = [back_lens, lens, stop]
    ray_group = tp.ray_plane(geo, [0., 0., 0.], 1.1, d=[0.,0.,1.], nrays=100)
    rms = tp.spotdiagram(geo, ray_group, optimizer=True)
    assert pytest.approx(rms, 0., abs=1e-14)

