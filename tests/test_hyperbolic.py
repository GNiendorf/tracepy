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

def test_rms_hyperbolic():
    geo = [back_lens, lens, stop]
    ray_group = tp.ray_plane(geo, [0., 0., 0.], 1.1, d=[0., 0., 1.], nrays=100)
    rms = tp.spot_rms(geo, ray_group)
    assert rms == 0.