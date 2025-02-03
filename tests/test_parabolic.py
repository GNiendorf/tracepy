import tracepy as tp

mirror = {
    'action': 'reflection',
    'P': 1.5,
    'kappa': 0.,
    'c': -0.5,
    'Diam': 2.2,
}

stop = {
    'action': 'stop',
    'P': 0.5,
    'Diam': 0.2
}

def test_rms_parabolic():
    geo = [mirror, stop]
    ray_group = tp.ray_plane(geo, [0., 0., -1.5], 1.1, d=[0., 0., 1.], nrays=100)
    rms = tp.spot_rms(geo, ray_group)
    assert rms == 0.