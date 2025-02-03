import numpy as np
import tracepy as tp

mirror = {
    'action': 'reflection',
    'D': np.array([0., np.pi, 0.]),
    'P': np.array([0., 0., 10.]),
    'kappa': 0.,
    'c': 0.05,
    'Diam': 3.6,
}

small_mirror = {
    'action': 'reflection',
    'D': np.array([0., -np.pi/4., 0.]),
    'P': np.array([0., 0., 2.]),
    'Diam': 1.4,
}

selector = {
    'action': 'stop',
    'P': np.array([0., 0., 0.1]),
    'Diam': 3.6,
    'diam': 1.4
}

stop = {
    'action': 'stop',
    'P': np.array([0., -4.6, 2.]),
    'D': np.array([0., np.pi/2., 0.]),
    'Diam': 0.8
}

lensy = -4.
thickness = 0.1

lens1 = {
    'action': 'refraction',
    'P': np.array([0., lensy, 2.]),
    'D': np.array([0., 3*np.pi/2., 0.]),
    'Diam': 0.8,
    'kappa': 0.,
    'glass': 'ZF4 cdgm',
    'c': 0.3
}

lens2 = {
    'action': 'refraction',
    'P': np.array([0., lensy - thickness, 2.]),
    'D': np.array([0., 3*np.pi/2., 0.]),
    'Diam': 0.8,
    'kappa': 0.,
    'c': -0.3
}

def test_plotting_functions():
    geo = [selector, mirror, small_mirror, lens1, lens2, stop]
    ray_group = tp.ray_plane(geo, [0., 0., 0.], 1.8, d=[0., 0., 1.], nrays=50, wvl=0.55)

    # Test plotyz
    try:
        tp.plotyz(geo, ray_group)
    except Exception as e:
        assert False, f"plotyz raised an exception: {e}"

    # Test plotxz
    try:
        tp.plotxz(geo, ray_group)
    except Exception as e:
        assert False, f"plotxz raised an exception: {e}"

    # Test spotdiagram
    try:
        tp.spotdiagram(geo, ray_group)
    except Exception as e:
        assert False, f"spotdiagram raised an exception: {e}"

    # Test plotobject
    try:
        tp.plotobject(geo, ray_group)
    except Exception as e:
        assert False, f"plotobject raised an exception: {e}"