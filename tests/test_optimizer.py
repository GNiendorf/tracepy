import numpy as np
import tracepy as tp

lens = {
    'name': 'lens',
    'action': 'refraction',
    'P': 2.,
    'kappa': -0.004,
    'c': -0.5,
    'Diam': 2.2
}

back_lens = {
    'name': 'back_lens',
    'action': 'refraction',
    'P': 1.,
    'c': 0.5,
    'kappa': -0.004,
    'N': 1.5,
    'Diam': 2.2
}

stop = {
    'name': 'stop',
    'action': 'stop',
    'P': 3.,
    'Diam': 2.2
}

# Initial geometry
geo_1 = [back_lens, lens, stop]

ray_group = tp.ray_plane(geo_1, [0., 0., 0.], 1.1, d=[0., 0., 1.], nrays=100)

# Define variable elements for optimization
vlens = {
    'name': 'stop',
    'vary': ['P']
}

vlens2 = {
    'name': 'back_lens',
    'vary': ['c', 'kappa', 'N']
}

vlens3 = {
    'name': 'lens',
    'vary': ['c', 'kappa']
}

vary_list = [vlens, vlens2, vlens3]

def test_optimizer():
    try:
        geo_opt = tp.optimize(geo_1, vary_list, typeof='least_squares', max_iter=10)
    except Exception as e:
        assert False, f"Optimization raised an exception: {e}"

    try:
        ray_group_3 = tp.ray_plane(geo_opt, [0., 0., 0.], 1.1, d=[0., 0., 1.], nrays=1000)
        tp.rayaberration(geo_1, ray_group_3)
    except Exception as e:
        assert False, f"rayaberration raised an exception: {e}"