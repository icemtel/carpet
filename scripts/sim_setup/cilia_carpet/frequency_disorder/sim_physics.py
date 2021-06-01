'''
- Loads a vector of frequencies `freq.npy`
'''

import numpy as np
from sim_geometry import  *
import carpet.physics.friction_pairwise as physics

## Parameters
# Physics
set_name = 'machemer_4' # which hydrodynamic coefficients to use
order_g11 = (4,0)
order_g12 = (4,4)

# Physics
# period = 31.25            # [ms] period of cilia beat
# freq = 2 * np.pi / period # [rad/ms] angular frequency
freq = np.load('freq.npy') # vector!
period = 2 * np.pi / freq
t_max = 2 * period.mean() # for the solver

gmat_glob, q_glob = physics.define_gmat_glob_and_q_glob(set_name, e1, e2, a, NN, TT, order_g11, order_g12, period)
right_side_of_ODE = physics.define_right_side_of_ODE(gmat_glob, q_glob)
solve_cycle = carpet.define_solve_cycle(right_side_of_ODE, t_max, phi_global_func=carpet.get_mean_phase)


# Define solve_cycle assuming symmetry classes - used to find fixed points faster.
def define_solve_cycle_class(NN_class, TT_class):
    gmat_glob_class, q_glob_class = physics.define_gmat_glob_and_q_glob(set_name, e1, e2, a, NN_class, TT_class,
                                                                        order_g11, order_g12, period)
    right_side_of_ODE_class = physics.define_right_side_of_ODE(gmat_glob_class, q_glob_class)
    return     carpet.define_solve_cycle(right_side_of_ODE_class, 2 * period,
                              phi_global_func=carpet.get_mean_phase)
