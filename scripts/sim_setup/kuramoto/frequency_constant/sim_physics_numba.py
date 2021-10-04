'''
NB: using the friction dataset for triangular lattice - for now.
'''

import numpy as np
from sim_geometry import  *
import carpet.physics.kuramoto_numba as physics

## Parameters
# Physics
period = 31.25            # [ms] period of cilia beat
freq = 2 * np.pi / period # [rad/ms] angular frequency
t_max = 2 * period # for the solver
sin_str = 3.45e-4

right_side_of_ODE = physics.define_right_side_of_ODE_kuramoto(NN, freq * np.ones(N), sin_str, use_numba=True)
solve_cycle = carpet.define_solve_cycle(right_side_of_ODE, t_max, carpet.get_mean_phase)


# Define solve_cycle assuming symmetry classes - used to find fixed points faster.
def define_solve_cycle_class(NN_class, TT_class):
    right_side_of_ODE_class = physics.define_right_side_of_ODE_kuramoto(NN_class, freq * np.ones(N), sin_str, use_numba=True)
    return carpet.define_solve_cycle(right_side_of_ODE_class, t_max, carpet.get_mean_phase)



