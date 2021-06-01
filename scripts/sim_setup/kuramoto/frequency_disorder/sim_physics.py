'''
To load frequencies

Renate+numba:
tol=1e-6 => ms per cycle
tol=e1-8 => ms per cycle
stochastic dt=0.01*T =>  ms per period
'''
from sim_geometry import  *
import carpet.physics.kuramoto as physics

# Physics
period0 = 31.25            # [ms] period of cilia beat
freq0 = 2 * np.pi / period0 # [rad/ms] angular frequency
sin_str = 0.0016 * freq0  # coupling strength

# Load frequencies
freq = np.load('freq.npy') # vector!
period = 2 * np.pi / freq
t_max = 2 * period.mean() # for the solver

coupling = physics.define_sine_coupling(sin_str)
right_side_of_ODE = physics.define_right_side_of_ODE(coupling, freq, NN, TT)
solve_cycle = carpet.define_solve_cycle(right_side_of_ODE, t_max, carpet.get_mean_phase)


# Define solve_cycle assuming symmetry classes - used to find fixed points faster.
def define_solve_cycle_class(NN_class, TT_class):
    right_side_of_ODE = physics.define_right_side_of_ODE(coupling, freq, NN_class, TT_class)
    return     carpet.define_solve_cycle(right_side_of_ODE, t_max,
                              phi_global_func=carpet.get_mean_phase)