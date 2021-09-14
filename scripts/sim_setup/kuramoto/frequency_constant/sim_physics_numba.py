from sim_geometry import  *
import carpet.physics.kuramoto_numba as physics

# Physics
period = 31.25            # [ms] period of cilia beat
freq = 2 * np.pi / period # [rad/ms] angular frequency
sin_str = 0.0016 * freq   # coupling strength

# Load frequencies
period = 2 * np.pi / freq
t_max = 2 * period # for the solver

right_side_of_ODE = physics.define_right_side_of_ODE_kuramoto(NN, freq, sin_str, use_numba=True)
solve_cycle = carpet.define_solve_cycle(right_side_of_ODE, t_max, carpet.get_mean_phase)


# Define solve_cycle assuming symmetry classes - used to find fixed points faster.
def define_solve_cycle_class(NN_class, TT_class):
    physics.define_right_side_of_ODE_kuramoto(NN_class, freq, sin_str, use_numba=True)
    return     carpet.define_solve_cycle(right_side_of_ODE, t_max,
                              phi_global_func=carpet.get_mean_phase)