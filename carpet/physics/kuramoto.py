"""
Kuramoto-type model.

- Sinusoidal coupling
- Cosine + sine coupling
- Possible to define any other
"""
import numpy as np
import math


def define_sine_coupling(sin_str):
    '''
    :sin_str: coefficient in front of sine, units of [rad/s]
             In case of two cilia: positive => in-phase synchronization is stable
    '''
    def coupling(phi1, phi2, translation=None):
        c = sin_str * math.sin(phi2 - phi1)
        return c

    return coupling


def define_sine_cosine_coupling(sin_str, cos_str):
    def coupling(phi1, phi2, translation=None):
        c = sin_str * math.sin(phi2 - phi1) + cos_str * math.sin(phi2 - phi1)
        return c

    return coupling


def define_right_side_of_ODE(coupling_func, freq, neighbours_list, translations_list):
    '''
    :param coupling_func: function of pairwise interactions with signature f(phi1,phi2, translation)
    :param freq: frequency of a single oscillator
    :param neighbours_list: list of lists of neighbours indices
    :param translations_list: UNUSED - keep for generalizations
    :return:
    '''

    def right_side_of_ODE(t, phi):
        # Self-driving (multiply by frequency later
        right_side = np.full(len(neighbours_list), fill_value=freq)
        # Interactions
        # MAYBE: Optimize -  sum array for; or numba
        for ix, (neighbours, translations) in enumerate(zip(neighbours_list, translations_list)):
            for ix2, translation in zip(neighbours, translations):
                right_side[ix] += coupling_func(phi[ix], phi[ix2], translation)
        return right_side

    return right_side_of_ODE


if __name__ == '__main__':
    import carpet.lattice.triangular as lattice
    import carpet.visualize as vis
    import matplotlib.pyplot as plt
    import carpet, time

    a = 18  # [um]
    period = 31.25  # [ms] period
    freq = 2 * np.pi / period
    nx = 3
    ny = 4  # must be even
    tol = 10 ** -6  # solver tolerance

    coords, lattice_ids = lattice.get_nodes_and_ids(nx, ny, a)  # get cilia (nodes) coordinates
    N1, T1 = lattice.get_neighbours_list(coords, nx, ny, a)  # get list of neighbours and relative positions

    sin_str = 0.001 * freq
    coupling = define_sine_coupling(sin_str)
    right_side_of_ODE = define_right_side_of_ODE(coupling, freq, N1, T1)
    solve_cycle = carpet.define_solve_cycle(right_side_of_ODE, 2 * period, phi_global_func=carpet.get_mean_phase)

    # Solve
    phi0 = np.zeros([len(coords)])  # initial condition
    phi0[5] += 3  # perturb
    sol = solve_cycle(phi0, tol)

    start = time.time()  # track CPU time
    solution = solve_cycle(phi0, tol)
    time_spent = time.time() - start
    print("Time spent", time_spent)

    # Get result
    phi1 = solution.y.T[-1]
    print("Change in phase after one cycle:", phi1 - phi0 - 2 * np.pi)
    # Visualize
    vals = (phi1 - phi0 - 2 * np.pi)
    plt.title("How much phases changed after one cycle")
    vis.plot_nodes(coords, phi=vals, vmin=vals.min(), vmax=vals.max(),cmap='jet')
    vis.plt.show()

    # Plot as a function of time
    ys = sol.y[5]
    ts = sol.t
    plt.plot(ts, ys  - freq * ts, 'o')
    plt.show()