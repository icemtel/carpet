"""
Kuramoto-type model.

- Sinusoidal coupling
- accelerated with numba

"""
import numpy as np
import math
import numba

def njit_wrapper(use_numba, *args, **kwargs):
    '''
    Wrapper for numba njit decorator
    If :param use_numba: is True  -> will use @njit with given args and kwargs.
                            False -> decorator does not do anything
    :return: a decorator
    '''
    if use_numba:
        from numba import njit

    def decorator(func):
        if not use_numba:
            # Return the function unchanged, not decorated.
            return func
        else:
            return njit(func, *args, **kwargs)

    return decorator


def define_right_side_of_ODE_kuramoto(neighbours_indices, omega, sin_str, use_numba=True):
    '''
    :param use_numba: if True, accelerate code with numba (5-50 times faster);
                      requires the same number of neighbours for every oscillator
    :return:
    '''
    N = len(omega)

    if use_numba:
        neighbours_indices = np.array(neighbours_indices)  # numba doesn't like list of lists
        # neighbours_indices = [tuple(x) for x in neighbours_indices]  # -> tested; does not work

    @njit_wrapper(use_numba)
    def q_glob(t, phi):
        """
        :param neighbours_indices: list of neighbours for each oscillator
        :param omega: vector of frequencies rad/s
        :param sin_str: coupling strength K/m in rad/s, float
        :param phi: vector of phases, corresponding to each of the cilia
        :param t: vector of time
        :return: right side of differential equation as a function of t and phi
        """
        q = np.zeros(N)
        for i in range(N):
            neighbours = neighbours_indices[i]
            coupling = 0
            for neighbour_index in neighbours:
                coupling += sin_str * math.sin(phi[i] - phi[neighbour_index])
            q[i] = omega[i] - coupling
        return q

    return q_glob

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
    N = nx * ny
    tol = 10 ** -6  # solver tolerance

    coords, lattice_ids = lattice.get_nodes_and_ids(nx, ny, a)  # get cilia (nodes) coordinates
    N1, T1 = lattice.get_neighbours_list(coords, nx, ny, a)  # get list of neighbours and relative positions

    sin_str = 0.001 * freq
    right_side_of_ODE = define_right_side_of_ODE_kuramoto(N1, freq * np.ones(N), sin_str, use_numba=False)
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