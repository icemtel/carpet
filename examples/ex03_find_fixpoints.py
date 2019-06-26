'''
Wasn't checked after copying from a notebook.
'''

import scipy as sp
from scipy.linalg import norm
import scipy.optimize as opt


def get_optimization_target(solve_cycle, tol):
    def get_distance_squared(phi0, info=None):
        '''
        :param info: Dictionary, must have following fields:
                     info["print"] - True/False
                     info['num_evals'] initialize with 0, then this parameter will change
                     info['initial'] - mtwist or other initial point, prints a distance to that point
        '''
        sol = solve_cycle(phi0, tol)
        phi1 = sol.y.T[-1]
        distance_squared = sp.sum((phi1 - phi0 - 2 * sp.pi) ** 2)
        if info is not None:
            if info['print'] == True:
                info['num_evals'] += 1
                distance_to_initial = norm(phi1 - info['initial'] - 2 * sp.pi)
                if info['num_evals'] == 1 or info['num_evals'] % 500 == 0:
                    print("{:<15}{:<20}{:<20}".format('num_evals', '||L(phi0) - phi0||', 'Distance from initial'))
                if info['num_evals'] == 1 or info['num_evals'] % 50 == 0:
                    print("{:<15}{:<20.3e}{:<20.3e}".format(
                        info['num_evals'], distance_squared ** (1 / 2), distance_to_initial))

        return distance_squared

    return get_distance_squared

if __name__ is "__main__":
    import time
    import carpet.parallel_with_threads as pwt
    import carpet.dynamics as dynamics
    import carpet.lattice_triangular as lattice

    a = 18
    nx = 6
    ny = 6
    tol = 10 ** -6
    mini_tol = 10 ** -6
    set_name = 'machemer_1'
    order_g11 = (8, 0)
    order_g12 = (4, 4)
    T = 31.25

    coords, lattice_ids = lattice.get_nodes_and_ids(nx, ny, a)
    N1, T1 = lattice.get_neighbours_list(coords, nx, ny, a)

    gmat_glob, q_glob = lattice.define_gmat_glob_and_q_glob(set_name, a, N1, T1, order_g11, order_g12, T)
    right_side_of_ODE = lattice.define_right_side_of_ODE(gmat_glob, q_glob)

    mean_phase = lambda phi: sp.mean(phi)

    solve_cycle = dynamics.define_solve_cycle(right_side_of_ODE, 2 * T, mean_phase)

    phi_global_func = dynamics.get_mean_phase
    distance_squared = get_optimization_target(phi_global_func, tol)

    fixpoint_dict = {}


    def function_to_run(k1, k2):  # job for parallel computing
        start = time.time()
        phi0 = lattice.get_mtwist_phi(k1, k2)

        res = opt.minimize(distance_squared, phi0, tol=mini_tol)
        time_spent = time.time() - start
        print("{:<7} {:<7} {:<15.3g} {:<15.3e}".format(k1, k2, time_spent, norm(phi0 - res.x)))

        fixpoint_dict[(k1, k2)] = sp.array(res.x)


    input_args = []
    for k1 in range(nx):
        for k2 in range(ny):
            input_args.append((k1, k2))

    print('k1\tk2\tTime Spent\tDistance from mtwist')
    pwt.run_parallel(3, function_to_run, list_of_args=input_args)
