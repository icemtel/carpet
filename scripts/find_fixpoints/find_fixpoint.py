import sys
import logging
import os
import pickle
import scipy as sp
import scipy.linalg as lin

import carpet
import carpet.classes as cc
import carpet.lattice_triangular2 as lattice

import scipy.optimize as opt


## Parameters
# Physics
set_name = 'machemer_1' # which hydrodynamic coefficients to use
order_g11 = (8,0)
order_g12 = (4,4)
period = 31.25 # [ms] period of (a single cilium) beat
# freq = 2 * sp.pi / period # [rad/ms] angular frequency

# Geometry
nx = 12 # even number
ny = 12
N = nx * ny
a = 18  # [um] lattice spacing
## Initialize
# Geometry
L1,L2 = lattice.get_domain_sizes(nx,ny ,a)
coords, lattice_ids = lattice.get_nodes_and_ids(nx, ny, a)
N1, T1 = lattice.get_neighbours_list(coords, nx, ny, a) # neighbours, translations
get_k = lattice.define_get_k(nx, ny, a)
get_mtwist = lattice.define_get_mtwist(coords, nx, ny, a)

# Physics
gmat_glob, q_glob = lattice.define_gmat_glob_and_q_glob(set_name, a, N1, T1, order_g11, order_g12, period)
right_side_of_ODE = lattice.define_right_side_of_ODE(gmat_glob, q_glob)
solve_cycle = carpet.define_solve_cycle(right_side_of_ODE,2 * period, carpet.get_mean_phase)



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
                distance_to_initial = lin.norm(phi1 - info['initial'] - 2 * sp.pi)
                if info['num_evals'] == 1 or info['num_evals'] % 500 == 0:
                    print("{:<15}{:<20}{:<20}".format('num_evals', '||L(phi0) - phi0||', 'Distance from initial'))
                if info['num_evals'] == 1 or info['num_evals'] % 50 == 0:
                    print("{:<15}{:<20.3e}{:<20.3e}".format(
                        info['num_evals'], distance_squared ** (1 / 2), distance_to_initial))

        return distance_squared

    return get_distance_squared


def find_fixpoint(phi0, tol, mini_tol):
    '''
    v2: optional phi0 as input; otherwise use m-twist
    2019-08-12: - subtract mean phase from the mtwist to work in the same Poincare plane
                - change dynamics for sine coupling
    '''
    phi0 = phi0 - carpet.get_mean_phase(phi0)  # only search for fixpoints with zero mean phase

    # Map to classes
    ix_to_class, class_to_ix = cc.get_classes(phi0)
    nclass = len(class_to_ix)
    # Get classes representatives
    # Get one cilium from each of cilia classes
    unique_cilia_ids = sp.array([class_to_ix[iclass][0] for iclass in range(nclass)], dtype=sp.int64)
    # Get neighbours
    N1_class, T1_class = cc.get_neighbours_list_class(unique_cilia_ids, ix_to_class, N1, T1)
    # Dynamics
    gmat_glob_class, q_glob_class = lattice.define_gmat_glob_and_q_glob(set_name, a, N1_class, T1_class,
                                                                        order_g11, order_g12, period)
    right_side_of_ODE_class = lattice.define_right_side_of_ODE(gmat_glob_class, q_glob_class)
    solve_cycle_class = carpet.define_solve_cycle(right_side_of_ODE_class, 2 * period, carpet.get_mean_phase)

    # Optimize!
    phi0_class = phi0[unique_cilia_ids]  # I.Co.
    distance_squared = get_optimization_target(solve_cycle_class, tol)
    res = opt.minimize(distance_squared, phi0_class, tol=mini_tol)
    fixpoint_class = sp.array(res.x)

    # Classes to cilia ix
    fixpoint = fixpoint_class[ix_to_class]

    return fixpoint


### Main
k1,k2 = int(sys.argv[1]), int(sys.argv[2])
dirname = os.path.dirname(__file__) # sys.argv[3]


def dump_object(obj, filename):
    filename = os.path.join(dirname, filename)
    print(filename)
    with open(filename, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_object(filename):
    filename = os.path.join(dirname, filename)
    with open(filename, 'rb') as f:
        obj = pickle.load(f)
    return obj




phi0 = get_mtwist(k1, k2)
phi0 = phi0 - carpet.get_mean_phase(phi0)  # only search for fixpoints with zero mean phase

if k1 == 0 and k2 == 0:
    # assuming that symmetries of fixpoint are preserved,
    # there is only one possibility up to a constant shift: (0,0,0..,0)
    fixpoint = phi0
else:
    tol = 10 ** - 4
    fixpoint = find_fixpoint(phi0, tol, tol)

    tol = 10 ** - 6
    fixpoint = find_fixpoint(fixpoint, tol, tol)

    tol = 10 ** - 8
    fixpoint = find_fixpoint(fixpoint, tol, tol)


filename = "fixpoint_k1={}_k2={}.pkl".format(k1,k2)
dump_object(fixpoint, filename)