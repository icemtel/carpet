'''
2020-02-03
- NEW procedure - use scipy.optimize.root - much faster
- fixed classes - before there were uneven classes due to jumps from 0 to 2pi!
2020-02-20: upd carpet, upd scipy->numpy

args: k1, k2

===Optimization method discussion===
Use now - `lm` with non-zero starting phase; If it fails - switch to other methods

`hybr` -  NOT OK: strange jump and breaks my solve_cycle
'lm' - OK if start with non-zero mean phase
'broyden1' - was OK, but need more testing

MAYBE: try scipy.optimize.fixed_point
MAYBE: try `krylov` - it is said to be good for large-scale problems.
'''

import sys
import os
import pickle
import numpy as np
import logging

import carpet
import carpet.classes as cc
import carpet.lattice.triangular2 as lattice
import carpet.physics.friction_pairwise as physics

import scipy.optimize as opt

carpet.setup_logging('master.log', mode='a')

## Parameters
# Physics
set_name = 'machemer_1'  # hydrodynamic friction coefficients data set
period = 31.25  # [ms] period
freq = 2 * np.pi / period  # [rad/ms] angular frequency
order_g11 = (8, 0)  # order of Fourier expansion of friction coefficients
order_g12 = (4, 4)
# Geometry
a = 18  # [um]
nx = 6  # number of cilia in x-direction
ny = 6  # must be even
N = nx * ny

## Initialize
# Geometry
L1, L2 = lattice.get_domain_sizes(nx, ny, a)
coords, lattice_ids = lattice.get_nodes_and_ids(nx, ny, a)  # get cilia (nodes) coordinates
N1, T1 = lattice.get_neighbours_list(coords, nx, ny, a)  # get list of neighbours and relative positions
e1, e2 = lattice.get_basis()
get_k = lattice.define_get_k(nx, ny, a)
get_mtwist = lattice.define_get_mtwist(coords, nx, ny, a)

# Physics
gmat_glob, q_glob = physics.define_gmat_glob_and_q_glob(set_name, e1, e2, a, N1, T1, order_g11, order_g12, period)
right_side_of_ODE = physics.define_right_side_of_ODE(gmat_glob, q_glob)
solve_cycle = carpet.define_solve_cycle(right_side_of_ODE, 2 * period, phi_global_func=carpet.get_mean_phase)


def get_classes_rmk(phi, eps=10 ** -8):
    """
    Which cilia have the same phase values?
    define classes of nodes, according to which nodes have the same phase for the given m-twist solution
    :param eps: a small number
    """
    phi = np.array(phi)

    unclassified_list = [(i, phii) for i, phii in enumerate(phi)]

    ix_to_class = np.full_like(phi, fill_value=np.nan, dtype=int)
    class_to_ix = []

    class_id = 0
    while unclassified_list:  # check if not empty
        i, phii = unclassified_list[0]  # take the next oscillator without a class

        class_list = []
        classified_ixs = []  # indices to remove from unclassified_list
        for ix, (j, phij) in enumerate(unclassified_list):
            diff = abs(np.exp(1j * (phii - phij)) - 1)
            if diff < eps:  # Add to the class if difference is small; take into account periodicity
                ix_to_class[j] = class_id
                class_list.append(j)
                classified_ixs.append(ix)

        class_to_ix.append(np.array(class_list, dtype=int))
        unclassified_list = [i for j, i in enumerate(unclassified_list) if j not in classified_ixs]
        class_id += 1
    return ix_to_class, class_to_ix


def get_vec_diff(solve_cycle, tol):
    def vec_diff(phi):
        #        phi_global_end = 0
        #        sol = solve_cycle(phi, tol, phi_global_end=phi_global_end)  # end at global phase = 0
        sol = solve_cycle(phi, tol)  # end at global phase = 0

        phi1 = sol.y.T[-1]

        diff = phi1 - phi - 2 * np.pi - phi.mean()  # force it to prefer zero mean phase
        return diff

    return vec_diff


def find_fixpoint(phi0, tol, mini_tol):
    '''
    v2: optional phi0 as input; otherwise use m-twist
    2019-08-12: - subtract mean phase from the mtwist to work in the same Poincare plane
                - change dynamics for sine coupling
    '''
    # Map to classes
    ix_to_class, class_to_ix = get_classes_rmk(phi0)  # cc.get_classes(phi0)
    nclass = len(class_to_ix)
    # Get classes representatives
    # Get one cilium from each of cilia classes
    unique_cilia_ids = np.array([class_to_ix[iclass][0] for iclass in range(nclass)], dtype=np.int64)
    # Get neighbours
    N1_class, T1_class = cc.get_neighbours_list_class(unique_cilia_ids, ix_to_class, N1, T1)
    # Dynamics
    gmat_glob_class, q_glob_class = physics.define_gmat_glob_and_q_glob(set_name, e1, e2, a, N1_class, T1_class,
                                                                        order_g11, order_g12, period)
    right_side_of_ODE_class = physics.define_right_side_of_ODE(gmat_glob_class, q_glob_class)
    solve_cycle_class = carpet.define_solve_cycle(right_side_of_ODE_class, 2 * period,
                                                  phi_global_func=carpet.get_mean_phase)

    # Optimize!
    phi0_class = phi0[unique_cilia_ids]  # I.Co.
    vec_diff = get_vec_diff(solve_cycle_class, tol)
    res = opt.root(vec_diff, phi0_class, tol=mini_tol, method='lm')

    if not res.success:
        logging.warning(f'Did not converge, k1,k2=({k1},{k2})')

    fixpoint_class = np.array(res.x)
    fixpoint = fixpoint_class[ix_to_class]  # from classes to cilia ix

    return fixpoint


### Main
k1, k2 = int(sys.argv[1]), int(sys.argv[2])
dirname = os.path.dirname(__file__)  # sys.argv[3]
tol = 10 ** -8


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
phi0 = phi0  # - carpet.get_mean_phase(phi0)  # Test without this line - it makes it worse at least in 1 case

if k1 == 0 and k2 == 0:
    # assuming that symmetries of fixpoint are preserved,
    # there is only one possibility up to a constant shift: (0,0,0..,0)
    fixpoint = phi0
else:
    fixpoint = find_fixpoint(phi0, tol, tol)

filename = "fixpoint_k1={}_k2={}.pkl".format(k1, k2)
dump_object(fixpoint, filename)
