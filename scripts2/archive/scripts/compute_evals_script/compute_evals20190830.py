"""
20190830: lattice_triangular2


script takes arguments
- k1, k2
- directory

- Load data needed from the directory
- Compute matrix L and eigenvalues for the fixpoint, located, corresponding to (k1,k2)-twist

TODO: more natural way to pass data than pickles?
TODO: use cwd instead of dirname?
"""
import sys
import logging
import os
import pickle
import scipy as sp
import scipy.linalg as lin
# import plotting routines

import carpet
import carpet.lattice_triangular2 as lattice

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


delta0 = 0.005
tol = 10 ** -8

set_name = 'machemer_1' # which hydrodynamic coefficients to use
order_g11 = (8,0)
order_g12 = (4,4)
T = 31.25 # [ms] period of cilia beat
freq = 2 * sp.pi / T # [rad/ms] angular frequency

# Geometry
nx = 12
ny = 12  # even number
N = nx * ny
a = 18  # [um] lattice spacing

#
# delta0 = load_object('delta0.pkl')
# tol = load_object('tol.pkl')
# set_name = load_object('set_name.pkl') # which hydrodynamic coefficients to use
# order_g11 = load_object('order_g11.pkl')
# order_g12 = load_object('order_g12.pkl')
# T = load_object('T.pkl')
# freq = 2 * sp.pi / T # [rad/ms] angular frequency
#
# # Geometry
# nx = load_object('nx.pkl')
# ny = load_object('ny.pkl')
# N = nx * ny
# a = load_object('a.pkl')


## Initialize
# Geometry
L1,L2 = lattice.get_domain_sizes(nx,ny ,a)
coords, lattice_ids = lattice.get_nodes_and_ids(nx, ny, a)
N1, T1 = lattice.get_neighbours_list(coords, nx, ny, a)
get_k = lattice.define_get_k(nx, ny, a)
get_mtwist = lattice.define_get_mtwist(coords, nx, ny, a)

# Physics
gmat_glob, q_glob = lattice.define_gmat_glob_and_q_glob(set_name, a, N1, T1,order_g11, order_g12, T)
right_side_of_ODE = lattice.define_right_side_of_ODE(gmat_glob, q_glob)
solve_cycle = carpet.define_solve_cycle(right_side_of_ODE,2 * T, carpet.get_mean_phase)



# Load fixpoints
filename = 'fixpoint_dict.pkl'
fixpoint_dict = load_object(filename)


def get_fixpoint(k1,k2):
    return sp.array(fixpoint_dict[k1,k2])

fixpoint = get_fixpoint(k1,k2)

## Eigenvalues functions
"""
2019-07-31: choose any set of N-1 perturbations with zero mean phase
2019-08-27: eigenvalues of L-I -> eigenvalues of ln(L) [logarithms of eigenvalues of L]
"""


def calc_sort_log_evals_evecs(L):
    '''
    :return: logairthm of eigenvalues, eigenvectors
    '''
    # Ljapunov spectrum
    evals, evecs = sp.linalg.eig(L)
    evals = sp.log(evals)
    # Sort eigenvalues
    idx = evals.argsort()
    evals = evals[idx]
    evecs = evecs[:, idx]
    return evals, evecs


def get_L3_general(k1, k2, Delta0s, tol):
    """
    Input: N-1 perturbation vectors with zero mean phase

    N-th perturbation will be added automatically:
    Since we have already tested that it's neutral - just set last Delta1 = Delta0 = (1,1,1...1)
    """

    Delta1s = []  # list of deviation from the fixpoint after one cycle
    fixpoint = get_fixpoint(k1, k2)
    for Delta0 in Delta0s:  # index of single node that was initially perturbed
        # Initial condition
        phi0 = fixpoint + Delta0

        assert abs(carpet.get_mean_phase(Delta0)) < tol

        # Solve
        solution = solve_cycle(phi0, tol)
        phi1 = solution.y.T[-1]
        # Fill matrix row
        Delta1 = (phi1 - fixpoint - 2 * sp.pi)
        Delta1s.append(Delta1)

    ## Out of Poincare-plane perturbation - make it stay neutral
    D0 = sp.array([Delta0 for Delta0 in Delta0s] + [sp.ones([N])])  # Truncated vectors as rows of a matrix
    D1 = sp.array([Delta1 for Delta1 in Delta1s] + [sp.ones([N])])  # Truncated vectors as rows of a matrix

    ## Get L
    L = lin.solve(D0, D1).transpose()  # Transpose - to work with vectors as columns again

    return L


def get_L_single(k1, k2, delta0, tol):
    """
    N-1 single-node perturbations (+ change other nodes to preserve mean phase)
    - Normalized by L_infinity norm
    """

    def get_single_node_perturbation(ix, delta0):
        Delta = sp.zeros(N) - 1 / (N - 1)
        Delta[ix] = 1
        return Delta * delta0

    Delta0s = [get_single_node_perturbation(ix, delta0) for ix in range(0, N - 1)]  # get basis, multiply by delta0
    return get_L3_general(k1, k2, Delta0s, tol)


def get_mtwist_trig_basis(delta0=1, phi0=0):
    '''
    ad hoc solution for 2D
    Go through all possible cosine,sine of mtwists, keep only the ones which are orthogonal to each other
    2019-08-28: phi0 - phase shift all mtwists
        Added this one from 2D: can shift by angle
        Checked: equivalent to the old one when phi0=0
    '''

    def zero_vec(vec, eps=10 ** -8):
        return lin.norm(vec) * N ** (-1 / 2) < eps

    def orthogonal(vec, basis, eps=10 ** -8):
        '''
        If vector is coaligned with a vector from the set - return True, else False
        '''
        for b in basis:
            if abs(vec @ b.conj()) > eps * lin.norm(vec) * lin.norm(b):
                return False
        return True

    basis = []
    for k1 in range(nx):
        for k2 in range(ny):
            if k1 == 0 and k2 == 0:
                continue

            cos_mtwist = sp.cos(get_mtwist(k1, k2) + phi0)
            sin_mtwist = sp.sin(get_mtwist(k1, k2) + phi0)

            if not zero_vec(cos_mtwist) and orthogonal(cos_mtwist, basis):
                basis.append(cos_mtwist)

            if not zero_vec(sin_mtwist) and orthogonal(sin_mtwist, basis):
                basis.append(sin_mtwist)

    Delta0s = [delta0 * b for b in basis]
    assert len(Delta0s) == N - 1
    return Delta0s


def get_L_mtwist(k1, k2, delta0, tol):
    """
    N-1 perturbation - cosines and sines of m-twists
    Nth - orthogonal to the Poincare plane; construct L s.t. this is a neutral perturbation
    """
    Delta0s = get_mtwist_trig_basis(delta0)
    return get_L3_general(k1, k2, Delta0s, tol)


def calc_evals_evecs_mtwist(k1, k2, delta0, tol):
    '''
    evecs: eigenvectors as columns!
    '''
    L = get_L_mtwist(k1, k2, delta0, tol)
    return calc_sort_log_evals_evecs(L)


def fill_evals_evecs_dict_mtwist(k1, k2, delta0, tol, evals_dict, evecs_dict):
    evals, evecs = calc_evals_evecs_mtwist(k1, k2, delta0, tol)

    evals_dict[(k1, k2)] = evals
    evecs_dict[(k1, k2)] = evecs
    logging.info("Finished: k1={} k2={}".format(k1, k2))

###


L_mtwist = get_L_mtwist(k1, k2, delta0, tol)

filename = "L_mtwist_k1={}_k2={}.pkl".format(k1,k2)
dump_object(L_mtwist, filename)

# L_log_lin = L_mtwist - sp.eye(N)
# evals, evecs =  calc_sort_evals_evecs(L_log_lin)
