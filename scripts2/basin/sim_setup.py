'''
6x6 triangular rotated

on Renate: (.. checked during other simulations running - might be faster)
Solve 1 cycle:
- 250 ms with tol=1e-6
- 575 ms with tol=1e-8

=> 1500 cycles, 1000 traj, 16 CPU => 6.5h
Taurus w/ 16 CPUS: 11h22m for 1200 traj; 1500 cycles, tol=1e-6
'''

import numpy as np
import pickle
import carpet
import carpet.lattice.triangular2 as lattice
import carpet.physics.friction_pairwise as physics

## Parameters
# Physics
set_name = 'machemer_1' # which hydrodynamic coefficients to use
order_g11 = (8,0)
order_g12 = (4,4)
period = 31.25            # [ms] period of cilia beat
freq = 2 * np.pi / period # [rad/ms] angular frequency

# Geometry
nx = 6
ny = 6  # even number
N = nx * ny
a = 18  # [um] lattice spacing

## Initialize
# Geometry

L1,L2 = lattice.get_domain_sizes(nx,ny ,a)
coords, lattice_ids = lattice.get_nodes_and_ids(nx, ny, a)  # get cilia (nodes) coordinates

distances = [1] #, 3 ** 0.5]
NN, TT = lattice.get_neighbours_list(coords, nx, ny, a, distances)
e1, e2 = lattice.get_basis()
get_k = lattice.define_get_k(nx, ny, a)
get_mtwist = lattice.define_get_mtwist(coords, nx, ny, a)

# Physics
gmat_glob, q_glob = physics.define_gmat_glob_and_q_glob(set_name, e1, e2, a, NN, TT, order_g11, order_g12, period)
right_side_of_ODE = physics.define_right_side_of_ODE(gmat_glob, q_glob)
solve_cycle = carpet.define_solve_cycle(right_side_of_ODE, 2 * period, phi_global_func=carpet.get_mean_phase)
solve_cycle_back = carpet.define_solve_cycle(right_side_of_ODE, 2 * period, phi_global_func=carpet.get_mean_phase,
                                             backwards=True)

# Load fixpoints
filename = 'fixpoint_dict.pkl' #'fixpoint_dict_nx={}_ny={}_tol={:.3E}.pkl'.format(nx,ny,tol)

with open(filename, 'rb') as f:
    fixpoint_dict = pickle.load(f)


def get_fixpoint(k1,k2):
    return np.array(fixpoint_dict[k1,k2])
