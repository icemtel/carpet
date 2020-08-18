'''
12x12 1-N triangular

Renate+numba:
tol=1e-6 => 300ms per cycle
tol=e1-8 => 687ms per cycle
stochastic dt=0.01*T => 241 ms per period
'''

import numpy as np

import carpet
import carpet.lattice.triangular as lattice
import carpet.physics.friction_pairwise as physics

## Parameters
# Physics
set_name = 'machemer_1' # which hydrodynamic coefficients to use
order_g11 = (4,0)
order_g12 = (4,4)
period = 31.25            # [ms] period of cilia beat
freq = 2 * np.pi / period # [rad/ms] angular frequency

# Geometry
nx = 12
ny = 12  # even number
N = nx * ny
a = 18  # [um] lattice spacing

## Initialize
# Geometry

L1,L2 = lattice.get_domain_sizes(nx,ny ,a)
coords, lattice_ids = lattice.get_nodes_and_ids(nx, ny, a)  # get cilia (nodes) coordinates

distances = [1]
NN, TT = lattice.get_neighbours_list(coords, nx, ny, a, distances)
e1, e2 = lattice.get_basis()
get_k = lattice.define_get_k_fbz(nx, ny, a)
get_mtwist = lattice.define_get_mtwist(coords, nx, ny, a)

# Physics
gmat_glob, q_glob = physics.define_gmat_glob_and_q_glob(set_name, e1, e2, a, NN, TT, order_g11, order_g12, period)
right_side_of_ODE = physics.define_right_side_of_ODE(gmat_glob, q_glob)
solve_cycle = carpet.define_solve_cycle(right_side_of_ODE, 2 * period, phi_global_func=carpet.get_mean_phase)


# Define solve_cycle assuming symmetry classes - used to find fixed points faster.
def define_solve_cycle_class(NN_class, TT_class):
    gmat_glob_class, q_glob_class = physics.define_gmat_glob_and_q_glob(set_name, e1, e2, a, NN_class, TT_class,
                                                                        order_g11, order_g12, period)
    right_side_of_ODE_class = physics.define_right_side_of_ODE(gmat_glob_class, q_glob_class)
    return     carpet.define_solve_cycle(right_side_of_ODE_class, 2 * period,
                              phi_global_func=carpet.get_mean_phase)
