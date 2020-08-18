'''
6x6 kuramoto with zero coupling (test)

'''

import numpy as np

import carpet
import carpet.lattice.triangular as lattice

## Parameters
# Physics
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
import carpet.physics.kuramoto as physics

sin_str = 0.0016 * freq
coupling = physics.define_sine_coupling(sin_str)
right_side_of_ODE = physics.define_right_side_of_ODE(coupling, freq, NN, TT)
solve_cycle = carpet.define_solve_cycle(right_side_of_ODE, 2 * period, carpet.get_mean_phase)
