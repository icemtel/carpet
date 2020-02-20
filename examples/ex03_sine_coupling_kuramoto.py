import numpy as np

import carpet
import carpet.lattice.triangular as lattice
import carpet.physics.kuramoto as physics

carpet.setup_logging('integrate_trajectory.log')

## Parameters
# Physics
period = 31.25 # [ms] period of cilia beat; freq = 2 * sp.pi / period [rad/ms]
freq = 2 * np.pi / period
sin_str = 0.003 * freq
# Geometry
nx = 6
ny = 6 # even number
N = nx * ny
a = 18  # [um] lattice spacing


## Initialize
# Geometry
L1,L2 = lattice.get_domain_sizes(nx,ny ,a)
coords, lattice_ids = lattice.get_nodes_and_ids(nx, ny, a)
N1, T1 = lattice.get_neighbours_list(coords, nx, ny, a)
get_k = lattice.define_get_k(nx, ny, a)
get_mtwist = lattice.define_get_mtwist(coords, nx, ny, a)

# Physics: sine-coupling
coupling = physics.define_sine_coupling(sin_str)
right_side_of_ODE = physics.define_right_side_of_ODE(coupling, freq, N1,T1)
solve_cycle = carpet.define_solve_cycle(right_side_of_ODE, 2 * period, carpet.get_mean_phase)


