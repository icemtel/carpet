"""
The example with nx=3, ny=4 should run ~1.5secs
"""
import time
import numpy as np
import carpet
import carpet.visualize as vis
import carpet.lattice.triangular as lattice
import carpet.physics.friction_pairwise as physics

## Parameters
# Physics
set_name = 'machemer_1'  # hydrodynamic friction coefficients data set
period = 31.25  # [ms] period
freq = 2 * np.pi / period  # [rad/ms] angular frequency
order_g11 = (8, 0)  # order of Fourier expansion of friction coefficients
order_g12 = (4, 4)
# Geometry
a = 18  # [um]
nx = 3  # number of cilia in x-direction
ny = 4  # must be even
N = nx * ny

## Initialize
# Geometry
L1, L2 = lattice.get_domain_sizes(nx, ny, a)
coords, lattice_ids = lattice.get_nodes_and_ids(nx, ny, a)  # get cilia (nodes) coordinates
N1, T1 = lattice.get_neighbours_list(coords, nx, ny, a)  # get list of neighbours and relative positions
e1, e2 = lattice.get_basis()
get_k = lattice.define_get_k_fbz(nx, ny, a)
get_mtwist = lattice.define_get_mtwist(coords, nx, ny, a)

# Physics
gmat_glob, q_glob = physics.define_gmat_glob_and_q_glob(set_name, e1, e2, a, N1, T1, order_g11, order_g12, period,
                                                        use_numba=False)
right_side_of_ODE = physics.define_right_side_of_ODE(gmat_glob, q_glob)
solve_cycle = carpet.define_solve_cycle(right_side_of_ODE, 2 * period, phi_global_func=carpet.get_mean_phase)

## Solve
phi0 = np.zeros([len(coords)])  # initial condition
phi0[5] += 3  # perturbation
tol = 10 ** -6  # solver tolerance

start = time.time()  # track wall-clock time
solution = solve_cycle(phi0, tol)
time_spent = time.time() - start
print("Time spent", time_spent)

## Get result
phi1 = solution.y.T[-1]
print("Change in phase after one cycle:", phi1 - phi0 - 2 * np.pi)
# Visualize
dphi = phi1 - phi0 - 2 * np.pi
vis.plot_nodes(coords, phi=dphi, vmin=dphi.min(), vmax=dphi.max(), cmap='jet')
vis.plt.show()
