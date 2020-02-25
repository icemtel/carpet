import numpy as np
import time
import carpet
import carpet.lattice.triangular as lattice
import carpet.visualize as vis

## Parameters

# Geometry
nx = 6
ny = 6  # even number
N = nx * ny
a = 18  # [um] lattice spacing

## Initialize
# Geometry
L1, L2 = lattice.get_domain_sizes(nx, ny, a)
coords, lattice_ids = lattice.get_nodes_and_ids(nx, ny, a)
N1, T1 = lattice.get_neighbours_list(coords, nx, ny, a)
e1, e2 = lattice.get_basis()
get_k = lattice.define_get_k(nx, ny, a)
get_mtwist = lattice.define_get_mtwist(coords, nx, ny, a)

# ====Physics====
## Sine-coupling
import carpet.physics.kuramoto as physics

period = 31.25  # [ms] period of cilia beat; freq = 2 * sp.pi / period [rad/ms]
freq = 2 * np.pi / period
sin_str = 0.003 * freq
coupling = physics.define_sine_coupling(sin_str)
right_side_of_ODE = physics.define_right_side_of_ODE(coupling, freq, N1, T1)
solve_cycle = carpet.define_solve_cycle(right_side_of_ODE, 2 * period, carpet.get_mean_phase)

# Solve
phi0 = np.zeros([len(coords)])  # initial condition
phi0[5] += 3  # perturbation
tol = 1e-6  # solver tolerance

start = time.time()  # track wall-clock time
solution = solve_cycle(phi0, tol)
time_spent = time.time() - start
print("Time spent", time_spent)

# Get result
phi1 = solution.y.T[-1]
print("Change in phase after one cycle:", phi1 - phi0 - 2 * np.pi)
# Visualize
dphi = phi1 - phi0 - 2 * np.pi
vis.plot_nodes(coords, phi=dphi, vmin=dphi.min(), vmax=dphi.max(), cmap='jet')
vis.plt.show()

## Hydrodynamic friction (cilia)
import carpet.physics.friction_pairwise as physics

set_name = 'machemer_1'  # hydrodynamic friction coefficients data set
period = 31.25  # [ms] period
freq = 2 * np.pi / period  # [rad/ms] angular frequency
order_g11 = (8, 0)  # order of Fourier expansion of friction coefficients
order_g12 = (4, 4)

gmat_glob, q_glob = physics.define_gmat_glob_and_q_glob(set_name, e1, e2, a, N1, T1, order_g11, order_g12, period)
right_side_of_ODE = physics.define_right_side_of_ODE(gmat_glob, q_glob)
solve_cycle = carpet.define_solve_cycle(right_side_of_ODE, 2 * period, phi_global_func=carpet.get_mean_phase)

# Solve
phi0 = np.zeros([len(coords)])  # initial condition
phi0[5] += 3  # perturbation
tol = 1e-6  # solver tolerance

start = time.time()  # track wall-clock time
solution = solve_cycle(phi0, tol)
time_spent = time.time() - start
print("Time spent", time_spent)

# Get result
phi1 = solution.y.T[-1]
print("Change in phase after one cycle:", phi1 - phi0 - 2 * np.pi)
# Visualize
dphi = phi1 - phi0 - 2 * np.pi
vis.plot_nodes(coords, phi=dphi, vmin=dphi.min(), vmax=dphi.max(), cmap='jet')
vis.plt.show()
