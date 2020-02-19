import time
import numpy as np
import carpet
import carpet.visualize as vis
import carpet.lattice.triangular as lattice
import carpet.physics.friction_pairwise as physics

a = 18  # [um]
period = 31.25  # [ms] period
nx = 3
ny = 4  # must be even
set_name = 'machemer_1'
order_g11 = (8, 0)
order_g12 = (4, 4)
tol = 10 ** -6  # solver tolerance

# Geometry
L1,L2 = lattice.get_domain_sizes(nx,ny ,a)
coords, lattice_ids = lattice.get_nodes_and_ids(nx, ny, a)  # get cilia (nodes) coordinates
N1, T1 = lattice.get_neighbours_list(coords, nx, ny, a)  # get list of neighbours and relative positions
e1, e2 = lattice.get_basis()
get_k = lattice.define_get_k(nx, ny, a)
get_mtwist = lattice.define_get_mtwist(coords, nx, ny, a)
# Physics
gmat_glob, q_glob = physics.define_gmat_glob_and_q_glob(set_name, e1, e2, a, N1, T1, order_g11, order_g12, period)
right_side_of_ODE = physics.define_right_side_of_ODE(gmat_glob, q_glob)
solve_cycle = carpet.define_solve_cycle(right_side_of_ODE, 2 * period, phi_global_func=carpet.get_mean_phase)


# Solve
phi0 = np.zeros([len(coords)])  # initial condition
phi0[5] += 3  # perturb

start = time.time()  # track wall-clock time
solution = solve_cycle(phi0, tol)
time_spent = time.time() - start
print("Time spent", time_spent)

# Get result
phi1 = solution.y.T[-1]
print("Change in phase after one cycle:", phi1 - phi0 - 2 * np.pi)
# Visualize
vis.plot_nodes(coords, phi=(phi1 - phi0 - 2 * np.pi) % (2 * np.pi))  # can't see phase difference on this scale
vis.plt.show()
