import time
import scipy as sp
import carpet
import carpet.visualize as vis
import carpet.lattice.triangular as lattice
import carpet.physics.friction_pairwise as physics

a = 18  # [um]
T = 31.25  # [ms] period
nx = 3
ny = 4  # must be even
set_name = 'machemer_1'
order_g11 = (8, 0)
order_g12 = (4, 4)
tol = 10 ** -6  # solver tolerance

# Geometry
coords, lattice_ids = lattice.get_nodes_and_ids(nx, ny, a)  # get cilia (nodes) coordinates
N1, T1 = lattice.get_neighbours_list(coords, nx, ny, a)  # get list of neighbours and relative positions
connections = lattice.get_connections()
e1, e2 = lattice.get_basis()

# Physics
gmat_glob, q_glob = physics.define_gmat_glob_and_q_glob(set_name, connections, e1, e2, a, N1, T1, order_g11, order_g12,
                                                        T)
right_side_of_ODE = physics.define_right_side_of_ODE(gmat_glob, q_glob)

mean_phase = lambda phi: sp.mean(phi)

solve_cycle = carpet.define_solve_cycle(right_side_of_ODE, 2 * T, mean_phase)

# Solve
phi0 = sp.zeros([len(coords)])  # initial condition
phi0[5] += 3  # perturb

start = time.time()  # track wall-clock time
solution = solve_cycle(phi0, tol)
time_spent = time.time() - start
print("Time spent", time_spent)

# Get result
phi1 = solution.y.T[-1]
print("Change in phase after one cycle:", phi1 - phi0 - 2 * sp.pi)
# Visualize
vis.plot_nodes(coords, phi=(phi1 - phi0 - 2 * sp.pi) % (2 * sp.pi))  # can't see phase difference on this scale
vis.plt.show()
