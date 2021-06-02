'''
Cilia classes are used to compute fixed points faster.

- Assume symmetry like in an m-twist (make a plot to see it)
- Assume that symmetries is not broken in time -> define classes of symmetry and interactions between them.

Done:
- Create a ring of cilia.
- Define symmetry classes
- Use classes to solve ODE
- Map back to cilia
'''
import numpy as np
import carpet
import carpet.lattice.ring1d as lattice
import carpet.physics.friction_pairwise as physics
import carpet.classes as cc
import carpet.visualize as vis
import matplotlib.pyplot as plt

## Parameters
# Physics
set_name = 'machemer_1'  # hydrodynamic friction coefficients data set
period = 31.25  # [ms] period
freq = 2 * np.pi / period  # [rad/ms] angular frequency
order_g11 = (4, 0)  # order of Fourier expansion of friction coefficients
order_g12 = (4, 4)
# Geometry
N = 6  # number of cilia
a = 18  # [um] lattice spacing
e1 = (1, 0)  # direction of the chain

## Initialize
# Geometry
L1 = lattice.get_domain_size(N, a)
coords, lattice_ids = lattice.get_nodes_and_ids(N, a, e1)  # get cilia (nodes) coordinates
NN, TT = lattice.get_neighbours_list(N, a, e1)  # get list of neighbours and relative positions
e1, e2 = lattice.get_basis(e1)
get_k = lattice.define_get_k(N, a, e1)
get_mtwist = lattice.define_get_mtwist(coords, N, a, e1)
# Physics
gmat_glob, q_glob = physics.define_gmat_glob_and_q_glob(set_name, e1, e2, a, NN, TT, order_g11, order_g12, period)
right_side_of_ODE = physics.define_right_side_of_ODE(gmat_glob, q_glob)
solve_cycle = carpet.define_solve_cycle(right_side_of_ODE, 2 * period, phi_global_func=carpet.get_mean_phase)

# k-twist
k1 = 2
phi0 = get_mtwist(k1)
vis.plot_nodes(coords, phi=phi0)  # visualize!
plt.ylim([-L1 / 10, L1 / 10])
plt.show()

## Solve regularly
tol = 1e-4
sol = solve_cycle(phi0, tol)
phi1 = sol.y.T[-1] - 2 * np.pi  # after one cycle

## Now solve with classes
# Map to classes
ix_to_class, class_to_ix = cc.get_classes(phi0)
nclass = len(class_to_ix)
# Get classes representatives
# Get one oscillator from each of cilia classes
unique_cilia_ids = cc.get_unique_cilia_ix(
    class_to_ix)  # equivalent to sp.array([class_to_ix[iclass][0] for iclass in range(nclass)], dtype=sp.int64)
# Get connections
N1_class, T1_class = cc.get_neighbours_list_class(unique_cilia_ids, ix_to_class, NN, TT)

# Define physics
gmat_glob_class, q_glob_class = physics.define_gmat_glob_and_q_glob(set_name, e1, e2, a, N1_class, T1_class,
                                                                    order_g11, order_g12, period)
right_side_of_ODE_class = physics.define_right_side_of_ODE(gmat_glob_class, q_glob_class)
solve_cycle_class = carpet.define_solve_cycle(right_side_of_ODE_class, 2 * period, carpet.get_mean_phase)

# Solve ODE
phi0_class = phi0[unique_cilia_ids]
sol = solve_cycle_class(phi0_class, tol)

phi1_class = sol.y.T[-1] - 2 * np.pi
# Map from classes back to cilia
phi1_mapped_from_class = phi1_class[ix_to_class]
## Print how much phase changed
print(phi1_mapped_from_class - phi1)  # difference between two - should be on the order of tolerance or smaller
