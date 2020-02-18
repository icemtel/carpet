'''
- Create a chain of cilia.
- Define symmetry classes
- Use classes to solve ODE
- Map back to cilia
'''
import scipy as sp
import carpet
import carpet.lattice.ring1d as lattice
import carpet.classes as cc

## Parameters
# Physics
set_name = 'machemer_1' # which hydrodynamic coefficients to use
order_g11 = (8,0)
order_g12 = (4,4)
T = 31.25 # [ms] period of cilia beat
freq = 2 * sp.pi / T # [rad/ms] angular frequency

# Geometry
N = 6
a = 18  # [um] lattice spacing
e1 = (1, 0)

## Initialize
# Geometry
L1 = lattice.get_domain_size(N, a)
coords, lattice_ids = lattice.get_nodes_and_ids(N, a, e1)
N1, T1 = lattice.get_neighbours_list(N, a, e1)
get_mtwist = lattice.define_get_mtwist(coords, N, a, e1)
get_k = lattice.define_get_k(N, a, e1)

# # Physics
# gmat_glob, q_glob = lattice.define_gmat_glob_and_q_glob(set_name, a, e1, N1, T1,order_g11, order_g12, T)
# right_side_of_ODE = lattice.define_right_side_of_ODE(gmat_glob, q_glob)
# solve_cycle = carpet.define_solve_cycle(right_side_of_ODE,2 * T, carpet.get_mean_phase)

# m-twist
k1 = 1
phi_k = get_mtwist(k1)

# Map to classes
ix_to_class, class_to_ix = cc.get_classes(phi_k)
nclass = len(class_to_ix)
# Get classes representatives
# Get one cilium from each of cilia classes
unique_cilia_ids = cc.get_unique_cilia_ix(class_to_ix)# equivalent to sp.array([class_to_ix[iclass][0] for iclass in range(nclass)], dtype=sp.int64)
# Get connections
N1_class, T1_class = cc.get_neighbours_list_class(unique_cilia_ids, ix_to_class, N1, T1)

# Physics
gmat_glob_class, q_glob_class = lattice.define_gmat_glob_and_q_glob(set_name, a, e1, N1_class, T1_class,
                                                                    order_g11, order_g12, T)
right_side_of_ODE_class = lattice.define_right_side_of_ODE(gmat_glob_class, q_glob_class)
solve_cycle_class = carpet.define_solve_cycle(right_side_of_ODE_class, 2 * T, carpet.get_mean_phase)

# Solve ODE
phi_k_class = phi_k[unique_cilia_ids]
sol = solve_cycle_class(phi_k_class, tol=10**-4)

phi_class_after_cycle = sol.y.T[-1]
# Map from classes back to cilia
phi_after_cycle = phi_class_after_cycle[ix_to_class]
## Print how much phase changed
print(phi_after_cycle - phi_k - 2 * sp.pi) # [ 1.91280547e-04  1.84610272e-05 -1.01082591e-04 -1.82536506e-04  -5.43819305e-05  1.28259454e-04]