"""
The example with nx=3, ny=4 should run ~1.5secs
"""
import numpy as np
import carpet.visualize as vis
import matplotlib.pyplot as plt

#=====Ring=====
import carpet.lattice.ring1d as lattice

# Geometry
N = 6   # number of cilia
a = 18  # [um] lattice spacing
e1 = (1, 0) # direction of the chain

## Initialize
# Geometry
L1 = lattice.get_domain_size(N, a)
coords, lattice_ids = lattice.get_nodes_and_ids(N, a, e1)  # get cilia (nodes) coordinates
N1, T1 = lattice.get_neighbours_list(N, a, e1) # get list of neighbours and relative positions
e1, e2 = lattice.get_basis(e1)
get_k = lattice.define_get_k(N, a, e1)
get_mtwist = lattice.define_get_mtwist(coords, N, a, e1)


phi = get_mtwist(2, 0)  # sp.zeros([len(coords)])
vis.plot_edges(coords, T1)
vis.plot_nodes(coords, phi=phi)
plt.ylim([-L1/10, L1/10])
plt.show()

#=====Chain=====
## Almost the same, as ring - only end nodes have only 1 neighbour each.

import carpet.lattice.ring1d as lattice

# Geometry
N = 6   # number of cilia
a = 18  # [um] lattice spacing
e1 = (1, 0) # direction of the chain

## Initialize
# Geometry
L1 = lattice.get_domain_size(N, a)
coords, lattice_ids = lattice.get_nodes_and_ids(N, a, e1)  # get cilia (nodes) coordinates
N1, T1 = lattice.get_neighbours_list_non_periodic(N, a, e1) # get list of neighbours and relative positions
e1, e2 = lattice.get_basis(e1)
get_k = lattice.define_get_k(N, a, e1)
get_mtwist = lattice.define_get_mtwist(coords, N, a, e1)


phi = get_mtwist(2, 0)  # sp.zeros([len(coords)])
vis.plot_edges(coords, T1)
vis.plot_nodes(coords, phi=phi)
plt.ylim([-L1/10, L1/10])
plt.show()



#=====Lattice 1=====
import carpet.lattice.triangular as lattice

# Geometry
a = 18  # [um]
nx = 3  # number of cilia in x-direction
ny = 4  # in y-direction:  must be even
N = nx * ny

## Initialize
# Geometry
L1,L2 = lattice.get_domain_sizes(nx,ny ,a)
coords, lattice_ids = lattice.get_nodes_and_ids(nx, ny, a)  # get cilia (nodes) coordinates
N1, T1 = lattice.get_neighbours_list(coords, nx, ny, a)  # get list of neighbours and relative positions
e1, e2 = lattice.get_basis()
get_k = lattice.define_get_k_fbz(nx, ny, a)
get_mtwist = lattice.define_get_mtwist(coords, nx, ny, a)

phi = get_mtwist(2, 0)  # sp.zeros([len(coords)])
vis.plot_edges(coords, T1)
vis.plot_nodes(coords, phi=phi)
plt.show()


#====Lattice 2=====
import carpet.lattice.triangular2 as lattice

# Geometry
a = 18  # [um]
nx = 4  # must be even
ny = 3  # number of cilia in y-direction
N = nx * ny

## Initialize
# Geometry
L1,L2 = lattice.get_domain_sizes(nx,ny ,a)
coords, lattice_ids = lattice.get_nodes_and_ids(nx, ny, a)  # get cilia (nodes) coordinates
N1, T1 = lattice.get_neighbours_list(coords, nx, ny, a)  # get list of neighbours and relative positions
e1, e2 = lattice.get_basis()
get_k = lattice.define_get_k(nx, ny, a)
get_mtwist = lattice.define_get_mtwist(coords, nx, ny, a)


phi = get_mtwist(2, 0)  # sp.zeros([len(coords)])
vis.plot_edges(coords, T1)
vis.plot_nodes(coords, phi=phi)
plt.show()
