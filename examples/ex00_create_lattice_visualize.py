import carpet
import carpet.lattice.triangular as lattice
import carpet.visualize as vis
import matplotlib.pyplot as plt

# Geometry
a = 18  # [um]
nx = 6  # number of cilia in x-direction
ny = 6  # must be even
N = nx * ny

L1, L2 = lattice.get_domain_sizes(nx, ny, a)
coords, lattice_ids = lattice.get_nodes_and_ids(nx, ny, a)  # get cilia (nodes) coordinates
NN, TT = lattice.get_neighbours_list(coords, nx, ny, a)   # get list of neighbours and relative positions
get_mtwist = lattice.define_get_mtwist(coords, nx, ny, a)  # metachronal (plane) wave

phi = get_mtwist(2, 1)  # sp.zeros([len(coords)])
vis.plot_edges(coords, TT)
vis.plot_nodes(coords, phi=phi)
plt.show()
