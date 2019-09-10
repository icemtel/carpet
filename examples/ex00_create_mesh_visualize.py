import scipy as sp
import carpet
import carpet.lattice_triangular as lattice
import carpet.visualize as vis
import matplotlib.pyplot as plt


a = 18 # [um]
nx = 9
ny = 6 # must be even

coords, lattice_ids = lattice.get_nodes_and_ids(nx,ny,a) # get cilia (nodes) coordinates
N1, T1 = lattice.get_neighbours_list(coords, nx,ny, a)   # get list of neighbours and relative positions

Phi = sp.zeros([len(coords)])
vis.plot_edges(coords, T1)
vis.plot_nodes(coords, phi=Phi)
plt.show()