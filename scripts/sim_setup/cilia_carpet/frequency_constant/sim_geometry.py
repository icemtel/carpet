'''

'''

import numpy as np
import carpet
import carpet.lattice.triangular as lattice

## Parameters
# Geometry
nx = 16
ny = 16  # even number
N = nx * ny
a = 18  # [um] lattice spacing

## Initialize
# Geometry
connections = [a * np.array([np.cos(psi), np.sin(psi)])  # First neighbours
               for psi in np.linspace(0, 2 * np.pi, 6, endpoint=False)] \
              + [np.sqrt(3) * a * np.array([np.cos(psi), np.sin(psi)])
                  for psi in [np.pi / 2, 3 * np.pi / 2]]  # 2nd neighbour (only 1 direction)

# pre save coordinates
try:
    coords = np.load('coords.npy')
    lattice_ids = np.load('lattice_ids.npy')
except:
    coords, lattice_ids = lattice.get_nodes_and_ids(nx, ny, a)  # get cilia (nodes) coordinates
    np.save('coords.npy', coords)
    np.save('lattice_ids.npy', lattice_ids)


assert len(coords) == len(lattice_ids)== N # check

try:
    NN = np.load('NN.npy')
    TT = np.load('TT.npy')
except:
    NN, TT = lattice.get_neighbours_list_general(coords, nx, ny, a, connections)
    np.save('NN.npy', NN)
    np.save('TT.npy', TT)

assert len(NN) == len(TT)== N # check



L1,L2 = lattice.get_domain_sizes(nx,ny ,a)
e1, e2 = lattice.get_basis()
get_k = lattice.define_get_k_fbz(nx, ny, a)
get_mtwist = lattice.define_get_mtwist(coords, nx, ny, a)
