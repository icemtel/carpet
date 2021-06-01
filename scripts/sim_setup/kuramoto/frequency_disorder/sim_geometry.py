'''
2021-03-12: separate geometry and physics in order to have 1 geometry file for different frequqncies sets
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

distances = [1]

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
    NN, TT = lattice.get_neighbours_list(coords, nx, ny, a, distances)
    np.save('NN.npy', NN)
    np.save('TT.npy', TT)

assert len(NN) == len(TT)== N # check
if len(distances) == 1:
    assert len(TT[0]) == 6 # check that there are 6 neighbours
else:
    print("WARNING: ! check sim geom -> distances")

L1,L2 = lattice.get_domain_sizes(nx,ny ,a)
e1, e2 = lattice.get_basis()
get_k = lattice.define_get_k_fbz(nx, ny, a)
get_mtwist = lattice.define_get_mtwist(coords, nx, ny, a)
