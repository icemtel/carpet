'''
(Non-periodic) chain of oscillators
'''

import numpy as np
import carpet
import carpet.lattice.ring1d as lattice

## Parameters
# Geometry
nx = 16
ny = 1  # even number
N = nx * ny
a = 18  # [um] lattice spacing
e1 = np.array([1, 0])
_, e2 = lattice.get_basis(e1)

## Initialize

# pre save coordinates
try:
    coords = np.load('coords.npy')
    lattice_ids = np.load('lattice_ids.npy')
except:
    coords, lattice_ids = lattice.get_nodes_and_ids(N, a, e1)  # get nodes coordinates
    np.save('coords.npy', coords)
    np.save('lattice_ids.npy', lattice_ids)

assert len(coords) == len(lattice_ids) == N  # check

try:
    NN = np.load('NN.npy')
    TT = np.load('TT.npy')
except:
    NN, TT = lattice.get_neighbours_list_non_periodic(N, a, e1)
    np.save('NN.npy', NN)
    np.save('TT.npy', TT)

assert len(NN) == len(TT) == N  # check

L1 = lattice.get_domain_size(N, a)
get_k = lattice.define_get_k(N, a, e1)
get_mtwist = lattice.define_get_mtwist(N)
