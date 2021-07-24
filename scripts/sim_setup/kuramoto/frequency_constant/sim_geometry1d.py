'''
Chain of oscillators with (optional) periodic boundary conditions (ring).
'''

import numpy as np
import carpet
import carpet.lattice.ring1d as lattice

## Parameters
# Geometry
N = 16
a = 18  # [um] lattice spacing

## Initialize
# Geometry
e1 = np.array([1,0]) # direction of the chain (for Kuramoto model - does not matter)

# pre save coordinates
try:
    coords = np.load('coords.npy')
    lattice_ids = np.load('lattice_ids.npy')
except:
    coords, lattice_ids = lattice.get_nodes_and_ids(N, a, e1)  # get oscillator positions
    np.save('coords.npy', coords)
    np.save('lattice_ids.npy', lattice_ids)


assert len(coords) == len(lattice_ids)== N # check

try:
    NN = np.load('NN.npy')
    TT = np.load('TT.npy')
except:
    NN, TT = lattice.get_neighbours_list(N, a, e1) # ring
    ## Uncomment line below for system without periodic boundary conditions
    # NN,TT = lattice.get_neighbours_list_non_periodic(N, a, e1) # chain
    np.save('NN.npy', NN)
    np.save('TT.npy', TT)

assert len(NN) == len(TT)== N # check
if len(distances) == 1:
    assert len(TT[0]) == 6 # check that there are 6 neighbours
else:
    print("WARNING: ! check sim geom -> distances")

L1 = lattice.get_domain_size(N ,a)
e1, e2 = lattice.get_basis(e1)
get_k = lattice.define_get_k(N, a, e1)
get_mtwist = lattice.define_get_mtwist(coords, N, a, e1)
