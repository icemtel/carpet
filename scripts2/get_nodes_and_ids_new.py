'''
- Much cleaner than previous version.
- But not using it, as I don't want to deal with compatibility issues
'''
import carpet.lattice.triangular as lattice
import numpy as np

def get_nodes_and_ids_new(nx, ny, a):
    '''
    NB: lattice_ids are different here from what defined previously
    '''
    coords = []
    lattice_ids = []
    a1 = np.array([1, 0]) * a
    a2 = np.array([0, 1]) * a * np.sqrt(3) / 2
    for i in range(nx):
        for j in range(ny):
            c = i * a1 + j * a2 + (j % 2) * a1 / 2  # extra term for triangular lattice
            coords.append(c)
            lattice_ids.append((i, j))

    return np.array(coords), np.array(lattice_ids, dtype=np.int)


def ix22D_to_ix2(ix2D):
    try:  # If array
        ix = np.array([j + i * ny for (i, j) in ix2D], dtype=np.int)
    except:  # if a tuple
        ix = ix2D[1] + ix2D[0] * ny
    return ix


def ix2_to_ix22D(ix):
    try:  # If array
        ix2D = np.array([[i // ny, i % ny] for i in ix], dtype=np.int)
    except:  # if integer
        ix2D = np.zeros(2, dtype=np.int)
        ix2D[0] = ix // ny
        ix2D[1] = ix % ny
    return ix2D


def get_window(corner_ix2, window_sizes):
    '''
    :param corner_ix2: index of the oscillator in the left bottom corner
    :param window_sizes: two integers
    :return: 1D indices in the original lattice of nodes in the window.
    '''
    if np.array(corner_ix2).size != 1:
        raise TypeError('corner ix should be an integer')
    assert window_sizes[1] % 2 == 0

    corner_ix2D = ix2_to_ix22D(corner_ix2)

    ix2D_list = []
    for i in range(window_sizes[0]):
        for j in range(window_sizes[1]):
            ix2D = corner_ix2D + (i, j)

            # If corner in an odd row: shift indices,
            # s.t. (1,0) node is to the upper-right from the corner
            ix2D[0] -= (corner_ix2D[1] % 2) * ((j + 1) % 2)  # shift

            ix2D = ix2D % (nx, ny)  # account for the finite lattice
            ix2D_list.append(ix2D)
    indices_window = ix22D_to_ix2(ix2D_list)
    return indices_window


nx,ny = 4, 32
a = 18

coords2, lattice_ids2 = get_nodes_and_ids_new(nx, ny, a)  # nodes, numbered by the new scheme
get_mtwist2 = lattice.define_get_mtwist(coords2, nx, ny, a)

print(coords2.shape)