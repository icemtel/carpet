'''
Triangular lattice, but rotated by 90 degrees. Because of rectangular domain, it is not trivial to introduce lattice for
any rotation.
'''
import math
import numpy as np
from scipy.linalg import norm

def get_cell_sizes(a):
    cell_length = a * 3 ** (1 / 2) / 2
    cell_height = a
    return cell_length, cell_height


def get_basis():
    e1 = np.array([3 ** (1 / 2) / 2, 0.5])
    e2 = np.array([0, 1])
    return e1, e2


def get_domain_sizes(nx, ny, a):
    cell_length, cell_height = get_cell_sizes(a)
    L1 = cell_length * nx
    L2 = cell_height * ny
    return L1, L2


def get_nodes_and_ids(nx, ny, a):
    '''
    Create a triangular lattice in rectangular domain
    :param nx,ny: - Number of cilia in each direction; ny must be even
    :param a: cilia spacing
    '''
    # unit vectors
    # spatial dimension of simulation domain

    e1, e2 = get_basis()

    cell_length, cell_height = get_cell_sizes(a)
    L1, L2 = get_domain_sizes(nx, ny, a)
    eps = 10 ** -4 * a  # better safe than sorry: small number to account for rounding errors below

    # generate triangular lattice
    coords = []  # list of position vectors for honeycomb lattice
    lattice_ids = []  # list of corresponding lattice indices (n,m)

    for n in range(-2 * nx, 2 * nx):
        for m in range(- 2 * ny, 2 * ny):
            x = n * a * e1 + m * a * e2  # position vector
            # position vector within bounds?
            if (x[0] >= 0 - eps) and (x[1] >= 0 - eps) and (x[0] < L1 - cell_length + eps) and (x[1] < L2 - eps):
                coords.append(x)
                lattice_ids.append((n, m))

    coords = np.array(coords)
    lattice_ids = np.array(lattice_ids)

    return coords, lattice_ids


def get_neighbours_list(coords, nx, ny, a):
    '''
    Find nodes which are located at the distance `a`
    :return: list of neighbours, list of relative neighbour positions
    '''
    eps = 10 ** -4 * a

    L1, L2 = get_domain_sizes(nx, ny, a)

    N = len(coords)

    ## find nearest neigbors
    def get_neighbours(i, j):
        '''
        If the distance between two points is equal to the lattice spacing, return vector connecting them, else None.
        Takes into account lattice periodicity
        OK: sign
        '''
        for a1 in range(-1, 2):
            for a2 in range(-1, 2):
                translation = coords[j, :] - coords[i, :] + [a1 * L1, 0] + [0, a2 * L2]
                if a - eps < norm(translation) < a + eps:
                    return translation
        return None

    ## TO change in case of rectangular domain
    N1 = [[] for _ in coords]  # list of lists of neighbours
    T1 = [[] for _ in coords]  # list of translation vectors between neighbours
    # loop over pairs of lattice points
    for i in range(N):
        for j in range(i + 1, N):
            # lattice points close?
            # (account for periodic boundary conditions)
            translation = get_neighbours(i, j)
            if translation is not None:
                N1[i].append(j)
                T1[i].append(translation)

                N1[j].append(i)
                T1[j].append(- translation)
    return N1, T1


### mtwist solutions ###

def get_dual_basis(a):
    '''
    Reciprocal lattice for rectangular unit cell
    Ben's method
    '''
    ax, ay = get_cell_sizes(a)
    a1 = ax * np.array([1, 0]) / a
    a2 = ay * np.array([0, 1]) / a
    R = np.array([[0, 1], [-1, 0]])  # rotation by 90deg
    a1dual = 2 * np.pi * (R @ a2) / (a1 @ (R @ a2)) / a
    a2dual = 2 * np.pi * (R @ a1) / (a2 @ (R @ a1)) / a
    return a1dual, a2dual  # [1/L]


def get_dual_basis2(a):
    '''
    Does the same as the first one, but maybe it's more clear what the method does.
    Google:
    B = (a1, a2)
    D = (a1dual, a2dual)
    -> D = inv(B.T)
    Below implemented the same thing, but with a factor of 2pi
    '''
    from scipy.linalg import inv, solve
    ax, ay = get_cell_sizes(a)
    a1 = ax * np.array([1, 0])
    a2 = ay * np.array([0, 1])
    # B = np.array([a1 ,a2]).T
    BT = np.array([a1, a2])
    # By definition D = inv(B.T)
    # Miss a step of transposing
    D = 2 * np.pi *  inv(BT)

    return D[:, 0], D[:, 1]


def define_get_k_naive(nx, ny, a):
    a1dual, a2dual = get_dual_basis(a)

    def get_k(k1, k2):  # get wave vector corresponding to wave numbers
        k = k1 * a1dual / nx + k2 * a2dual / ny
        return k

    return get_k


def define_get_k(nx, ny, a):
    '''
    Checked: get_k is equivalent to get_k_naive: gives the same mtwists mod 2pi
    The same as in `lattice_triangular`, but indices swapped, and a1dual swapped with a2dual.
    '''
    a1dual, a2dual = get_dual_basis(a)

    def get_k(k1, k2):  # get wave vector corresponding to wave numbers
        k = k1 * a1dual / nx + k2 * a2dual / ny
        if k[1] >= a2dual[1] / 2:
            k[1] -= a2dual[1]
            k[0] -= a1dual[0] / 2
        if k[0] >= a1dual[0] / 2:
            k[0] -= a1dual[0]
        return k

    return get_k


def define_get_mtwist(coords, nx, ny, a):
    get_k = define_get_k(nx, ny, a)

    def mod(x):
        '''
        fmod(x,y) is not equivalent to (x % y): https://docs.python.org/3/library/math.html and
        is preferred when working with floats
        :return: a value in interval from 0 to 2pi
        '''
        x = math.fmod(x, 2 * np.pi)
        if x < 0:
            x += 2 * np.pi
        return x

    # Fill mtwist array
    mtwist_phi = np.zeros((nx, ny, nx * ny))

    for k1 in range(nx):
        for k2 in range(ny):
            # wave vector
            k = get_k(k1, k2)  # k1 * a1dual / nx + k2 * a2dual / ny
            for ix in range(nx * ny):
                mtwist_phi[k1, k2, ix] = mod(- np.dot(k, coords[ix, :]))

    def get_mtwist(k1, k2):
        return np.array(mtwist_phi[k1, k2])

    return get_mtwist


### Friction coefficients - ver 1. ###
def get_connections():
    '''
    :return: Relative positions of neighbouring cilia in lattice coordinates
             First order neighbours
    '''
    return [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, 1), (1, -1)]



# Cilia coupling. Must be removed in the future - call module with coupling


def define_gmat_glob_and_q_glob(set_name, a, neighbours_indices, neighbours_rel_positions,
                                order_g11, order_g12, T):
    '''
    Shortcut for cilia coupling - keep for backwards compatibility
    :param set_name: e.g. 'machemer_1'
    :param a: lattice spacing
    :param neighbours_indices: list of neighbours
    :param neighbours_rel_positions: list of relative neighbour positions
    :param order_g11: Fourier expansion order
    :param order_g12:
    :param T:
    :return: gmat, q_glob  - functions
    '''
    from carpet.physics.friction_pairwise import define_gmat_glob_and_q_glob as define_gmat_glob_and_q_glob0

    import warnings

    warnings.warn("To be depricated! Import 'coupling' instead", DeprecationWarning)

    connections = get_connections()
    e1, e2 = get_basis()
    return define_gmat_glob_and_q_glob0(set_name, connections, e1, e2, a,
                                       neighbours_indices, neighbours_rel_positions,
                                       order_g11, order_g12, T)


def define_right_side_of_ODE(gmat_glob, q_glob):
    import carpet.physics.friction_pairwise as coupling
    import warnings

    warnings.warn("To be depricated! Import 'coupling' instead", DeprecationWarning)

    return coupling.define_right_side_of_ODE(gmat_glob, q_glob)


if __name__ == '__main__':
    # OK: nodes plot
    # OK: cell sizes: the same as in lattice_triangular, but inversed components
    # OK: dual basis: the same as in lattice_triangular, but inversed components
    # OK: mtwists

    import matplotlib.pyplot as plt
    import carpet.visualize as visualize

    a = 18
    nx = 8
    ny = 12  # must be even

    coords, lattice_ids = get_nodes_and_ids(nx, ny, a)
    N1, T1 = get_neighbours_list(coords, nx, ny, a)

    ## Visualize
    # visualize.plot_edges(coords, T1)
    # visualize.plot_nodes(coords)
    # visualize.plt.show()

    ### Check mtwists
    get_mtwist_phi = define_get_mtwist(coords, nx, ny, a)
    # Check: all m-twist solutions are indeed different from each other
    small_number = 1e-8
    for k1 in range(0, nx):
        for k2 in range(0, ny):
            for m1 in range(0, nx):
                for m2 in range(0, ny):
                    if max(abs(get_mtwist_phi(k1, k2) - get_mtwist_phi(m1, m2))) < small_number:
                        assert (k1 == m1) and (k2 == m2)
    print("OK: m-twists checked")

    # Visualize mtwist
    # Visualization
    # ... nodes
    # wave numbers
    k1 = 3  # [] integer
    k2 = 1
    phi = get_mtwist_phi(k1, k2)

    # fig, ax = visualize.plot_nodes(coords, phi=phi, colorbar=False)
    ## Duplicate
    # L1, L2 = get_domain_sizes(nx, ny, a)
    # visualize.plot_nodes(coords + np.array([L1, 0])[np.newaxis, :], phi=phi, colorbar=False)
    # visualize.plot_nodes(coords + np.array([0, L2])[np.newaxis, :], phi=phi, colorbar=True)
    # visualize.plot_node_numbers(coords, a)
    #
    # ax.set_title('m-twist: (' + str(k1) + ',' + str(k2) + ')')
    # ax.get_xaxis().set_visible(False)
    # ax.get_yaxis().set_visible(False)
    # fig.set_size_inches(6, 8)
    # plt.show()

    ## Friction
    # order_g11 = (8, 0)
    # order_g12 = (4, 4)
    # T = 31.25
    # gmat, qglob = define_gmat_glob_and_q_glob('machemer_1', a, N1, T1, order_g11, order_g12, T)

    ### Check get_k vs get_k_naive:
    # Result: they are equivalent
    get_k = define_get_k(nx, ny, a)
    get_k_naive = define_get_k_naive(nx, ny, a)

    for k1 in range(nx):
        for k2 in range(ny):
            print(k1, k2)

            k = get_k(k1, k2)
            k_naive = get_k_naive(k1, k2)

            for coord in coords:
                if abs(np.exp(1j * k @ coord) - np.exp(1j * k_naive @ coord)) > 10 ** -8:
                    print('WHOOPS')
                # print("whoops", k, k_naive)

    ### Dual basis test?
    print(get_cell_sizes(a))
    print(get_dual_basis(a))
    print(get_dual_basis2(a))
