'''
- Generate rectangular lattice
- In rectangular domain, assuming periodic boundary conditions.

MAYBE: rewrite get_neighbour_list - there's much easier way for the rectangular lattice
'''
import math
import numpy as np
from scipy.linalg import norm

from carpet.various import get_basis_dual, mod2pi


def get_basis():
    e1 = np.array([1, 0])
    e2 = np.array([0, 1])
    return e1, e2


def get_cell_sizes(a):
    cell_length = a
    cell_height = a
    return cell_length, cell_height


def get_domain_sizes(nx, ny, a):
    cell_length, cell_height = get_cell_sizes(a)
    L1 = cell_length * nx
    L2 = cell_height * ny
    return L1, L2


def get_nodes_and_ids(nx, ny, a):
    '''
    '''
    e1, e2 = get_basis()
    coords = []
    lattice_ids = []
    a1 = e1 * a
    a2 = e2 * a
    for i in range(nx):
        for j in range(ny):
            c = i * a1 + j * a2  # extra term for triangular lattice
            coords.append(c)
            lattice_ids.append((i, j))

    return np.array(coords), np.array(lattice_ids, dtype=np.int)


def get_neighbours_list(coords, nx, ny, a, distances=(1,)):
    '''
    For each node looks for other nodes at specified distances (multiplied by lattice edge length `a`).
    Those nodes are saved in `N1` list. Relative positions are saved in `T1` list.
    :distances: list of expected distances to neighbours (normalized by a)
                Examples: 1-neighbours: [1]
                          2-neighbours*: [1, 3 ** 0.5]
                          2-neighbours:  [1, 3 ** 0.5, 2]
                Assumption: d * a < max(L1,L2)
    :return: list of neighbours, list of relative neighbour positions
    '''
    if nx == 2 or ny == 2:
        import warnings
        warnings.warn("nx=2 or ny=2 => wrong number of neighbours (5 or 4)\n"
                      "some oscillators were supposed to be connected twice, but this is not implemented")

    eps = 10 ** -4 * a
    L1, L2 = get_domain_sizes(nx, ny, a)
    if max(distances) * a >= max([L1, L2]):
        raise NotImplementedError("Assumption: d * a < max(L1,L2) is not satisfied")

    N = len(coords)

    ## find nearest neighbors
    def get_neighbours(i, j):
        '''
        If the distance between two points is equal to the lattice spacing, return vector connecting them, else None.
        Takes into account lattice periodicity
        OK: sign
        '''
        for a1 in range(-1, 2):
            for a2 in range(-1, 2):
                translation = coords[j, :] - coords[i, :] + [a1 * L1, 0] + [0, a2 * L2]
                for d in distances:
                    if d * a - eps < norm(translation) < d * a + eps:
                        return translation
        return None

    N1 = [[] for _ in coords]  # list of lists of neighbours indices
    T1 = [[] for _ in coords]  # list of lists of translation vectors between neighbours
    # loop over pairs of lattice points
    for i in range(N):
        for j in range(i + 1, N):
            translation = get_neighbours(i, j)  # check if neighbours
            if translation is not None:  # is yes - add to the list
                N1[i].append(j)
                T1[i].append(translation)

                N1[j].append(i)  # save some iterations by using that being neighbors is symmetrical relation
                T1[j].append(- translation)
    return N1, T1


### Wave vectors and reciprocal lattice ###

def get_basis_dual_domain(nx, ny, a):
    '''
    Reciprocal vectors for rectangular domain
    '''
    e1, e2 = get_basis()
    d1, d2 = get_cell_sizes(a)
    a1 = nx * d1 * e1
    a2 = ny * d2 * e2
    b1, b2 = get_basis_dual(a1, a2)
    return b1, b2  # [rad/L]


def get_basis_dual_cell(a):
    '''
    Reciprocal vectors for the triangular unit
    '''
    e1, e2 = get_basis()
    a1 = a * e1
    a2 = a * e2
    b1, b2 = get_basis_dual(a1, a2)
    return b1, b2  # [rad/L]


def define_get_k_naive(nx, ny, a):
    '''
    The simplest way to get a wave vector from the dual basis.
    Other functions shift wave vector k to a different unit cell of reciprocal lattice.
    :return: wave vector k [rad/L]
    '''
    a1dual, a2dual = get_basis_dual_domain(nx, ny, a)

    def get_k(k1, k2):  # get wave vector corresponding to wave numbers
        k = k1 * a1dual + k2 * a2dual
        return k

    return get_k


#
# def define_shift_k_to_fbz(a):
#     '''
#     Defines a function, which get any wave vector, and shifts it into the first Brillouin zone
#     :param a:
#     :return:
#     '''
#
#     def project(vec, basis_vec):
#         basis_vec = np.asarray(basis_vec)
#         return vec @ basis_vec / (basis_vec @ basis_vec)
#
#     def decompose_recip(k):
#         # decompose a vector to vectors b1,b2,b3 (double of normals of the hexagon cell)
#         ms = np.array([project(k, bi) for bi in bs])
#         return ms
#
#     b1, b2 = get_basis_dual_cell(a)
#     b3 = b1 + b2
#     bs = [b1, b2, b3]
#
#     def k_to_fbz(k, eps=1e-8):
#         k = np.array(k)
#         num_iter = 0
#         ms = decompose_recip(k)
#
#         while np.amax(abs(ms)) > 0.5 + eps and num_iter < 10:  # +eps to acccount for numerical edge case
#             i = int(np.argmax(abs(ms)))  # start with the direction with the biggest projection
#             mi = ms[i]
#             bi = bs[i]
#             k -= bi * np.round(mi)  # shift by integer value
#             ms = decompose_recip(k)
#             num_iter += 1
#         if num_iter == 10:
#             raise ValueError("Didn't converge to a unit cell - check algorithm!")
#         return k
#
#     return k_to_fbz


def define_get_k_fbz(nx, ny, a):
    assert ny % 2 == 0  # check that ny is even
    get_k_naive = define_get_k_naive(nx, ny, a)
    b1, b2 = get_basis_dual_cell(a)
    size1 = b1[0]
    size2 = b2[1]

    def get_k(k1, k2):
        k = get_k_naive(k1, k2)
        k[0] = (k[0] + size1 / 2) % size1 - size1 / 2
        k[1] = (k[1] + size2 / 2) % size2 - size2 / 2
        return k

    return get_k


def define_get_k_fbz_all(nx, ny, a):
    assert ny % 2 == 0  # check that ny is even
    get_k_fbz = define_get_k_fbz(nx, ny, a)

    b1, b2 = get_basis_dual_cell(a)
    b3 = b1 + b2
    bs = [b1, b2, b3]

    def get_k_all(k1, k2, eps=1e-8):
        """
        Return a list of all possible representations of wave vector with wave numbers (k1,k2)
        in the first Brillouin zone: 1 - if inside, 2 - if on the edge, 3 - in vertex
        """
        k = get_k_fbz(k1, k2)
        ks = [k]
        knorm = norm(k)
        # Shift k in all lattice directions; if it still whithin FBZ (<=> equal norm)
        # => add it to the list
        for b in bs:
            for sign in [-1, 1]:
                if abs(norm(k + sign * b) - knorm) < eps:
                    ks.append(k + sign * b)

        return ks

    return get_k_all


def define_get_mtwist(coords, nx, ny, a):
    get_k = define_get_k_naive(nx, ny, a)

    # Fill mtwist array
    mtwist_phi = np.zeros((nx, ny, nx * ny))

    for k1 in range(nx):
        for k2 in range(ny):
            # wave vector
            k = get_k(k1, k2)  # k1 * a1dual / nx + k2 * a2dual / ny
            for ix in range(nx * ny):
                mtwist_phi[k1, k2, ix] = mod2pi(- np.dot(k, coords[ix, :]))

    def get_mtwist(k1, k2):
        return np.array(mtwist_phi[k1, k2])

    return get_mtwist


if __name__ == '__main__':
    # OK: rectangular lattice
    # OK: neighbours and translations
    # OK: get_k_fbz and m-twist: two versions of get_mtwist based on different get_k match
    # OK: get_k_fbz_all - correct length of lists

    import matplotlib.pyplot as plt
    import carpet.visualize as vis

    nx = 6
    ny = 6
    a = 10
    coords, lattice_ids = get_nodes_and_ids(nx, ny, a)

    # N1, T1 = get_neighbours_list(coords, nx, ny, a, distances=[1])
    # print(N1, T1)
    # print("Neighbours as array shape:", np.array(N1).shape)
    # ## Visualize
    # vis.plot_edges(coords, T1)
    # vis.plot_nodes(coords)
    # vis.plt.show()

    print(get_basis_dual_domain(nx, ny, a))
    print(get_basis_dual_cell(a))

    get_k_naive = define_get_k_naive(nx, ny, a)
    get_k = define_get_k_fbz(nx, ny, a)

    # k1,k2 = 3,0
    # k_naive = get_k_naive(k1,k2)
    # k = get_k(k1,k2)
    # plt.scatter(*k_naive, color='blue')
    # plt.scatter(*k, color='red')
    # plt.gca().set_aspect('equal')
    # plt.show()
    #
    # for k1 in range(nx):
    #     for k2 in range(ny):
    #         k_naive = get_k_naive(k1, k2)
    #         k = get_k(k1, k2)
    #         print(k_naive)
    #
    #         plt.scatter(*k_naive, alpha=0.5, color='blue')
    #         plt.scatter(*k, color='red', alpha=0.5)
    # plt.gca().set_aspect('equal')
    # plt.show()
    #
    #
    # # test get_k_fbz
    # def define_get_mtwist2(coords, nx, ny, a):
    #     get_k = define_get_k_fbz(nx, ny, a)
    #
    #     # Fill mtwist array
    #     mtwist_phi = np.zeros((nx, ny, nx * ny))
    #
    #     for k1 in range(nx):
    #         for k2 in range(ny):
    #             # wave vector
    #             k = get_k(k1, k2)  # k1 * a1dual / nx + k2 * a2dual / ny
    #             for ix in range(nx * ny):
    #                 mtwist_phi[k1, k2, ix] = mod2pi(- np.dot(k, coords[ix, :]))
    #
    #     def get_mtwist(k1, k2):
    #         return np.array(mtwist_phi[k1, k2])
    #
    #     return get_mtwist
    #
    #
    # get_mtwist = define_get_mtwist(coords, nx, ny, a)
    # get_mtwist2 = define_get_mtwist2(coords, nx, ny, a)
    #
    # for k1 in range(nx):
    #     for k2 in range(ny):
    #         assert np.allclose(get_mtwist(k1,k2), get_mtwist2(k1,k2))
    #
    # Test get_k_fbz_all

    get_k_all_fbz = define_get_k_fbz_all(nx, ny, a)
    k1,k2 = 3, 0
    print(len(get_k_all_fbz(k1,k2)))

    k1, k2 = 3,3
    print(len(get_k_all_fbz(k1, k2)))