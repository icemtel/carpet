'''
- Generate triangular lattice
- Define friction matrix, etc.
'''
import math
import scipy as sp
from scipy.linalg import norm
import carpet.friction as friction


def get_domain_size(N, a):
    return N * a


def get_nodes_and_ids(num_cilia, a, direction):
    direction = sp.array(direction) / norm(direction)
    translation = a * direction

    coords = []
    lattice_ids = []
    position = sp.array([0., 0])
    for ix in range(num_cilia):
        coords.append(sp.array(position))
        lattice_ids.append((ix, 0))
        position += translation

    return sp.array(coords), sp.array(lattice_ids)


def get_neighbours_list(num_cilia, a, direction):
    direction = sp.array(direction) / norm(direction)
    translation = a * direction

    N1 = []
    T1 = []
    for ix in range(num_cilia):
        neighbours = sp.array([(ix - 1) % num_cilia, (ix + 1) % num_cilia], dtype=sp.int64)
        translations = [translation, - translation]
        N1.append(neighbours)
        T1.append(translations)
    return N1, T1


def get_neighbours_list_non_periodic(num_cilia, a, direction):
    direction = sp.array(direction) / norm(direction)
    translation = a * direction

    translation_neighbours = [-translation, translation]
    N1 = []
    T1 = []
    for ix in range(num_cilia):
        neighbours = []
        translations = []
        for ix_neighbour, translation_neighbour in zip([(ix - 1), (ix + 1)], translation_neighbours):
            if 0 <= ix_neighbour < num_cilia:
                neighbours.append(ix_neighbour)
                translations.append(translation_neighbour)
        N1.append(neighbours)
        T1.append(translations)
    return N1, T1


### mtwist solutions ###
#
# def define_get_k_naive(N, a, direction):
#     L = get_domain_size(N, a)
#
#     def get_k_naive(k1):  # get wave vector corresponding to wave numbers
#         return k1 * 2 * sp.pi / L * direction
#
#     return get_k_naive


def define_get_k(N, a, direction):
    direction = sp.array(direction) / norm(direction)
    L = get_domain_size(N, a)

    def get_k(k1, k2=0):  # get wave vector corresponding to wave numbers
        if k1 >= N // 2 + 1:
            return (k1 - N) * 2 * sp.pi / L * direction

        else:
            return k1 * 2 * sp.pi / L * direction

    return get_k


def define_get_mtwist(coords, N, a, direction):
    get_k = define_get_k(N, a, direction)

    def mod(x):
        '''
        fmod(x,y) is not equivalent to (x % y): https://docs.python.org/3/library/math.html and
        is preferred when working with floats
        :return: a value in interval from 0 to 2pi
        '''
        x = math.fmod(x, 2 * sp.pi)
        if x < 0:
            x += 2 * sp.pi
        return x

    # Fill mtwist array
    mtwist_phi = sp.zeros((N, 1, N * 1))

    for k1 in range(N):
        for k2 in range(1):
            # wave vector
            k = get_k(k1)  # k1 * a1dual / nx + k2 * a2dual / ny
            for ix in range(N * 1):
                mtwist_phi[k1, k2, ix] = mod(- sp.dot(k, coords[ix, :]))

    def get_mtwist(k1, k2=0):  # two arguments for compatibility
        return sp.array(mtwist_phi[k1, k2])

    return get_mtwist


### Friction coefficients - ver 1. ###
def _get_connections():
    '''
    :return: Relative positions of neighbouring cilia in lattice coordinates
    '''
    return [(-1, 0), (1, 0)]


def define_gmat_glob_and_q_glob(set_name, a, direction, neighbours_indices, neighbours_rel_positions,
                                order_g11, order_g12, T):
    '''
    :param set_name: e.g. 'machemer_1'
    :param a: lattice spacing
    :param neighbours_indices: list of neighbours
    :param neighbours_rel_positions: list of relative neighbour positions
    :param order_g11: Fourier expansion order
    :param order_g12:
    :param T:
    :return: gmat, q_glob  - functions
    '''
    connections = _get_connections()
    e1 = sp.array(direction) / norm(direction)
    e2 = sp.array([[0,-1],[1, 0]]) @ e1
    return friction.define_gmat_glob_and_q_glob0(set_name, connections, e1, e2, a,
                                                 neighbours_indices, neighbours_rel_positions,
                                                 order_g11, order_g12, T)


import scipy.sparse.linalg as splin


def define_right_side_of_ODE(gmat_glob, q_glob):
    def right_side_of_ODE(t, phi):
        '''
        Which method to use?
        One test on sparse and iterative sparse methods gave the same result ~270s of computation
        (maybe solving the linear equation is not the bottleneck?)
        '''
        ## Dense
        # phidot = lin.inv(gmat_glob(phi)) @ q_glob(phi)
        ## Sparse
        phidot = splin.spsolve(gmat_glob(phi), q_glob(phi))
        ## Sparse - iterative
        # phidot, info = splin.bicgstab(gmat_glob(phi), q_glob(phi), tol=tol) #Check if this tolerance is sufficient
        return phidot

    return right_side_of_ODE


# TODO: Optimize:
#  - dont calculate gii twice in gmat and qglob


if __name__ == '__main__':
    # OK: chain lattice
    # OK: translations lengths (plot)
    # OK: translations direction (plot)
    # REDO:
    # OK: mtwist no duplicates
    # OK: mtwist (plot)
    # OK: get_k (since mtwist is OK)

    import matplotlib.pyplot as plt
    import carpet.visualize as visualize

    N = 10
    a = 14
    e1 = (1, 0)
    L1 = get_domain_size(N, a)
    coords, lattice_ids = get_nodes_and_ids(N, a, e1)
    N1, T1 = get_neighbours_list(N, a, e1)  # get_neighbours_list_chain_non_periodic(N, a * e1) #

    # Passes the check
    # eps = 0.001 * a
    # #If neighbours are not close: raise error
    # for ix, (neighbours, translations) in enumerate(zip(N1,T1)):
    #     for ii in range(2):
    #         delta = norm(coords[ix] - coords[neighbours[ii]] - translations[ii])
    #         if  not(delta < eps  or abs(delta - L1) < eps):
    #             print(ix, ii, delta)
    #             print(coords[ix], coords[neighbours[ii]], translations[ii])
    #             assert False
    #
    # ## Visualize
    #
    # # visualization
    # # ... plot edges
    # for i in range(N):
    #     neighbours = N1[i]
    #
    #     if len(neighbours) != 2:  # For debugging
    #         print(i, len(neighbours))
    #
    #     for j in neighbours:
    #         if (norm(coords[i, :] - coords[j, :]) < a + eps):
    #             plt.plot(coords[[i, j], 0], coords[[i, j], 1], 'b')
    #         else:
    #             # plot edges that "cut" the boundaries of the simulation domain in different color
    #             plt.plot(coords[[i, j], 0], coords[[i, j], 1], 'y:')
    #
    # # TEST: ... plot edges based on translations list
    # # OK: translations
    # # OK: translation direction
    # for i, (neighbours, translations) in enumerate(zip(N1, T1)):
    #     for j, translation in zip(neighbours, translations):
    #         points_to_plot = sp.array([coords[i], coords[i] + translation])
    #         angle = sp.arctan2(translation[1], translation[0])
    #         if abs(angle - sp.pi / 3) < eps:
    #             code = 'g-.'
    #         elif abs(angle + sp.pi * 2 / 3) < eps:
    #             code = 'r:'
    #         else:
    #             code = 'b'
    #         plt.plot(points_to_plot[:, 0], points_to_plot[:, 1], code)
    #
    # # ... plot nodes
    # plt.plot(coords[:, 0], coords[:, 1], 'r.', markersize=12)
    # plt.title('Cilia positions')
    # plt.ylim([-L1 / 3, L1 / 3])
    # plt.gca().set_aspect('equal')
    # plt.show()

    ### Check mtwist

    # Visualization
    # ... nodes
    # wave numbers
    get_mtwist_phi = define_get_mtwist(coords, N, a, e1)
    k1 = 9  # [] integer
    Phi = get_mtwist_phi(k1)

    plt.ylim([-L1 / 3, L1 / 3])
    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    plt.gcf().set_size_inches(8, 6)
    visualize.plot_nodes(coords, Phi)  # Doesn't work after update to matplotlib 3.1.0 ?!
    plt.show()
