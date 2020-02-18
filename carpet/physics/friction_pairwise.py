'''
Define coupling and driving force based on hydrodynamic simulations of cilia.

- Sum up pair-wise interactions
- Assume that self-friction doesn't depend on phases of neighbouring cilia.
'''

import numpy as np
from scipy import linalg as lin, sparse as sparse
from scipy.sparse.linalg import spsolve
import carpet.friction as friction


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
        phidot = spsolve(gmat_glob(phi), q_glob(phi))
        ## Sparse - iterative
        # phidot, info = splin.bicgstab(gmat_glob(phi), q_glob(phi), tol=tol) #Check if this tolerance is sufficient
        return phidot

    return right_side_of_ODE


def define_gmat_glob_and_q_glob(set_name, connections, e1, e2, a, neighbours_indices, neighbours_rel_positions,
                                order_g11, order_g12, period):
    '''
    Define friction matrix and active driving force

    Sparse matrix representation (csr)
    https://docs.scipy.org/docs/scipy/reference/generated/scipy.sparse.csr_matrix.html
    https://docs.scipy.org/docs/scipy-0.14.0/reference/sparse.linalg.html

    :param connections: neighbours positions in lattice coordinates
    :param e1,e2: lattice basis vectors
    :param neighbours_indices: list of list: first list has indices of the first cilium, second - neighbours of second cilium, etc.
    :param neighbours_rel_positions: the same structure as `neighbours_indices`; contains relative positions of cilia (taking into account periodicity)
    MAYBE: Optimize:
          - dont calculate gii twice in gmat and q_glob
    '''

    transformation_matrix = lin.inv((np.array([e1, e2]).T)) / a  # from Euclidean to lattice

    def euclidean_to_lattice_coords(translation):
        nm = transformation_matrix @ translation
        return nm

    self_friction_dict, interactions_dict = friction.load_self_friction_and_interactions(set_name, connections,
                                                                                e1, e2, a, order_g11, order_g12)

    gii = self_friction_dict[(1, 0)]

    rows = []
    cols = []
    gij_func_list = []

    for i, (neighbours, translations) in enumerate(zip(neighbours_indices, neighbours_rel_positions)):
        # Add self-friction
        rows.append(i)
        cols.append(i)
        gij_func_list.append(gii)
        # Add interactions
        for j, t in zip(neighbours, translations):
            nm = euclidean_to_lattice_coords(t)  # Detetermine relative positions of cilia in lattice
            nm_int = np.around(nm)

            if not np.allclose(nm, nm_int, rtol=10 ** -3, atol=10 ** -3):  # Assure correct rounding
                raise ValueError

            gij = interactions_dict[(nm_int[0], nm_int[1])]
            rows.append(i)
            cols.append(j)
            gij_func_list.append(gij)

    rows = np.array(rows)
    cols = np.array(cols)
    N = len(neighbours_indices)  # matrix size

    def gmat_glob(phi):
        '''
        Compute friction matrix component; return a sparse matrix object.
        '''
        vals = []
        for i, j, gij_func in zip(rows, cols, gij_func_list):
            vals.append(gij_func(phi[i], phi[j]))
        return sparse.csr_matrix((vals, (rows, cols)), shape=(N, N))

    freq = 2 * np.pi / period

    ## Define Q(phi)
    ## First - calibrate only with self-friction
    def q_glob(phi):
        '''
        :param phi: vector of phases, corresponding to each of the cilia
        '''
        return np.array([freq * gii(phii, 0) for phii in phi])

    return gmat_glob, q_glob




if __name__ == "__main__":
    import matplotlib.pyplot as plt

    set_name = "machemer_1"
    connections = [(-1, 0), (1, 0)]

    a = 18
    e1 = np.array([1, 0])
    e2 = np.array([0, 1])
    order_g11 = (8, 0)
    order_g12 = (4, 4)
    T = 31.25
    N1 = [np.array([2, 1]), np.array([0, 2]), np.array([1, 0])]
    T1 = [[np.array([18, 0]), np.array([-18, 0])],
          [np.array([18, 0]), np.array([-18, 0])],
          [np.array([18, 0]), np.array([-18, 0])]]

    gmat_glob, q_glob = define_gmat_glob_and_q_glob(set_name, connections, e1, e2, a, N1, T1, order_g11, order_g12, T)

    phis = np.array([[0, 0, 0], [0.5, 0.5, 0.5], [np.pi, np.pi, np.pi],
                     [3 / 2 * np.pi, 3 / 2 * np.pi, 3 / 2 * np.pi], [2 * np.pi, 2 * np.pi, 2 * np.pi]])

    plt.plot(phis.T[0], 1 / 2 / np.pi * T * np.array([q_glob(phi) for phi in phis]))
    plt.show()