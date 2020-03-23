'''
Define coupling and driving force based on hydrodynamic simulations of cilia.

- Sum up pair-wise interactions
- Assume that self-friction doesn't depend on phases of neighbouring cilia.
- This version is about x17 times faster than the old one, if max(order_g11) = max(order_g12), otherwise about x12 times.

Can possibly change function `define_gmat_glob_and_q_glob` s.t. it doesn't need e1,e2,a:
- Load what translations are present in dataset + add opposite vectors, exploiting symmetry
- Compare them with translations in neighbours_rel_positions
- Save intersection of two. If some are present in neighburs_rel_positions, but not in data - raise Error
- Use list index instead of `connections` in dictionary index
(Not doing it as right now it's not worth time investment); will only be useful if
I start studying networks rather than lattices

MAYBE: move load_gii_coeff_matrix, load_gij_coeff_matrix to friction
'''
import os
import numpy as np
from scipy import linalg as lin, sparse as sparse
from scipy.sparse.linalg import spsolve
from carpet.friction import get_basis_function1D, \
    get_translation_rotation_folder, load_coeffs_from_file, get_friction_coeffs_path


def define_right_side_of_ODE(gmat_glob, q_glob, linear_solver=spsolve):
    '''
    \Gamma . \Phi' = Q => \Phi' = right_side_of_ODE = \Gamma^-1 \Q
    :param gmat_glob: function of phase vector which returns friction matrix \Gamma
    :param q_glob: function of phase vector which returns active driving force Q
    :param linear_solver: A x = b => computes x.
           Recommended: scipy.sparse.linalg.solve for sparse 2D-array \Gamma or scipy.linalg.solve for dense matrix
    :return: instantaneous phase frequency \Phi'
    '''

    def right_side_of_ODE(t, phi):
        '''
        Which method to use?
        One test on sparse and iterative sparse methods gave the same result ~270s of computation
        (maybe solving the linear equation is not the bottleneck?)
        '''
        ## Dense
        # phidot = lin.inv(gmat_glob(phi)) @ q_glob(phi)
        ## Sparse
        phidot = linear_solver(gmat_glob(phi), q_glob(phi))
        ## Sparse - iterative
        # phidot, info = splin.bicgstab(gmat_glob(phi), q_glob(phi), tol=tol) #Check if this tolerance is sufficient
        return phidot

    return right_side_of_ODE


def load_gii_coeff_matrix(friction_coeffs_root, translation, order_g11, eps=1e-8):
    """
    Load Fourier expansion coefficients of friction coefficients gii and gij.
    Represent result as a matrix (numpy.array)
    :param translation: means relative position of the second cilium
    """

    # Determine if translation vector is in lower half-plane - then use coeffs from upper part of the plane
    # but change indexes 1 and 2
    if translation[1] < 0 or (abs(translation[1]) < eps and translation[0] < 0):
        in_lower_halfplane = True
        translation = - translation
    else:
        in_lower_halfplane = False

    # Replace very small coordinates with zero - to correctly compile the folder name
    for i, ti in enumerate(translation):
        if abs(ti) < eps:
            translation[i] = 0

    translation3D = (*translation, 0.)
    translation_folder = get_translation_rotation_folder(translation3D, (0, 0))

    if in_lower_halfplane:  # fix symmetry
        filename = os.path.join(friction_coeffs_root, translation_folder, 'g22_ft.dat')
        order = order_g11[::-1]  # reverse order
    else:
        filename = os.path.join(friction_coeffs_root, translation_folder, 'g11_ft.dat')
        order = order_g11

    # Load friction as a function (fourier sum); swap order of input when needed
    df = load_coeffs_from_file(filename, order_max=order, truncate_triangular=False)
    coeffs = np.array(df['coeff'])
    coeff_ids = [(m, n) for (m, n) in zip(df['n1'], df['n2'])]

    order = max(order_g11)  # determine matrix size
    coeff_mat = np.zeros([2 * order + 1, 2 * order + 1])
    for (i, j), coeff in zip(coeff_ids, coeffs):
        imat = i + order
        jmat = j + order
        coeff_mat[imat, jmat] = coeff

    # If in lower half-plane, we used data from the upper half-plane;
    # swap arguments (transpose the fourier coeffs matrix)
    if in_lower_halfplane:
        coeff_mat = coeff_mat.T

    return coeff_mat


def load_gij_coeff_matrix(friction_coeffs_root, translation, order_g12, eps=1e-8):
    """
    Load Fourier expansion coefficients of friction coefficients gii and gij.
    Represent result as a matrix (numpy.array)
    :param translation: means relative position of the second cilium
    """

    # Determine if translation vector is in lower half-plane - then use coeffs from upper part of the plane
    # but swap cilia indices
    if translation[1] < 0 or (abs(translation[1]) < eps and translation[0] < 0):
        in_lower_halfplane = True
        translation = - translation
    else:
        in_lower_halfplane = False

    # Replace very small coordinates with zero - to correctly compile the folder name
    for i, ti in enumerate(translation):
        if abs(ti) < eps:
            translation[i] = 0

    # Get file name and load g_ij as a function
    translation3D = (*translation, 0.)
    translation_folder = get_translation_rotation_folder(translation3D, (0, 0))

    if in_lower_halfplane:  # fix symmetry
        filename = os.path.join(friction_coeffs_root, translation_folder, 'g21_ft.dat')  # should be the same as g12
        order = order_g12[::-1]  # reverse order
    else:
        filename = os.path.join(friction_coeffs_root, translation_folder, 'g12_ft.dat')
        order = order_g12

    # Load friction as a function (Fourier sum); swap order of input when needed
    df = load_coeffs_from_file(filename, order_max=order, truncate_triangular=False)
    coeffs = np.array(df['coeff'])
    coeff_ids = [(m, n) for (m, n) in zip(df['n1'], df['n2'])]

    order = max(order_g12)  # determine matrix size
    coeff_mat = np.zeros([2 * order + 1, 2 * order + 1])
    for (i, j), coeff in zip(coeff_ids, coeffs):
        imat = i + order
        jmat = j + order
        coeff_mat[imat, jmat] = coeff

    # If in lower half-plane, we used data from the upper half-plane;
    # swap arguments (transpose the fourier coeffs matrix)
    if in_lower_halfplane:
        coeff_mat = coeff_mat.T

    return coeff_mat


def load_self_friction_coeffs_dict(set_name, connections, e1, e2, a, order_g11, eps=1e-8):
    '''
    A shortcut to load friction coeffs for a number of different relative positions of cilia.
    :param set_name: e.g. 'machemer_1'
    :param connections: [(n1,m1),(n2,m2),..] - relative positions of the second cilia in lattice coordinates - integer numbers
    :param e1,e2: directors of the lattice - unit vectors
    :param a: lattice spacing
    :return: a dictionary. keys: (n,m), values:  coeffs matrices gii or gij; size depends on order_g11 and order_g12
    '''
    friction_coeffs_root = get_friction_coeffs_path(set_name)

    ##  Load self-friction
    self_friction_dict = {}
    for n, m in connections:
        translation = a * e1 * n + a * e2 * m
        self_friction_dict[(n, m)] = load_gii_coeff_matrix(friction_coeffs_root, translation, order_g11, eps)

    return self_friction_dict


def load_interactions_coeffs_dict(set_name, connections, e1, e2, a, order_g12, eps=1e-8):
    '''
    A shortcut to load friction coeffs for a number of different relative positions of cilia.
    :param set_name: e.g. 'machemer_1'
    :param connections: [(n1,m1),(n2,m2),..] - relative positions of the second cilia in lattice coordinates - integer numbers
    :param e1,e2: directors of the lattice - unit vectors
    :param a: lattice spacing
    :return: a dictionary. keys: (n,m), values:  coeffs matrices gii or gij; size depends on order_g11 and order_g12
    '''
    friction_coeffs_root = get_friction_coeffs_path(set_name)

    interactions_dict = {}
    for n, m in connections:  # relative cilia positions in lattice space
        translation = a * e1 * n + a * e2 * m
        interactions_dict[(n, m)] = load_gij_coeff_matrix(friction_coeffs_root, translation, order_g12, eps)

    return interactions_dict


def define_get_basis_matrix(order):
    basis_functions = [get_basis_function1D(j) for j in range(-order, order + 1)]

    def get_basis_matrix(phi):
        '''
        This is a transpose of Alternant matrix - https://en.wikipedia.org/wiki/Alternant_matrix
        Matrix which entries are each of (2*order+1) basis functions, evaluated at each of N phases.
        :param phi: vector of N phases
        :return: Represented by list of length N with elements - numpy.array of size (2 * order +1)
        '''
        CS = [np.array([b(phij) for b in basis_functions], dtype=np.float) for phij in phi]
        return CS

    return get_basis_matrix


def define_gmat_glob_and_q_glob(set_name, e1, e2, a, neighbours_indices, neighbours_rel_positions,
                                order_g11, order_g12, period, eps=1e-8):
    '''
    Define friction matrix and active driving force.

    - Collect matrices Kij, representing Fourier expansion coefficients of gmat entries (gij).
        gij = b(phi_i) @ Kij @ b(phi_j), (*)
        where b(x) - vector of Fourier basis functions (cos(m x), 1, sin(m x)), evaluated at x
     - Define a function get_AT, which calculates a matrix with columns b(phi_i), i=1..N
     - Then function gmat_glob would calculate gij as (*),
       calculating AT matrix and getting Kij from a pre-saved list K3D.
       This proved to be rather faster, than previous implementation.
     - Output is a *sparse* matrix. See

      Sparse matrix representation (csr)
      https://docs.scipy.org/docs/scipy/reference/generated/scipy.sparse.csr_matrix.html
      https://docs.scipy.org/docs/scipy-0.14.0/reference/sparse.linalg.html

    :param connections: neighbours positions in lattice coordinates
    :param e1,e2: lattice basis vectors
    :param neighbours_indices: list of list: first list has indices of the first cilium, second - neighbours of second cilium, etc.
    :param neighbours_rel_positions: the same structure as `neighbours_indices`; contains relative positions of cilia (taking into account periodicity)
    :param period: period of single cilium beat - used to calibrate active driving force
    :param eps: small number; used to check that cilia positions are at lattice nodes
    '''
    N = len(neighbours_indices)  # number of cilia
    # Define a helper functions - transforms from Euclidean to lattice coordinates
    transformation_matrix = lin.inv((np.array([e1, e2]).T)) / a  # from Euclidean to lattice

    def euclidean_to_lattice_coords(translation):
        nm = transformation_matrix @ translation
        return nm

    # Determine which translations are present - this might be slow, but we usually call this function only once
    # Store them in lattice coordinates (decomposition in form n * a * e1 + m * a * e2)
    connections = set()
    for i, (neighbours, translations) in enumerate(zip(neighbours_indices, neighbours_rel_positions)):
        for j, t in zip(neighbours, translations):
            nm = euclidean_to_lattice_coords(t)  # Detetermine relative positions of cilia in lattice
            nm_int = np.around(nm)
            if not np.allclose(nm, nm_int, rtol=10 ** -3, atol=10 ** -3):  # Assure correct rounding
                raise ValueError

            nm_int = tuple(np.array(nm_int, dtype=np.int16))  # for hashing; and keeping int
            connections.add(nm_int)

    gii_coeff_dict = load_self_friction_coeffs_dict(set_name, connections, e1, e2, a, order_g11, eps)
    gii_coeff_mat = gii_coeff_dict[(1, 0)]
    gij_coeff_dict = load_interactions_coeffs_dict(set_name, connections, e1, e2, a, order_g12, eps)

    K3D = []
    for i in range(N):
        K3D.append(gii_coeff_mat)  # put .copy() if want to edit entries!

    for i, (neighbours, translations) in enumerate(zip(neighbours_indices, neighbours_rel_positions)):
        for j, t in zip(neighbours, translations):
            nm = euclidean_to_lattice_coords(t)  # get lattice coordinates
            nm_int = np.around(nm)
            gij_coeff_mat = gij_coeff_dict[tuple(nm_int)]
            K3D.append(gij_coeff_mat)

    rows = []
    cols = []
    for i in range(N):
        rows.append(i)
        cols.append(i)

    for i, (neighbours, translations) in enumerate(zip(neighbours_indices, neighbours_rel_positions)):
        for j, t in zip(neighbours, translations):
            rows.append(i)
            cols.append(j)

    order1 = max(order_g11)
    order2 = max(order_g12)

    get_AT = define_get_basis_matrix(order1)  # alternation matrix transposed

    if order1 == order2:
        def gmat_glob(phi):
            AT = get_AT(phi)

            vals = []
            ## Collect non-zero entries of gij matrix
            for (i, j, Kij) in zip(rows, cols, K3D):
                vals.append(AT[i].dot(Kij.dot(AT[j])))
            # Equivalent, but slower: np.einsum('kl,k,l',Kij,AT[:,i],AT[:,j]) # AT[:,i] @ Kij @ AT[:,j]
            vals = np.array(vals, dtype=np.float)
            return sparse.csr_matrix((vals, (rows, cols)), shape=(N, N))


    else:  # If orders are different, have to calculate separate matrices for gii and gij
        get_AT2 = define_get_basis_matrix(order2)

        def gmat_glob(phi):
            vals = []
            ## First N elements - diagonal gii
            AT = get_AT(phi)

            for (i, j, Kij) in zip(rows, cols, K3D[:N]):
                vals.append(AT[i].dot(Kij.dot(AT[j])))

            ## Then gij
            AT = get_AT2(phi)
            for (i, j, Kij) in zip(rows[N:], cols[N:], K3D[N:]):
                vals.append(AT[i].dot(Kij.dot(AT[j])))
            # slower: np.einsum('kl,k,l',Kij,CS[:,i],CS[:,j]) # CS[:,i] @ Kij @ CS[:,j]
            vals = np.array(vals, dtype=np.float)

            return sparse.csr_matrix((vals, (rows, cols)), shape=(N, N))

    ## Define Q(phi)
    ## In this version - calibrate only with self-friction
    ## Could optimize this and calculate right_side_of_ODE straight away
    ## But having q_glob and gmat decoupled will allow easier changes in the future
    ## And the speed gain is not big (not tested)
    freq = 2 * np.pi / period

    def get_gii(phi1, phi2):
        AT = get_AT([phi1, phi2])
        return AT[0].dot(gii_coeff_mat.dot(AT[1]))

    def q_glob(phi):
        '''
        :param phi: vector of phases, corresponding to each of the cilia
        '''
        return np.array([freq * get_gii(phii, 0) for phii in phi])

    return gmat_glob, q_glob


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import carpet.physics.friction_pairwise_v1 as physics_test

    set_name = "machemer_1"
    connections = [(-1, 0), (1, 0)]

    a = 18
    e1 = np.array([1, 0])
    e2 = np.array([0, 1])
    order_g11 = (4, 0)
    order_g12 = (4, 4)
    T = 31.25
    N1 = [np.array([2, 1]), np.array([0, 2]), np.array([1, 0])]
    T1 = [[np.array([18, 0]), np.array([-18, 0])],
          [np.array([18, 0]), np.array([-18, 0])],
          [np.array([18, 0]), np.array([-18, 0])]]

    gmat_glob, q_glob = define_gmat_glob_and_q_glob(set_name, e1, e2, a, N1, T1, order_g11, order_g12, T)
    gmat_glob_test, q_glob_test = physics_test.define_gmat_glob_and_q_glob(set_name, e1, e2, a, N1, T1, order_g11,
                                                                           order_g12, T)

    # phis = np.array([[0, 0, 0], [0.5, 0.5, 0.5], [np.pi, np.pi, np.pi],
    #                  [3 / 2 * np.pi, 3 / 2 * np.pi, 3 / 2 * np.pi], [2 * np.pi, 2 * np.pi, 2 * np.pi]])
    # plt.plot(phis.T[0], 1 / 2 / np.pi * T * np.array([q_glob(phi) for phi in phis]))
    # plt.show()

    for _ in range(10):
        phi = 2 * np.pi * np.random.rand(3)

        assert np.allclose(gmat_glob(phi).todense(), gmat_glob_test(phi).todense(), rtol=1e-12, atol=1e-12)
        assert np.allclose(q_glob(phi), q_glob_test(phi), rtol=1e-12, atol=1e-12)

    print('pass')
