'''
Functions to load Fourier series of hydrodynamic friction coefficients from a file
'''
import os
import scipy as sp
import scipy.linalg as lin
import scipy.sparse as sparse
import pandas as pd
import math

script_path = os.path.dirname(os.path.realpath(__file__))
number_format = '.4G'


# Basis functions for one-dimensional Fourier series - used to compute 2D basis
def get_basis_function1D(j):
    if j > 0:
        return lambda x: math.sin(j * x)
    elif j == 0:
        return lambda x: 1
    elif j < 0:
        return lambda x: math.cos(j * x)  # cos(-n x) = cos(n x)


# Basis functions for two-dimensional Fourier series
def get_basis_function(j1, j2):
    a = get_basis_function1D(j1)
    b = get_basis_function1D(j2)
    return lambda x1, x2: a(x1) * b(x2) #basis1D_dict[j1](x1) * basis1D_dict[j2](x2) ! dont make one-liner, it runs slower

# Compute value of two-dimensional Fourier series, given its list of coefficients
def fourier_series2D(coeffs, coeff_ids, swap_axes=False):
    coeffs = sp.array(coeffs)
    # Save copies of basis functions - gains ~20% to speed, including parallel jobs; but not sure due to variation
    basis_functions = [get_basis_function(j1,j2) for j1, j2 in coeff_ids]

    if swap_axes:
        def series(x, y):
            b = [basis_function(y,x) for basis_function in basis_functions]
            summ = coeffs.dot(b)
            return summ
    else:
        def series(x, y):
            b = [basis_function(x,y) for basis_function in basis_functions]
            summ = coeffs.dot(b)
            return summ

    return series


# Prepare list of mode indices for two-dimensional Fourier series
def get_coeff_ids(order, sample_size=None, truncate_triangular=False):
    '''
    Get ids of the coefficients, get basis function ids up to specified order in each direction.
    :param order:   Max mode number to extract; if a tuple given - separate number for x and y direction.
                    Can be `None` - then the sample_size will be used to determine the order.
    :param sample_size: Tuple (len(xs),len(ys). If given - Sin( N_x x) won't be added to the set, if 2 * N_x = len(xs).
                        Reason: if xs are equally spaced, starting from 0,
                        then Sin(N_x x) will evaluate to zero at every sample point.
    :param truncate_triangular: If True: truncate series, s.t. n_x/N_x + n_y/N_y =< 1
    '''
    if sample_size == None:
        sample_size = (sp.inf, sp.inf)
    len_xs, len_ys = sample_size
    if order is None:
        if len_xs == sp.inf or len_ys == sp.inf:
            raise ValueError("Either order or sample_size should be specified.")
        order = len_xs // 2, len_ys // 2

    try:
        order1, order2 = order  # If order is a tuple
    except:
        order1 = order2 = order  # If order is a number

    if 2 * order1 > len_xs or 2 * order2 > len_ys:
        raise ValueError('Requested more coefficients than number of sample points')

    coeff_ids = []
    for j1 in range(-order1, order1 + 1):
        for j2 in range(-order2, order2 + 1):
            if 2 * j1 == len_xs or 2 * j2 == len_ys:  # Skip Sin(N_x x) or Sin(N_y y)
                continue
            elif truncate_triangular and abs(j1) * order2 + abs(j2) * order1 > order1 * order2:
                continue  # Skip if outside the triangleq
            else:
                coeff_ids.append((j1, j2))
    return coeff_ids


def get_translation_rotation_folder(translation, angle, format=number_format):
    translation_rotation_folder_template = 't={0[0]:{format}}_{0[1]:{format}}_{0[2]} a={1[0]:{format}}_{1[1]:{format}}'
    return translation_rotation_folder_template.format(translation, angle, format=format)


def load_coeffs_from_file(filename, order_max=None, truncate_triangular=False):
    '''
    Load Fourier series coefficients from a file.
    :param order_max: see order in get_coeff_ids()
    '''

    df = pd.read_csv(filename, dtype={'n1': sp.int32, 'n2': sp.int32, 'coeff': sp.float64})
    df.set_index(['n1', 'n2'], drop=False, inplace=True)

    if order_max is not None:
        coeffs_ids_set = get_coeff_ids(order_max, sample_size=None, truncate_triangular=truncate_triangular)
        ids_in_df = df.index.values
        # Drop those values, which are not returned by the get_coeff_ids function
        ids_to_truncate = [idx for idx in ids_in_df if
                           idx not in coeffs_ids_set]
        df.drop(index=ids_to_truncate, inplace=True)
    return df


def load_function_from_file(filename, order_max=None, truncate_triangular=False):
    '''
    Load Fourier series from a file with coefficients.
    :param order_max: see order in get_coeff_ids()
    '''
    df = load_coeffs_from_file(filename, order_max, truncate_triangular)
    coeffs = sp.array(df['coeff'])
    coeff_ids = [(m, n) for (m, n) in zip(df['n1'], df['n2'])]

    series = fourier_series2D(coeffs, coeff_ids)
    return series


###
# connections = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, 1), (1, -1)]

def load_self_friction_and_interactions(set_name, connections, e1, e2, a, order_g11, order_g12):
    '''
    Load as a function of translation?
    But then I will have troulbe searching for the right entry in the dictionary
    Load as a function of (n,m)?
    - Can work, just need to specify e1,e2,a and connection numbers
    :param set_name: e.g. 'machemer_1'
    '''
    if set_name == 'machemer_1':
        friction_coeffs_root = os.path.join(script_path, 'friction_coeffs', 'machemer_1')

    small_number = 1e-8
    interactions_dict = {}
    for n, m in connections:  # relative cilia positions in lattice space
        translation = a * e1 * n + a * e2 * m
        # Determine if translation vector is in lower half-plane - then use coeffs from upper part of the plane
        # but swap cilia indices
        if translation[1] < 0 or (abs(translation[1]) < small_number and translation[0] < 0):
            in_lower_halfplane = True
            translation = - translation
        else:
            in_lower_halfplane = False

        # Replace very small coordinates with zero - to correctly compile the folder name
        for i, ti in enumerate(translation):
            if abs(ti) < small_number:
                translation[i] = 0

        # Get file name and load g_ij as a function
        translation3D = (*translation, 0.)
        translation_folder = get_translation_rotation_folder(translation3D, (0, 0))
        filename = os.path.join(friction_coeffs_root, translation_folder, 'g12_ft.dat')

        # Load friction as a function (Fourier sum); swap order of input when needed
        g_ij = load_function_from_file(filename, order_max=order_g12, truncate_triangular=False)

        df = load_coeffs_from_file(filename, order_max=order_g12, truncate_triangular=False)
        coeffs = sp.array(df['coeff'])
        coeff_ids = [(m, n) for (m, n) in zip(df['n1'], df['n2'])]

        # swap variable in fouerier series if in lower half plane - fix of symmetry
        interactions_dict[(n, m)] = fourier_series2D(coeffs, coeff_ids, swap_axes=in_lower_halfplane)

    ##  Load self-friction
    self_friction_dict = {}
    for n, m in connections:
        translation = a * e1 * n + a * e2 * m
        # Determine if translation vector is in lower half-plane - then use coeffs from upper part of the plane
        # but change indexes 1 and 2
        if translation[1] < 0 or (translation[1] == 0 and translation[0] < 0):
            in_lower_halfplane = True
            translation = - translation
        else:
            in_lower_halfplane = False

        # Replace very small coordinates with zero - to correctly compile the folder name
        for i, ti in enumerate(translation):
            if abs(ti) < 10 ** -8:
                translation[i] = 0

        translation3D = (*translation, 0.)
        translation_folder = get_translation_rotation_folder(translation3D, (0, 0))

        if in_lower_halfplane:  # fix symmetry
            filename = os.path.join(friction_coeffs_root, translation_folder, 'g22_ft.dat')  # TODO: Check
            order = order_g11[::-1]  # reverse order
        else:
            filename = os.path.join(friction_coeffs_root, translation_folder, 'g11_ft.dat')
            order = order_g11

        # Load friction as a function (fourier sum); swap order of input when needed
        df = load_coeffs_from_file(filename, order_max=order, truncate_triangular=False)
        coeffs = sp.array(df['coeff'])
        coeff_ids = [(m, n) for (m, n) in zip(df['n1'], df['n2'])]

        # swap variable in fouerier series if in lower half plane - fix of symmetry
        self_friction_dict[(n, m)] = fourier_series2D(coeffs, coeff_ids, swap_axes=in_lower_halfplane)

    return self_friction_dict, interactions_dict


# Simple self-friction


###
def define_gmat_glob_and_q_glob0(set_name, connections, e1, e2, a, neighbours_indices, neighbours_rel_positions,
                                 order_g11, order_g12, period):
    '''
    Define friction matrix and active driving force

    Sparse matrix representation (csr)
    https://docs.scipy.org/docs/scipy/reference/generated/scipy.sparse.csr_matrix.html
    https://docs.scipy.org/docs/scipy-0.14.0/reference/sparse.linalg.html

    :param neighbours_indices: list of list: first list has indices of the first cilium, second - neighbours of second cilium, etc.
    :param neighbours_rel_positions: the same structure as `neighbours_indices`; contains relative positions of cilia (taking into account periodicity)
    '''

    transformation_matrix = lin.inv((sp.array([e1, e2]).T)) / a  # from Euclidean to lattice

    def euclidean_to_lattice_coords(translation):
        nm = transformation_matrix @ translation
        return nm

    self_friction_dict, interactions_dict = load_self_friction_and_interactions(set_name, connections,
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
            nm_int = sp.around(nm)

            if not sp.allclose(nm, nm_int, rtol=10 ** -3, atol=10 ** -3):  # Assure correct rounding
                raise ValueError

            gij = interactions_dict[(nm_int[0], nm_int[1])]
            rows.append(i)
            cols.append(j)
            gij_func_list.append(gij)

    rows = sp.array(rows)
    cols = sp.array(cols)
    N = len(neighbours_indices)  # matrix size

    def gmat_glob(phi):
        '''
        Compute friction matrix component; return a sparse matrix object.
        '''
        vals = []
        for i, j, gij_func in zip(rows, cols, gij_func_list):
            vals.append(gij_func(phi[i], phi[j]))
        return sparse.csr_matrix((vals, (rows, cols)), shape=(N, N))

    freq = 2 * sp.pi / period

    ## Define Q(phi)
    ## First - calibrate only with self-friction
    def q_glob(phi):
        '''
        :param phi: vector of phases, corresponding to each of the cilia
        '''
        return sp.array([freq * gii(phii, 0) for phii in phi])

    return gmat_glob, q_glob


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    set_name = "machemer_1"
    connections = [(-1, 0), (1, 0)]

    a = 18
    e1 = sp.array([1, 0])
    e2 = sp.array([0, 1])
    order_g11 = (8, 0)
    order_g12 = (4, 4)
    T = 31.25
    N1 = [sp.array([2, 1]), sp.array([0, 2]), sp.array([1, 0])]
    T1 = [[sp.array([18, 0]), sp.array([-18, 0])],
          [sp.array([18, 0]), sp.array([-18, 0])],
          [sp.array([18, 0]), sp.array([-18, 0])]]

    gmat_glob, q_glob = define_gmat_glob_and_q_glob0(set_name, connections, e1, e2, a, order_g11, order_g12, N1, T1, T)

    phis = sp.array([[0, 0, 0], [0.5, 0.5, 0.5], [sp.pi, sp.pi, sp.pi],
                     [3 / 2 * sp.pi, 3 / 2 * sp.pi, 3 / 2 * sp.pi], [2 * sp.pi, 2 * sp.pi, 2 * sp.pi]])

    plt.plot(phis.T[0], 1 / 2 / sp.pi * T * sp.array([q_glob(phi) for phi in phis]))
    plt.show()

## Checks
# # Ben checked: mtwist solution -> classes of cilia with different inital phase have the same phase after N cycles
# # Is self-friction different if cilia are located differently? Check
# # Result: difference ~0.2%
# phis = sp.linspace(0, 2 * sp.pi, 50, endpoint=False)
# test_vals = []
# for key, g_ii in self_friction_dict.items():
#     gii_vals = [g_ii(phi,0) for phi in phis]
#     test_vals.append(gii_vals)

# test_vals = sp.array(test_vals)
# mean_vals = sp.mean(test_vals, axis=0)
# for i,vals in enumerate(test_vals):
#     delta = (vals - mean_vals ) / mean_vals
#     plt.plot(phis, delta)
# plt.show()

# # Test 2: g_ii(phi,phi); compare with mean_vals obtained earlier
# test_vals = []
# for key, g_ii in self_friction_dict.items():
#     gii_vals = [g_ii(phi,phi) for phi in phis]
#     test_vals.append(gii_vals)

# test_vals = sp.array(test_vals)
# for i,vals in enumerate(test_vals):
#     delta = (vals - mean_vals ) / mean_vals
#     plt.plot(phis, delta)
# plt.show()

# # If different - Check one by one
# for key, g_ii in self_friction_dict.items():
#     vals = [g_ii(phi,0) for phi in phis]
#     plt.plot(phis, vals)
# plt.show()

# for key, g_ii in self_friction_dict.items():
#     vals = [g_ii(phi,phi) for phi in phis]
#     plt.plot(phis, vals)
# plt.show()
