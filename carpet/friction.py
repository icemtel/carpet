'''
Functions to load Fourier series of hydrodynamic friction coefficients from a file

- Almost not used in `physics.friction_pairwise.py` anymore

MAYBE:
- rewrite define_right_side_of_ODE - s.t. gii is not calculated twice
  in Gamma_glob and Q_glob - implemented, but the speed remained equal
- Cubic interpolation of gij
  Does give speed boost 2-10 times. (slower in parallel jobs, faster in sequential)
  But using it might give numerical problems: we might need higher order derivatives for Runge-Kutta-4(5) method
'''
import os
import numpy as np
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
    return lambda x1, x2: a(x1) * b(x2)
    # basis1D_dict[j1](x1) * basis1D_dict[j2](x2) -> shorter, but runs slower


# Compute value of two-dimensional Fourier series, given its list of coefficients
def fourier_series2D(coeffs, coeff_ids, swap_axes=False):
    coeffs = np.array(coeffs)
    # Save copies of basis functions - gains ~20% to speed, including parallel jobs; but not sure due to variation
    basis_functions = [get_basis_function(j1, j2) for j1, j2 in coeff_ids]

    if swap_axes:
        def series(x, y):
            b = [basis_function(y, x) for basis_function in basis_functions]
            summ = coeffs.dot(b)
            return summ
    else:
        def series(x, y):
            b = [basis_function(x, y) for basis_function in basis_functions]
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
        sample_size = (np.inf, np.inf)
    len_xs, len_ys = sample_size
    if order is None:
        if len_xs == np.inf or len_ys == np.inf:
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

    df = pd.read_csv(filename, dtype={'n1': np.int32, 'n2': np.int32, 'coeff': np.float64})
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
    coeffs = np.array(df['coeff'])
    coeff_ids = [(m, n) for (m, n) in zip(df['n1'], df['n2'])]

    series = fourier_series2D(coeffs, coeff_ids)
    return series


def get_friction_coeffs_path(set_name):
    if set_name in ['machemer_1', 'machemer_2', 'machemer_3']:
        friction_coeffs_root = os.path.join(script_path, '../data/friction_coeffs', set_name)
    else:
        raise KeyError("Unknown hydr. friction coefficients set name")
    return friction_coeffs_root


def load_gii(friction_coeffs_root, translation, order_g11, eps=1e-8):
    """
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

    # swap variable in fouerier series if in lower half plane - fix of symmetry
    return fourier_series2D(coeffs, coeff_ids, swap_axes=in_lower_halfplane)


def load_gij(friction_coeffs_root, translation, order_g12, eps=1e-8):
    """
    translation means relative position of the second cilium
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

    # swap variable in fouerier series if in lower half plane - fix of symmetry
    return fourier_series2D(coeffs, coeff_ids, swap_axes=in_lower_halfplane)


def load_self_friction_and_interactions(set_name, connections, e1, e2, a, order_g11, order_g12, eps=1e-8):
    '''
    :param set_name: e.g. 'machemer_1'
    :param connections: [(n1,m1),(n2,m2),..] - relative positions of the second cilia in lattice coordinates - integer numbers
    :param e1,e2: directors of the lattice - unit vectors
    :param a: lattice spacing
    :return: a tuple of dictionaries. keys: (n,m), values: function gii or gij of two arguments - phases of two cilia
    '''
    friction_coeffs_root = get_friction_coeffs_path(set_name)

    ##  Load self-friction
    self_friction_dict = {}
    for n, m in connections:
        translation = a * e1 * n + a * e2 * m
        self_friction_dict[(n, m)] = load_gii(friction_coeffs_root, translation, order_g11, eps)

    # Load interactions
    interactions_dict = {}
    for n, m in connections:  # relative cilia positions in lattice space
        translation = a * e1 * n + a * e2 * m
        interactions_dict[(n, m)] = load_gij(friction_coeffs_root, translation, order_g12, eps)

    return self_friction_dict, interactions_dict


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from carpet.physics.friction_pairwise import define_gmat_glob_and_q_glob

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

## Checks
# # Ben checked: mtwist solution -> classes of cilia with different inital phase have the same phase after N cycles
# # Is self-friction different if cilia are located differently? Check
# # Result: difference ~0.2%
# phis = np.linspace(0, 2 * np.pi, 50, endpoint=False)
# test_vals = []
# for key, g_ii in self_friction_dict.items():
#     gii_vals = [g_ii(phi,0) for phi in phis]
#     test_vals.append(gii_vals)

# test_vals = np.array(test_vals)
# mean_vals = np.mean(test_vals, axis=0)
# for i,vals in enumerate(test_vals):
#     delta = (vals - mean_vals ) / mean_vals
#     plt.plot(phis, delta)
# plt.show()

# # Test 2: g_ii(phi,phi); compare with mean_vals obtained earlier
# test_vals = []
# for key, g_ii in self_friction_dict.items():
#     gii_vals = [g_ii(phi,phi) for phi in phis]
#     test_vals.append(gii_vals)

# test_vals = np.array(test_vals)
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
