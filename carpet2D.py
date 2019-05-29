'''
Generate triangular lattice
'''
import scipy as sp
from scipy.linalg import norm


def get_cell_dimensions(e1, e2, a):
    cell_length = a * max(e1[0], e2[0])
    cell_height = a * max(e1[1], e2[1])
    return cell_length, cell_height


def get_domain_dimensions(e1, e2, l1, l2, a):
    cell_length, cell_height = get_cell_dimensions(e1, e2, a)
    L1 = cell_length * (l1 + 1)
    L2 = cell_height * (l2 + 1)
    return L1, L2


def get_nodes_and_ids(e1, e2, l1, l2, a):
    '''
    Create a triangular lattice in rectangular domain
    '''
    # unit vectors
    # spatial dimension of simulation domain
    cell_length, cell_height = get_cell_dimensions(e1, e2, a)
    L1, L2 = get_domain_dimensions(e1, e2, l1, l2, a)
    eps = 0.001 * a  # better safe than sorry: small number to account for rounding errors below

    # generate triangular lattice
    coords = []  # list of position vectors for honeycomb lattice
    lattice_ids = []  # list of corresponding lattice indices (n,m)

    for n in range(-4 * l1, 4 * l1):
        for m in range(0, 4 * l2):
            x = n * a * e1 + m * a * e2  # position vector
            # position vector within bounds?
            if (x[0] >= 0 - eps) and (x[1] >= 0 - eps) and (x[0] < L1 - eps) and (x[1] < L2 - cell_height + eps):
                coords.append(x)
                lattice_ids.append((n, m))

    # discard entries not used
    coords = sp.array(coords)
    lattice_ids = sp.array(lattice_ids)

    return coords, lattice_ids


def get_neighbours_list(coords, e1, e2, l1, l2, a):
    eps = 10 ** -3 * a
    L1_full, L2_full = get_domain_dimensions(e1, e2, l1, l2, a)

    N = len(coords)

    ## find nearest neigbors
    def get_neighbours(i, j):
        '''
        If the distance between two points is equal to the lattice spacing, return vector connecting them, else None.
        Takes into account lattice periodicity
        OK: sign
        TODO: obtain the same information from lattice indices instead
        '''
        for a1 in range(-1, 2):
            for a2 in range(-1, 2):
                translation = coords[j, :] - coords[i, :] + [a1 * L1_full, 0] + [0, a2 * L2_full]
                if a - eps < norm(translation) < a + eps:
                    return translation
        return None

    ## TO change in case of rectangular domain
    N1 = [[] for _ in coords]  # list of lists of neighbours
    T1 = [[] for _ in coords]  # list of translation vectors between neighbours
    # loop over pairs of lattice points
    for i in range(N):  # TODO: optimize
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


# TODO: test functions above


def get_gmat_glob():
    pass


def get_q_glob():
    pass


if __name__ == '__main__':
    # OK: both triangular and rectangular lattice
    # Plot:
    # OK: translations
    # OK: translation direction

    import matplotlib.pyplot as plt
    e1 = sp.array([1, 0])
#    e2 = sp.array([0,1])
    e2 = sp.array([1 / 2, 3 ** (1 / 2) / 2])
    a = 10
    l1 = 5
    l2 = 5

    coords, lattice_ids = get_nodes_and_ids(e1, e2, l1, l2, a)
    N1, T1 = get_neighbours_list(coords, e1, e2, l1, l2, a)

    # TEST: ... plot edges based on translations list
    eps = 0.001 * a
    for i, (neighbours, translations) in enumerate(zip(N1, T1)):
        for j, translation in zip(neighbours, translations):
            points_to_plot = sp.array([coords[i], coords[i] + translation])
            angle = sp.arctan2(translation[1], translation[0])
            if  abs(angle - sp.pi / 3) < eps:
                code = 'g-.'
            elif abs(angle + sp.pi * 2 / 3) < eps:
                code = 'r:'
            else:
                code = 'b'
            plt.plot(points_to_plot[:,0],points_to_plot[:,1], code)

    # ... plot nodes
    plt.plot(coords[:, 0], coords[:, 1], 'r.', markersize=12)
    plt.title('Triangular lattice')
    plt.gca().set_aspect('equal')
    plt.show()