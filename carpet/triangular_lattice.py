'''
Generate triangular lattice
'''
import scipy as sp
from scipy.linalg import norm
import carpet.friction as friction

def get_cell_sizes(a):
    cell_length = a
    cell_height = a * 3 ** (1 / 2) / 2
    return cell_length, cell_height

def get_basis():
    e1 = sp.array([1, 0])
    e2 = sp.array([0.5, sp.sqrt(3) / 2])
    return e1,e2

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

    e1 ,e2 = get_basis()

    cell_length, cell_height = get_cell_sizes(a)
    L1, L2 = get_domain_sizes(nx, ny, a)
    eps = 10 ** -4 * a  # better safe than sorry: small number to account for rounding errors below

    # generate triangular lattice
    coords = []  # list of position vectors for honeycomb lattice
    lattice_ids = []  # list of corresponding lattice indices (n,m)

    for n in range(-4 * nx, 4 * nx):
        for m in range(0, 4 * ny):
            x = n * a * e1 + m * a * e2  # position vector
            # position vector within bounds?
            if (x[0] >= 0 - eps) and (x[1] >= 0 - eps) and (x[0] < L1 - eps) and (x[1] < L2 - cell_height + eps):
                coords.append(x)
                lattice_ids.append((n, m))

    # discard entries not used
    coords = sp.array(coords)
    lattice_ids = sp.array(lattice_ids)

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

### Friction coefficients - ver 1. ###


def get_connections():
    '''
    :return: Relative positions of neighbouring cilia in lattice coordinates
    '''
    return [(-1,0),(1,0),(0,-1),(0,1),(-1,1),(1,-1)]

def define_gmat_glob_and_q_glob_triangular(set_name, a, neighbours_indices, neighbours_rel_positions,
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
    connections = get_connections()
    e1, e2 = get_basis()
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
#  - Cubic interpolation for gamma_ij


if __name__ == '__main__':
    # OK: both triangular and rectangular lattice
    # Plot:
    # OK: translations lengths
    # OK: translations direction

    import visualize

    a = 18
    nx = 6
    ny = 6

    coords, lattice_ids = get_nodes_and_ids(nx, ny, a)
    N1, T1 = get_neighbours_list(coords, nx, ny, a)

    ## Visualize
    # visualize.plot_edges(coords, T1)
    # visualize.plot_nodes(coords)
    # visualize.plt.show()


    ## Friction
    order_g11 = (8,0)
    order_g12 = (4,4)
    T = 31.25
    gmat, qglob = define_gmat_glob_and_q_glob_triangular('machemer_1', a, N1, T1,order_g11, order_g12,T)

