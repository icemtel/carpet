'''
- Generate triangular lattice
- In rectangular domain, assuming periodic boundary conditions.

Comment: It would be more natural to rewrite most of the functions here as methods of one class, but I don't do it to
(a) stay flexible, (b) avoid rewriting a lot of working code.
'''
import math
import numpy as np
from numpy.linalg import norm

from carpet.various import get_basis_dual, mod2pi


def get_cell_sizes(a):
    cell_length = a
    cell_height = a * 3 ** (1 / 2) / 2
    return cell_length, cell_height


def get_basis():
    e1 = np.array([1, 0])
    e2 = np.array([0.5, 3 ** (1 / 2) / 2])
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

    nmax = max(nx, ny)
    for n in range(- nmax, nmax):
        for m in range(0, nmax):
            x = n * a * e1 + m * a * e2  # position vector
            # position vector within bounds?
            if (x[0] >= 0 - eps) and (x[1] >= 0 - eps) and (x[0] < L1 - eps) and (x[1] < L2 - cell_height + eps):
                coords.append(x)
                lattice_ids.append((n, m))

    coords = np.array(coords)
    lattice_ids = np.array(lattice_ids)

    return coords, lattice_ids


def get_neighbours_list(coords, nx, ny, a, distances=(1,)):
    '''
    For each node looks for other nodes at specified distances (multiplied by lattice edge length `a`).
    Those nodes are saved in `N1` list. Relative positions are saved in `T1` list.
    WARNING: very slow for a large lattice; -> compute once and save to disk
           - 64x64: ~90 minutes
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


def get_neighbours_list_general(coords, nx, ny, a, connections):
    '''
    For each node looks for other nodes at specified translation vectors.
    Those nodes are saved in `N1` list. Relative positions are saved in `T1` list.
    :param connections: list of vectors - relative oscillator positions [units of length];
                        - if reciprocal, then for each vector `t`, the vector `-t` has to be included as well;
                        - an oscillator can be its own neighbour.
    :return: list of neighbours, list of relative neighbour positions
    '''
    if nx == 2 or ny == 2:
        import warnings
        warnings.warn("nx=2 or ny=2 => wrong number of neighbours (5 or 4)\n"
                      "some oscillators were supposed to be connected twice, but this is not implemented")

    eps = 10 ** -4 * a
    L1, L2 = get_domain_sizes(nx, ny, a)
    N = len(coords)
    # Check:
    distances = norm(connections, axis=1)
    if max(distances) >= max([L1, L2]):
        raise NotImplementedError("Assumption: d * a < max(L1,L2) is not satisfied")

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
                for t in connections:
                    if norm(t - translation) < eps:
                        return np.array(t)
        return None

    N1 = [[] for _ in coords]  # list of lists of neighbours indices
    T1 = [[] for _ in coords]  # list of lists of translation vectors between neighbours
    # loop over pairs of lattice points
    for i in range(N):
        for j in range(N):
            translation = get_neighbours(i, j)  # check if neighbours
            if translation is not None:  # is yes - add to the list
                N1[i].append(j)
                T1[i].append(translation)
    return N1, T1


### Wave vectors and reciprocal lattice ###

def get_basis_dual_domain(nx, ny, a):
    '''
    Reciprocal vectors for rectangular domain
    '''
    e1, e2 = get_basis()
    a1 = nx * a * e1
    a2 = ny * a * (2 * e2 - e1) / 2
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


def define_get_k(nx, ny, a):
    '''
    Checked: get_k is equivalent to get_k_naive: gives the same phase vectors mod 2pi
    '''
    import warnings
    warnings.warn("To be depricated! Returns vectors from a rectangular, rather than hexagonal cell (FBZ)")

    assert ny % 2 == 0  # check that ny is even

    def shift_integer(k, n, s):
        '''
        :param k,n,s: integers
        Assuming that k in Z/nZ (n-periodic integer),
        shift it to interval [-s,n-s]

        e.g. if s= n//2, the interval is:
        even: - n1 / 2      to  n1 / 2 - 1
        odd:  - (n1-1) / 2  to  (n1-1)/2
        '''
        return (k + s) % n - s

    a1dual, a2dual = get_basis_dual_domain(nx, ny, a)

    nxhalf = nx // 2
    nyhalf = ny // 2

    def get_k(k1, k2, shiftx=nxhalf, shifty=nyhalf):
        '''
        :param k1,k2: wavemodes
        :param shiftx,shifty: define additionally a shift
                            Get k in region, x: [- shiftx * a1dual / nx, a1dual - shiftx * a1dual / nx]
                                             y: ...
        :return:
        '''
        div = (k1 + shiftx) // nx
        k2 = k2 - ny * div // 2  # shift k2 as a response on shift in k1 - needed in triangular lattice

        k = shift_integer(k1, nx, shiftx) * a1dual + shift_integer(k2, ny, shifty) * a2dual

        return k

    return get_k


def define_shift_k_to_fbz(a):
    '''
    Defines a function, which get any wave vector, and shifts it into the first Brillouin zone
    :param a:
    :return:
    '''

    def project(vec, basis_vec):
        basis_vec = np.asarray(basis_vec)
        return vec @ basis_vec / (basis_vec @ basis_vec)

    def decompose_recip(k):
        # decompose a vector to vectors b1,b2,b3 (double of normals of the hexagon cell)
        ms = np.array([project(k, bi) for bi in bs])
        return ms

    b1, b2 = get_basis_dual_cell(a)
    b3 = b1 + b2
    bs = [b1, b2, b3]

    def k_to_fbz(k, eps=1e-8):
        k = np.array(k)
        num_iter = 0
        ms = decompose_recip(k)

        while np.amax(abs(ms)) > 0.5 + eps and num_iter < 10:  # +eps to acccount for numerical edge case
            i = int(np.argmax(abs(ms)))  # start with the direction with the biggest projection
            mi = ms[i]
            bi = bs[i]
            k -= bi * np.round(mi)  # shift by integer value
            ms = decompose_recip(k)
            num_iter += 1
        if num_iter == 10:
            raise ValueError("Didn't converge to a unit cell - check algorithm!")
        return k

    return k_to_fbz


def define_get_k_fbz(nx, ny, a):
    assert ny % 2 == 0  # check that ny is even
    get_k_naive = define_get_k_naive(nx, ny, a)
    k_to_fbz = define_shift_k_to_fbz(a)

    def get_k(k1, k2, eps=1e-8):
        k = get_k_naive(k1, k2)
        return k_to_fbz(k, eps)

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
    '''
    :return: a function -> see its description
    '''
    get_k = define_get_k_naive(nx, ny, a)
    mtwist_phi = np.zeros((nx, ny, nx * ny))

    for k1 in range(nx):
        for k2 in range(ny):
            # wave vector
            k = get_k(k1, k2)  # k1 * a1dual / nx + k2 * a2dual / ny
            for ix in range(nx * ny):
                mtwist_phi[k1, k2, ix] = mod2pi(- np.dot(k, coords[ix, :]))

    def get_mtwist(k1, k2):
        '''
        :return: phase vector, corresponding to a wave with vector k with wave numbers k1,k2
        '''
        if k1 < 0 or k1 >= nx:
            raise ValueError("k1 out of allowed interval")
        if k2 < 0 or k2 >= ny:
            raise ValueError("k2 out of allowed interval")
        return np.array(mtwist_phi[k1, k2])

    return get_mtwist


def define_get_mtwist_slow(coords, nx, ny, a):
    '''
    Returns a slower version of :get_mtwist: function
    Advantage: can give any positive or negative k1 and k2
    '''
    get_k = define_get_k_naive(nx, ny, a)
    coords_cp = np.array(coords)

    def get_mtwist(k1, k2):
        '''
        :return: phase vector, corresponding to a wave with vector k with wave numbers k1,k2
        '''
        k = get_k(k1, k2)
        mtwist = np.zeros(nx * ny)
        for ix in range(nx * ny):
            mtwist[ix] = mod2pi(- np.dot(k, coords_cp[ix, :]))
        return mtwist

    return get_mtwist


### Visualize reciprocal lattice: plot the First Brillouin Zone (FBZ) border

def plot_fbz(a, ax=None, **kwargs):
    """
    For triangular lattice FBZ is a hexagon
    :param a: lattice spacing
    :param kwargs: passed to Polygon
    :return:
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon

    rot_matrix = lambda angle: np.array(
        [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])  # rotation matrix

    hexagon_edge = 4 * np.pi / 3 / a  # checked via numerical experiment

    # Define hexagon
    def get_hexagon(scale=hexagon_edge, center=np.array([0, 0]), **kwargs):
        angles = np.linspace(0, 2 * np.pi, 6, endpoint=False)
        r = np.array((scale, 0))

        polygon_kwargs = dict(fc="none", ec="black", linestyle='--', alpha=0.5, lw=2)  # default params
        polygon_kwargs.update(kwargs)  # update with input kwargs
        return Polygon([center + rot_matrix(angle) @ r for angle in angles], **polygon_kwargs)

    # Plot hexagon
    if ax is None:
        ax = plt.gca()
    ax.add_artist(get_hexagon(hexagon_edge, **kwargs))
    return ax


### Friction coefficients - ver 1. ###
## Use functions in "physics" instead.
def get_connections():
    '''
    :return: Relative positions of neighbouring cilia in lattice coordinates
             First order neighbours
    '''
    return [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, 1), (1, -1)]


# Cilia coupling. Must be removed in the future - call module with coupling


def define_gmat_glob_and_q_glob(set_name, a, neighbours_indices, neighbours_rel_positions,
                                order_g11, order_g12, T, use_numba=True):
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

    warnings.warn("To be depricated! Import 'physics' instead")

    e1, e2 = get_basis()
    return define_gmat_glob_and_q_glob0(set_name, e1, e2, a,
                                        neighbours_indices, neighbours_rel_positions,
                                        order_g11, order_g12, T, use_numba=use_numba)


def define_right_side_of_ODE(gmat_glob, q_glob):
    import carpet.physics.friction_pairwise as coupling
    import warnings

    warnings.warn("To be depricated! Import 'physics' instead")

    return coupling.define_right_side_of_ODE(gmat_glob, q_glob)


if __name__ == '__main__':
    # OK: both triangular and rectangular lattice
    # OK: translations lengths (plot)
    # OK: translations direction (plot)
    # OK: mtwist no duplicates
    # OK: mtwist (plot)
    # OK: new get_neighbours

    import matplotlib.pyplot as plt
    import carpet.visualize as visualize

    a = 18
    nx = 8
    ny = 8  # must be even

    coords, lattice_ids = get_nodes_and_ids(nx, ny, a)
    # N1, T1 = get_neighbours_list(coords, nx, ny, a, distances=(1, np.sqrt(3), 2))
    # print("Neighbours as array shape:", np.array(N1).shape)

    ## New get_neighbours:
    translations = [a * np.array([np.cos(psi), np.sin(psi)])  # First neighbours
                    for psi in np.linspace(0, 2 * np.pi, 6, endpoint=False)] \
                   + [np.sqrt(3) * a * np.array([np.cos(psi), np.sin(psi)])
                      for psi in [np.pi / 2, 3 * np.pi / 2]]  # 2nd neighbour (only 1)
    N1, T1 = get_neighbours_list_general(coords, nx, ny, a, translations)
    print("Neighbours as array shape:", np.array(N1).shape)
    ## Visualize
    visualize.plot_edges(coords, T1)
    visualize.plot_nodes(coords)
    visualize.plt.show()

    ## Non-reciprocal neighbours: (interactions in one direction)
    translations = [a * np.array([np.cos(psi), np.sin(psi)])  # First neighbours
                    for psi in np.linspace(0, np.pi, 3, endpoint=False)] \
                   + [np.sqrt(3) * a * np.array([np.cos(psi), np.sin(psi)])
                      for psi in [np.pi / 2]]  # 2nd neighbour (only 1)
    N1, T1 = get_neighbours_list_general(coords, nx, ny, a, translations)
    print("Neighbours as array shape:", np.array(N1).shape)
    ## Visualize
    visualize.plot_edges(coords, T1)
    visualize.plot_nodes(coords)
    visualize.plt.show()

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
    ## Result: they are equivalent
    get_k = define_get_k_fbz(nx, ny, a)
    get_k_naive = define_get_k_naive(nx, ny, a)

    for k1 in range(nx):
        for k2 in range(ny):
            print(k1, k2)

            k = get_k(k1, k2)
            k_naive = get_k_naive(k1, k2)

            for c in coords:
                assert np.allclose(np.exp(1j * k_naive @ c), np.exp(1j * k @ c))

    # Check: get_k_fbz maps to FBZ (hexagon)
    b1, b2 = get_basis_dual_cell(a)
    n1 = b1 / 2  # hexagon_from_edge_to_origin *np.array([0,1])
    n2 = b2 / 2
    n3 = n1 + n2
    ns = [n1, n2, n3]

    for _ in range(100):
        k1, k2 = np.random.randint(100, size=2)  # test on random wave numbers
        k_naive = get_k_naive(k1, k2)
        k = get_k(k1, k2)

        # Check that wave vectors are equivalent
        for c in coords:
            assert np.allclose(np.exp(1j * k_naive @ c), np.exp(1j * k @ c))

        # Check that it's whithin the hexagon borders
        for nvec in ns:
            assert abs(k @ nvec) / (nvec @ nvec) <= 1 + 1e-8  # small number to account for numerical edge cases

    ## Dual basis?
    print(get_cell_sizes(a))
    print(get_basis_dual_domain(nx, ny, a))


    plot_fbz(a)
    plt.xlim([-0.25, 0.25])
    plt.ylim([-0.25, 0.25])
    plt.show()
