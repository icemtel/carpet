'''
Take advantage of carpet symmetries by implementing "classes" procedure.
If a number of cilia has the same phase - map them into one class.
Then we can solve ODE taking only 1 cilium from each class, saving some resources.

Important:
- If two cilia have the same phase, it doesn't in general imply that they will have the same phase in every moment of time.
  This will be True, only if those cilia have identical neighbours, and those neighbours have identical neighbours, etc..
  For m-twist solutions this holds True, and as a quick check below there is a function check_class_neighbours(), which
  checks phases of the first level neighbours.
'''
import scipy as sp
from scipy.linalg import norm
import copy


def get_classes(phi, eps=1e-8):
    """
    Which cilia have the same phase values?
    define classes of nodes, according to which nodes have the same phase for the given m-twist solution
    :param eps: a small number
    """
    phi = sp.array(phi)

    unclassified_list = [(i, phii) for i, phii in enumerate(phi)]

    ix_to_class = sp.full_like(phi, fill_value=sp.nan, dtype=int)
    class_to_ix = []

    class_id = 0
    while unclassified_list:  # check if not empty
        i, phii = unclassified_list[0]  # take the next oscillator without a class

        class_list = []
        classified_ixs = []  # indices to remove from unclassified_list
        for ix, (j, phij) in enumerate(unclassified_list):
            diff = abs(sp.exp(1j * (phii - phij)) - 1)
            if diff < eps:  # Add to the class if difference is small; take into account periodicity
                ix_to_class[j] = class_id
                class_list.append(j)
                classified_ixs.append(ix)

        class_to_ix.append(sp.array(class_list, dtype=int))
        unclassified_list = [i for j, i in enumerate(unclassified_list) if j not in classified_ixs]
        class_id += 1
    return ix_to_class, class_to_ix


def get_neighbours_list_class(unique_cilia_ids, ix_to_class, N1, T1):
    # Define a class number as the position in `unique_cilia_ids` array
    # Rewrite neighbours of cilia in terms of classes
    N1_class = []
    T1_class = []
    for id in unique_cilia_ids:
        translations = T1[id]
        T1_class.append(copy.copy(translations))

        # identify to which class the neighbour belongs to
        neighbours = N1[id]
        neighbour_classes = []
        for neighbour_ix in neighbours:
            neighbour_class = ix_to_class[neighbour_ix]
            neighbour_classes.append(neighbour_class)
        N1_class.append(neighbour_classes)
    return N1_class, T1_class

def get_unique_cilia_ix(class_to_ix):
    """
    Get one cilium from each class. Return their indices.
    """
    nclass = len(class_to_ix)
    unique_cilia_ix = sp.array([class_to_ix[iclass][0] for iclass in range(nclass)], dtype=sp.int64)
    return unique_cilia_ix

def check_class_neighbours(phi_k, class_to_ix, N1, T1, eps=10**-8):
    """
    This one is not very clean, but it checks if each cilium in a class have identical neighbours:
    - Take one cilium from a class
    - For every other cilium in this class check that it has neighbours with the same phases and with the same relative positions
    as neighbours of the first cilium
    - Raises a error message if that's not True
    - Otherwise prints a success message
    """
    def mod(x):
        '''
        fmod(x,y) is not equivalent to (x % y): https://docs.python.org/3/library/math.html and
        is preferred when working with floats
        :return: a value in interval from 0 to 2pi
        '''
        import math

        x = math.fmod(x, 2 * sp.pi)
        if x < 0:
            x += 2 * sp.pi
        return x

    phi_k = [mod(phi) for phi in phi_k]
    phi_k = sp.array(phi_k)

    for c in class_to_ix:
        # Get neighbours parameters of the first cilium in a class - compare others with it
        first_neighbours_ix = N1[c[0]]
        first_neighbours_translations = sp.array(T1[c[0]])  # position of the neighbour, relative to the cilium
        first_neighbours_phases = phi_k[first_neighbours_ix]

        # Go through cilia from the same class
        for ix in c[1:]:
            neighbours_ix = N1[ix]
            neighbours_translations = T1[ix]
            neighbours_phases = phi_k[neighbours_ix]

            checked_flag = sp.full_like(neighbours_ix, dtype=bool, fill_value=False)
            for i, (phase, translation) in enumerate(zip(neighbours_phases, neighbours_translations)):
                for phase1, translation1 in zip(first_neighbours_phases, first_neighbours_translations):
                    if abs(sp.exp(1j *(phase-phase1))-1) < eps:
                        if norm(translation - translation1) < eps:
                            if checked_flag[i] == False:
                                checked_flag[i] = True
                            else:
                                raise ValueError("Found a second matching neighbour")

            if not sp.all(checked_flag):
                print(ix)
                print(first_neighbours_phases)
                print(neighbours_phases)
                print(first_neighbours_translations)
                print(neighbours_translations)
                raise ValueError("Not all neighbours are matched")

    print("Test passed successfully!")

if __name__ == '__main__':
    # phi_k = (0, 2.5, 0, 2.5)
    # ix_to_class, class_to_ix = get_classes(phi_k)
    # print(ix_to_class) # [0, 1, 0, 1]
    #
    # phi_k = (0, 2.5, 0, 2.5 + 2 * sp.pi)
    # ix_to_class, class_to_ix = get_classes(phi_k)
    # print(ix_to_class) # [0, 1, 0, 1]

    #phi_k = (0, 2.5, 3.5, 2 * sp.pi, 2.5, 3.5)
    # N1 = sp.array([[5,1], [0,2],[1,3],[2, 4],[3,5],[4,0]])
    # T1 = sp.array([[-1,1],[-1,1],[-1,1],[-1,1],[-1,1],[-1,1]])
    #
    """
    Tested on all m-twists of lattice_triangular and lattice_triangular2 6x6
    """
    import carpet.lattice_triangular2 as lattice
    a = 18
    nx = 6
    ny = 6 # must be even
    N = nx * ny

    # for k1 in range(nx):
    #     for k2 in range(ny):
    #         print(k1, k2)

    k1,k2 = 1,2

    coords, lattice_ids = lattice.get_nodes_and_ids(nx, ny, a)
    N1, T1 = lattice.get_neighbours_list(coords, nx, ny, a)

    get_mtwist = lattice.define_get_mtwist(coords, nx, ny, a)
    phi_k = get_mtwist(k1, k2)
    phi_k -= phi_k.mean()

    ix_to_class, class_to_ix = get_classes(phi_k)

    # Test if each cilium is classified only once
    assert len(set(sp.concatenate(class_to_ix))) == N

    # Test if classes are balanced
    ncl = [len(cl) for cl in class_to_ix]
    print("Classes lengths:", ncl)

    num_classes = len(ncl)
    # Check symmetry of class
    print(num_classes, set(get_mtwist(k1,k2)))
    # Run a second test:
    check_class_neighbours(phi_k, class_to_ix, N1, T1)