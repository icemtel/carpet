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



def my_argwhere(boolean_array):
    """
    Feed in some condition, e.g. my_argwhere(A == 1), return indices, where it is True.
    Extra line of code makes the dimensions the way I expect them:
    e.g. [[1],[2]] -> [1,2]
    Bulky because of edge cases: [[1]] -> [1]; [1] -> [1]
    Tested for 1D arrays.
    """
    ind_where = sp.argwhere(boolean_array)
    if ind_where.ndim > 1:
        ind_where = sp.array([sp.squeeze(ind) for ind in ind_where], dtype=sp.int64)
    return ind_where


def get_sorted_ind_periodic(phi, eps=10**-8):
    """
    (1) Sort phases in ascending order
    (2) If there are phases, equal to 2pi up to `small_number`, move them to the beginning of the array
    (3) Returns new index: phi_new = phi[sorted_ind]
    :param eps: a small number
    """
    # (1)
    presorted_ind = sp.argsort(phi)
    # (2)
    close_to_2pi = abs(phi[presorted_ind] - 2 * sp.pi) < eps  # Array of True/False
    ind1 = my_argwhere(close_to_2pi)  # indices of those elements which are close to 2pi
    ind2 = my_argwhere(sp.logical_not(close_to_2pi))
    rearranged_ind = sp.concatenate((ind1, ind2))  # rearranged indices after the sorting
    # (3)
    sorted_ind = presorted_ind[rearranged_ind]  # rearranged indices before the sorting; Checked: correct order
    return sorted_ind


def get_classes(phi_k, eps=10**-8):
    """
    Which cilia have the same phase values?
    define classes of nodes, according to which nodes have the same phase for the given m-twist solution
    :param eps: a small number
    """
    phi_k = sp.array(phi_k)
    sorted_ind = get_sorted_ind_periodic(phi_k)
    phi_k_sorted = phi_k[sorted_ind]
    # Find jumps:
    # Check if phi_j is equal to phi_(j-1) up to modulus 2pi
    diffs = sp.diff(phi_k_sorted)
    jumps = abs(sp.sin(diffs / 2)) > eps # Array of True/False

    # ... compute list of class IDs for permuted set of nodes; usage: class_ids[j] = "class of node ind_sorted_[j]"
    class_ids_sorted = sp.cumsum(jumps)
    class_ids_sorted = sp.concatenate(([0], class_ids_sorted))

    # compute look-up table of class IDs for nodes; usage: ix_to_class[ix] = "class of node ix"
    # ... construct inverse permutation for ind
    N = len(phi_k) # number of cilia
    ind_inv = sp.zeros(N, dtype=int)
    for j in range(N):
        ind_inv[sorted_ind[j]] = j
    # ... look-up table of class for each node: reshuffle class_ids
    ix_to_class = class_ids_sorted[ind_inv]

    # compute look-up table for list of nodes belonging to each class
    # number of distinct classes; equal to the last element + 1 since we sorted it; and start from 0
    nclass = class_ids_sorted[-1] + 1
    class_to_ix = []
    for iclass in range(nclass):
        # class_to_ix.append( np.argwhere(np.array(ix_to_class) == iclass) )
        iclass_to_ix = my_argwhere(ix_to_class == iclass)
        class_to_ix.append(iclass_to_ix)
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

    phi_k = sp.array(phi_k)

    for c in class_to_ix:
        # Get neighbours parameters of the first cilium in a class - compare others with it
        first_neighbours_ix = N1[c[0]]
        first_neighbours_translations = T1[c[0]]  # position of the neighbour, relative to the cilium
        first_neighbours_phases = phi_k[first_neighbours_ix]
        first_sorted_ind = get_sorted_ind_periodic(first_neighbours_phases)
        first_phases_sorted = sp.array([mod(x) for x in first_neighbours_phases[first_sorted_ind]])
        first_translations_sorted = first_neighbours_translations[first_sorted_ind]
        N_neighbours = len(first_neighbours_ix)

        for ix in c[1:]:
            neighbours_ix = N1[ix]
            neighbours_translations = T1[ix]
            neighbours_phases = phi_k[neighbours_ix]
            sorted_ind = get_sorted_ind_periodic(neighbours_phases)
            phases_sorted = sp.array([mod(x) for x in neighbours_phases[sorted_ind]])
            translations_sorted = neighbours_translations[sorted_ind]

            if norm(first_phases_sorted - phases_sorted)  * N_neighbours ** (-1/2) > eps:

                raise ValueError("Classes check is not passed: phase difference is too big: {:.3e}"
                                 .format(norm(first_phases_sorted - phases_sorted)  * N_neighbours ** (-1/2)))
            if norm(first_translations_sorted - translations_sorted) * N_neighbours ** (-1/2) > eps:
                print(ix)
                print(first_phases_sorted)
                print(phases_sorted)
                print(first_translations_sorted)
                print(translations_sorted)
                raise ValueError("Classes check is not passed: neighbours have different relative positions, diff: {:.3e}"
                                 .format(norm(first_translations_sorted - translations_sorted) * N_neighbours ** (-1/2)))

    print("Test passed successfully!")

if __name__ == '__main__':
    # phi_k = (0, 2.5, 0, 2.5)
    # ix_to_class, class_to_ix = get_classes(phi_k)
    # print(ix_to_class) # [0, 1, 0, 1]
    #
    # phi_k = (0, 2.5, 0, 2.5 + 2 * sp.pi)
    # ix_to_class, class_to_ix = get_classes(phi_k)
    # print(ix_to_class) # [0, 1, 0, 1]

    phi_k = (0, 2.5, 3.5, 2 * sp.pi, 2.5, 3.5)
    ix_to_class, class_to_ix = get_classes(phi_k)
    N1 = sp.array([[5,1], [0,2],[1,3],[2, 4],[3,5],[4,0]])
    T1 = sp.array([[-1,1],[-1,1],[-1,1],[-1,1],[-1,1],[-1,1]])

    print(check_class_neighbours(phi_k, class_to_ix, N1, T1))