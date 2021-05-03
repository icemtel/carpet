import math
import numpy as np

def get_basis_dual(a1, a2):
    '''
    For the lattice
    :param a1, a2: basis vectors
    :return: dual basis vectors; https://en.wikipedia.org/wiki/Reciprocal_lattice
    '''
    R = np.array([[0, 1], [-1, 0]])  # rotation by 90deg
    b1 = 2 * np.pi * (R @ a2) / (a1 @ (R @ a2))
    b2 = 2 * np.pi * (R @ a1) / (a2 @ (R @ a1))
    return b1, b2

def mod2pi(x):
    '''
    fmod(x,y) is not equivalent to (x % y): https://docs.python.org/3/library/math.html and
    is preferred when working with floats
    :return: a value in interval from 0 to 2pi
    '''
    x = math.fmod(x, 2 * np.pi)
    if x < 0:
        x += 2 * np.pi
    return x