"""
Kuramoto-type model.

Sinusoidal coupling or any other

NB:
- sign of coupling
- dimensionless coupling strength
"""
import numpy as np
import math


def define_sine_coupling(sin_str):
    '''
    NB: in papers usually the opposite sign. For consistency, I keep the sign like in my previous simulations
    '''
    def coupling(phi1, phi2):
        c = sin_str * math.sin(phi1 - phi2)
        return c

    return coupling

def define_sine_cosine_coupling(sin_str, cos_str):
    def coupling(phi1, phi2):
        c = sin_str * math.sin(phi1 - phi2) + cos_str * math.sin(phi1 - phi2)
        return c

    return coupling

def define_right_side_of_ODE(coupling_func, freq,  N1, T1):
    '''
    :param coupling_func: pairwise interactions: f(phi1,phi2). Dimensionless! Will be multiplied by frequency
    :param freq: frequency of a single oscillator
    :param N1: list of lists of neighbours indices
    :param T1: UNUSED - keep for generalizations
    :return:
    '''
    def right_side_of_ODE(t, phi):
        # Self-driving (multiply by frequency later
        right_side = np.ones([len(N1)])
        # Interactions
        # MAYBE: Optimize -  sum array for; or numba
        for ix, (neighbours, translations) in enumerate(zip(N1, T1)):
            for ix2 in neighbours:
                right_side[ix] += coupling_func(phi[ix], phi[ix2])  # neighbour_sign(ix,ix2) *

        right_side *= freq  # proper dimensions
        return right_side

    return right_side_of_ODE
