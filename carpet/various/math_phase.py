import numpy as np
import math
from .math_lattice import mod2pi


def get_phase_random(N):
    '''
    Uniformly distributed random phase vector - on the interval from 0 to 2pi.
    '''
    phi = 2 * np.pi * np.random.rand(N)
    return phi


def rms(vec):
    '''
    Root mean square; Quadratic averarge
    '''
    return np.sqrt(abs(np.mean(vec * np.conj(vec))))


def cexp_dist(phi1, phi2=0):
    '''
    Distance (normalized by sqrt(ndim), i.e. rms)
    between complex exponents (as vectors in 2N-Euclidean space);
    Equivalent to rms(np.exp(1j * phi1) - np.exp(1j * phi2))
    '''
    return rms(2 * np.sin(0.5 * (phi1 - phi2)))


def cexp(vec):
    '''
    Calculate complex exponent
    '''
    return np.exp(1j * vec)



def phases_to_interval(phi):
    '''
    Take phase vector.
    Return equivalent phase vector with
    phases in the +-2pi interval, centered at the mean of the initial phase vector (mean unchanged)

    Tested  - doesn't break.
    Problem: one vector has more than 1 representation in that interval
    '''
    x = np.array(phi)
    xmin = x.mean() - 2 * np.pi
    xmax = x.mean() + 2 * np.pi
    flag = (x.max() > xmax) or (x.min() < xmin)
    while flag:
        imax = x.argmax()
        imin = x.argmin()
        x[imax] -= 2 * np.pi
        x[imin] += 2 * np.pi

        flag = (x.max() > xmax) or (x.min() < xmin)

    return x



## Global phase
### Circular Average
from scipy.stats import circmean, circstd, circvar


def circ_dist(phi, phi0=0, axis=None):
    '''
    If two phase vectors age given, calculate for dphi=phi-phi0
    Without axis option is equivalent to math.sqrt(- 2 * math.log(abs(np.exp(1j *(phi0 - phi)).mean())))
    '''
    return circstd(phi - phi0, axis=axis)


def complex_exp_mean(phi, phi0=0):
    '''
    Mean of complex exponents of array of phases phi, relative to phi0
    '''
    return np.mean(np.exp(1j * (phi - phi0)))


def define_circular_mean_phase(phi0):
    '''
    :return: function which computes circular mean phase of phi, relative to phi0.
    '''

    def get_phi_global(phi):  # equivalent up to mod 2pi to np.angle(np.mean(np.exp(1j * (phi - phi0))))
        return circmean(phi - phi0)

    return get_phi_global


def order_parameter(phi, phi0=0):
    return abs(complex_exp_mean(phi, phi0))


## Centralized phase difference
def phase_diffs_centralized(phis, phi0):
    '''
    :param phis: list of phase vectors
    :param phi0: phase vector to subtract
    :return: list of phase-vector differences; centralized,
            i.e. circular mean of the phase vector differences equals to zero
            [For movies of phase deviations]
    '''
    dphis = (phis - phi0) % (2  *np.pi)
    for i, dphi in enumerate(dphis): # Centralzie
        dphis[i] -= circmean(dphi)
    return dphis

### Mean phase
def get_mean_phase(phi):  # keep for compatibility
    return np.mean(phi)


def mean_phase(phi):
    return np.mean(phi)


if __name__ == '__main__':
    ## Test RMS
    a = np.array([1, 1])
    print(rms(a))
    print(phases_to_interval(a))
    print(get_phase_random(3))
    ## Test Logger:
    #   logger = setup_logging('a.log')
    # #  logger = setup_logging()
    #
    #   logging.info("THIS IS INFO!!!")
    #   logging.warning("THIS IS WARNING!")
    #   logging.debug("THIS IS A DEBUG MESSAGE! IT SHOULDN't BE DISPLAYED!!")
    #   logger.info("Can also log this")
