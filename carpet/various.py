'''
- root mean square
- logging
- phases_to_interval: map phases to  mean+-2pi interval without changing mean
- get_phase_random
'''

import logging
import sys
import numpy as np


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


def setup_logging(filename=None, mode='a', print_log_messages=True, level=logging.INFO):
    """
    Setup logging: use `logging.info(msg)` or `logger.info(msg)` to log something into a file and console.
    """
    logger = logging.getLogger("")  # Use root logger   #  logging.getLogger(name) # or use a named one?
    # Clear previous handlers, if any
    if (logger.hasHandlers()):
        logger.handlers.clear()
    # Formatter
    _formatter = logging.Formatter('{asctime} {levelname} {threadName} {name} {message}', style='{')
    # Handler
    if filename is not None:
        log_handler = logging.FileHandler(filename, mode=mode)
        logger.addHandler(log_handler)
        log_handler.setFormatter(_formatter)
        log_handler.setLevel(level)
    if print_log_messages:
        stream_handler = logging.StreamHandler(sys.stdout)
        stream_handler.setFormatter(_formatter)
        logger.addHandler(stream_handler)
        stream_handler.setLevel(level)
    # Set level - pass all messages with level INFO or higher (WARNING, ERROR, CRITICAL)
    logger.setLevel(level)

    return logger


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


def get_phase_random(N):
    '''
    Uniformly distributed random phase vector - on the interval from 0 to 2pi.
    '''
    phi = 2 * np.pi * np.random.rand(N)
    return phi


def define_dump_object(default_output_folder='.'):
    '''
    Define a function-shortcut for pickling objects.
    - Pickling is a way to save python objects to the disk space.
    - complex objects might get broken (e.g. matplotlib figure)
    - works well with dictionaries, numpy.arrays
    - There is a chance that won't be able to load objects on a different python/libraries version.
    '''
    import pickle, os

    def dump_object(obj, filename, path=default_output_folder):
        filename = os.path.join(path, filename)
        print(filename)
        with open(filename, 'wb') as f:
            pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

    return dump_object


def define_load_object(default_output_folder='.'):
    '''
    Define a function-shortcut for loading (un-pickling) objects.
    - There is a chance that won't be able to load objects on a different python/libraries version.
    '''
    import pickle, os

    def load_object(filename, path=default_output_folder):
        filename = os.path.join(path, filename)
        with open(filename, 'rb') as f:
            obj = pickle.load(f)
        return obj

    return load_object




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


def get_basis_dual(a1,a2):
    '''
    :param a1, a2: basis vectors
    :return: dual basis vectors; https://en.wikipedia.org/wiki/Reciprocal_lattice
    '''
    R = np.array([[0, 1], [-1, 0]])  # rotation by 90deg
    b1 = 2 * np.pi * (R @ a2) / (a1 @ (R @ a2))
    b2 = 2 * np.pi * (R @ a1) / (a2 @ (R @ a1))
    return b1, b2