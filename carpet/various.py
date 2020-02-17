'''
- root mean square
- logging
'''

import logging
import sys
import numpy as np


def rms(vec):
    '''
    Root mean square; Quadratic averarge
    '''
    return np.sqrt(abs(np.mean(vec * np.conj(vec))))


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




if __name__ == '__main__':
    ## Test RMS
    a = np.array([1, 1])
    print(rms(a))
    print(phases_to_interval(a))

    ## Test Logger:
    #   logger = setup_logging('a.log')
    # #  logger = setup_logging()
    #
    #   logging.info("THIS IS INFO!!!")
    #   logging.warning("THIS IS WARNING!")
    #   logging.debug("THIS IS A DEBUG MESSAGE! IT SHOULDN't BE DISPLAYED!!")
    #   logger.info("Can also log this")
