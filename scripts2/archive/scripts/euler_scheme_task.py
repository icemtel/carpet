'''
2019-09-23 Anton Solovev

TODO: Euler scheme to integrate
Note: below y is a vector. use numpy arrays to represent vectors
Note: numpy and scipy arrays are the same thing, but scipy library provides more functionality


TODO: Test on example
Test the integration on a simple ODE example
You can make plots if you find it useful

TODO: Do you already know what to change in the `integrate_euler` function if you want to add noise?

You can send me back either a python file, or a jupyter notebook. Both are fine.
I use Python 3.6, but other versions should work as well.
'''

import scipy as sp
import math
import matplotlib.pyplot as plt


def integrate_euler(y0, fun, dt, t_span, eps=10 ** -8):
    """
    y' = f(t,y)
    y(t0) = y0
    :param y0: Initial state, array
    :param fun: Right-hand side of the system f(t,y)
    :param dt: Time step
    :param t_span: tuple (t0, tf) - start and end of integration interval
    :param eps: If t_span can't be divided in integer number of steps of size `dt`.
                The last time step will end at time `t_span[1]`, and the last step will have a different length,
                bigger than or equal to `eps`, but smaller than `dt + eps`.
    :return: (ys, ts)
             Where ys: list of states y(t_i), ts: list of times t_i
    """
    pass
