import scipy as sp


def gaussian(N):  # returns random values, distributed as gaussian
    return sp.randn(N)  # with the same dimension as number of oscillators


def integrate_euler(y0, fun, D, dt, t_span, eps=10 ** -8):
    """
    y' = f(t,y)
    y(t0) = y0
    :param y0: Initial state, array
    :param fun: Right-hand side of the system f(t,y)
    :param D: Diffusion coefficient; <xi(t),xi(t')> = 2 D delta(t-t')
    :param dt: Time step
    :param t_span: tuple (t0, tf) - start and end of integration interval
    :param eps: If t_span can't be divided in integer number of steps of size `dt`.
                The last time step will end at time `t_span[1]`, and the last step will have a different length,
                bigger than or equal to `eps`, but smaller than `dt + eps`.
    :return: (ys, ts)
             Where ys: array of states y(t_i), ts: array of times t_i
    """
    t = t_span[0]
    t_end = t_span[1]
    y = sp.array(y0)
    N = len(y0)  # number of oscillators
    noise_coeff = (2 * D * dt) ** (1 / 2)

    ys = [y]
    ts = [t]
    while t < t_end - dt - eps:
        dy = fun(t, y) * dt + noise_coeff * gaussian(N)
        y = y + dy  # don't  use += on vectors!
        t += dt

        ys.append(y)
        ts.append(t)

    # The last step
    dt = t_end - t
    noise_coeff = (2 * D * dt) ** (1 / 2)
    dy = fun(t, y) * noise_coeff * gaussian(N)
    y = y + dy
    t += dt

    ys.append(y)
    ts.append(t)

    return sp.array(ys), sp.array(ts)


# # TODO: Correct? Check
# # TODO: try, maybe this version is faster?
# def integrate_euler_TEST(y0, fun, D, dt, t_span, eps=10 ** -8):
#     """
#     y' = f(t,y)
#     y(t0) = y0
#     :param y0: Initial state, array
#     :param fun: Right-hand side of the system f(t,y)
#     :param D: Diffusion coefficient; <xi(t),xi(t')> = 2 D delta(t-t')
#     :param dt: Time step
#     :param t_span: tuple (t0, tf) - start and end of integration interval
#     :param eps: If t_span can't be divided in integer number of steps of
#         size `dt`.
#                 The last time step will end at time `t_span[1]`, and the last
#                 step will have a different length,
#                 bigger than or equal to `eps`, but smaller than `dt + eps`.
#     :return: (ys, ts)
#              Where ys: array of states y(t_i), ts: array of times t_i
#     """
#     # Creating in array, which contains the times at which a function value
#     # shall be calculated:
#     N = len(y0)  # number of oscillators
#     noise_coeff = (2 * D * dt) ** (1 / 2)
#
#     ts = sp.arange(t_span[0], t_span[1], dt)
#     if ts[-1] > (t_span[1] - eps):
#         ts[-1] = t_span[1]
#     else:
#         ts = sp.append(ts, t_span[1])
#
#     # first an empty array is created to later replace the entries by the
#     # calculated function values of the numeric solution:
#     ys = sp.zeros((len(ts), N))
#     ys[0] = y0  # setting y(t0) to y0
#     # using the Euler method numeric function values are calculated and enterd
#     # in the array of y-values:
#     y = sp.array(y0)
#     for i, t in enumerate(ts[:-1]):
#         y += dt * fun(t, ys[i]) + noise_coeff * gaussian(N)
#         ys[i + 1] = y
#         # Or maybe direct is faster?
#         # ys[i + 1] = ys[i] + dt * fun(t, ys[i])
#         #
#     return (ys, ts)


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    N = 5
    y0 = sp.zeros(N)

    fun = lambda t, y: sp.zeros(N)  # no coupling

    D =1
    dt = 0.1
    tspan = (0, 3)

    ys, ts=  integrate_euler(y0, fun, D, dt, tspan)

    plt.plot(ts, ys)
    plt.show()