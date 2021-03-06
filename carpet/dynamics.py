import numpy as np
from scipy.integrate import solve_ivp

def define_solve_cycle(right_side_of_ODE, t_max, phi_global_func, backwards=False):
    '''
    2019-09-19: only absolute tolerance
                parameter ncycle - solve n cycles
                raise an error if end_event wasn't triggered
    Create a function which would solve
        phi' = f(t,phi)
    given phi(0)=phi0 and solver tolerance with a termination criterium:
       phi_global(phi) = phi_global(phi0) + 2pi

    :param right_side_of_ODE: f(t,y)
    :param t_max: another termination criterium: t=t_max;  also used to determine maximum time-step size;
                  for safety one can choose t_max equal to double beat period
    :param backwards: if True: solve backwards in time; assuming now that global phase will decay in time;
                      In output time is still > 0; but the solution behaves as if times has the opposite sign
    :return: solve_cycle function
    - atol and rtol are chosen with an assumption that phases are close to interval [0, 2pi]
    '''

    def solve_cycle(phi_init, tol, phi_global_end=None, ncycle=1, **kwargs):
        """
        - Before end of cycle event is triggered the following condition is checked
          `glob_phase_increment > 2 * np.pi * (ncycle - 0.5)`
          Therefore this function is not suitable for running a small fraction of a cycle. TODO: fix it?
          Tried: setting ncycle to 0.5 + eps, but it doesn't seem to work (might work for bigger eps)
        :param phi_global_end: if not given, the cycle will finish at the initial phase, otherwise at phase_finish up to 2pi
        :param ncycle: solve several cycles; default: 1 cycle
               Works only if global phase increases over 2pi rather than  resetting to zero
               ex: mean phase works; circular mean phase - doesn't
        :param kwargs: key-word arguments for solve_ivp;
                       parameters `rtol` and `atol` will override tolerance specified by `tol` in `solve_cycle`
        :return:
        """
        # ==============================================================
        # Determine when to end computation (end cycle event)
        if phi_global_end is None:
            phi_global_end = phi_global_func(phi_init)  # end at starting phase (up to 2pi)

        def end_cycle_event(t, phi):
            '''
            Event triggered when this function returns zero.
            Defined in such a way that it returns zero when global phase increment = 2 * np.pi * ncycle
            '''
            glob_phase_increment = (phi_global_func(phi) - phi_global_end) * increment_sign
            if glob_phase_increment > 2 * np.pi * (ncycle - 0.5):
                return np.sin(glob_phase_increment / 2)
            else:
                return np.inf

        if backwards:
            increment_sign = -1
            right_side = lambda t, phi: - right_side_of_ODE(t, phi)
        else:
            increment_sign = +1
            right_side = right_side_of_ODE

        end_cycle_event.direction = -1  # event only triggered if return variable passes through zero from positive to negative values
        end_cycle_event.terminal = True  # tell solver to terminate the process in case of the event

        # ==============================================================
        # Local error estimates are kept less than `atol + rtol * abs(y)`
        if 'atol' not in kwargs.keys():
            kwargs['atol'] = tol  # absolute tolerance
        if 'rtol' not in kwargs.keys():
            kwargs['rtol'] = 3 * 10 ** -14  # = 0, but scipy doesn't like zero relative tolerance
        if 'max_step' not in kwargs.keys():
            kwargs['max_step'] = t_max / 20  # max step size - to be safe

        t_span = (0, ncycle * t_max)

        sol = solve_ivp(right_side, t_span, phi_init, events=end_cycle_event, **kwargs)  # returns a solution class
        # Check that we ended cycle;
        # if not - the list of times when event was triggered will be empty and we raise an error
        if sol.t_events[0].size == 0:
            raise RuntimeError("solve_cycle: end of cycle event was not triggered")

        return sol

    return solve_cycle


## Solve with noise


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
             Where ys: list of states y(t_i), ts: list of times t_i
    """
    N = len(y0)

    def gaussian():  # returns random values, distributed normally
        return np.random.randn(N)  # with the same dimension as number of oscillators

    t = t_span[0]
    t_end = t_span[1]
    y = np.array(y0)
    noise_coeff = (2 * D * dt) ** (1 / 2)

    ys = [y] # MAYBE: Create an array
    ts = [t]
    while t < t_end - dt - eps:
        dy = fun(t, y) * dt + noise_coeff * gaussian()
        y = y + dy  # don't  use += on vectors!
        t += dt

        ys.append(y)
        ts.append(t)

    # The last step
    dt = t_end - t
    noise_coeff = (2 * D * dt) ** (1 / 2)
    dy = fun(t, y) * dt + noise_coeff * gaussian()
    y = y + dy
    t += dt

    ys.append(y)
    ts.append(t)

    return np.array(ys), np.array(ts)


