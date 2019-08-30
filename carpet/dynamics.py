import scipy as sp
from scipy.integrate import  solve_ivp


def define_solve_cycle(right_side_of_ODE, t_max, phi_global_func):
    '''
    Create a function which would solve
        phi' = f(t,phi)
    given phi(0)=phi0 and solver tolerance with a termination criterium:
       phi_global(phi) = phi_global(phi0) + 2pi

    :param right_side_of_ODE: f(t,y)
    :param t_max: another termination criterium: t=t_max;  also used to determine maximum time-step size;
                  for safety one can choose t_max equal to double beat period
    :return: solve_cycle function
    - atol and rtol are chosen with an assumption that phases are close to interval [0, 2pi]
    '''
    def solve_cycle(phi_init, tol, phi_global_end=None, max_step=t_max / 20):
        """
        :param phi_global_end: if not given, the cycle will finish at the initial phase, otherwise at phase_finish up to 2pi
        :return:
        """
        if phi_global_end is None:
            phi_global_end = phi_global_func(phi_init) # end at starting phase (up to 2pi)

        def end_cycle_event(t, phi):
            '''
            Event triggered when this function returns zero.
            By definition, it returns zero when global phase gets increment of 2pi.
            '''
            if t > 0:
                return 2 * sp.sin((phi_global_func(phi) - phi_global_end) / 2)
            else:
                return sp.inf

       # end_cycle_event.direction = +1  # event only triggered if return variable passes through zero from negative to positive values
        end_cycle_event.terminal = True  # tell solver to terminate the process in case of the event

        # Local error estimates are kept less than `atol + rtol * abs(y)`
        atol = tol / 2  # absolute tolerance
        rtol = atol / 2 / sp.pi  # corresponding relative tolerance (since phi is bounded by 2 pi)

        t_span = (0, t_max)

        return solve_ivp(right_side_of_ODE, t_span, phi_init, rtol=rtol, atol=atol, max_step=max_step,
                         events=end_cycle_event)  # returns a solution class


    return solve_cycle


## Global phase
### Circular Average
def get_circular_average(phi, phi0):
    return sp.mean(sp.exp(1j * (phi - phi0)))


def define_circular_average_phase(phi0):
    '''
    k1,k2 - m-twist parameters;
    Global phase;
    :return: a function which gives global phase
    '''

    def get_phi_global(phi):
        return sp.angle(get_circular_average(phi, phi0))

    return get_phi_global


def get_order_parameter(phi, phi0):
    return abs(get_circular_average(phi, phi0))


### Mean phase
def get_mean_phase(phi):
    return sp.mean(phi)




