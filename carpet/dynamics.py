import scipy as sp
from scipy.integrate import solve_ivp
import scipy.stats as stats


def define_solve_cycle(right_side_of_ODE, t_max, phi_global_func, backwards=False):
    '''
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

    def solve_cycle(phi_init, tol, phi_global_end=None, max_step=t_max / 20):
        """
        :param phi_global_end: if not given, the cycle will finish at the initial phase, otherwise at phase_finish up to 2pi
        :return:
        """
        if phi_global_end is None:
            phi_global_end = phi_global_func(phi_init)  # end at starting phase (up to 2pi)

        def end_cycle_event(t, phi):
            '''
            Event triggered when this function returns zero.
            By definition, it returns zero when global phase gets increment of 2pi.
            '''
            if t > 0:
                return sp.sin((phi_global_func(phi) - phi_global_end) / 2)
            else:
                return sp.inf

        if backwards:
            right_side = lambda t, phi: - right_side_of_ODE(t, phi)
            end_cycle_event.direction = +1
        else:
            right_side = right_side_of_ODE
            end_cycle_event.direction = -1  # event only triggered if return variable passes through zero from positive to negative values

        end_cycle_event.terminal = True  # tell solver to terminate the process in case of the event

        # Local error estimates are kept less than `atol + rtol * abs(y)`
        atol = tol / 2  # absolute tolerance
        rtol = atol / 2 / sp.pi  # corresponding relative tolerance (since phi is bounded by 2 pi)

        t_span = (0, t_max)

        return solve_ivp(right_side, t_span, phi_init, rtol=rtol, atol=atol, max_step=max_step,
                         events=end_cycle_event)  # returns a solution class

    return solve_cycle


## Global phase
### Circular Average

def circmean(phi, phi0=0, high=2 * sp.pi, low=0, axis=None):
    '''
    If two phase vectors age given, calculate for dphi=phi-phi0
    '''
    return stats.circmean(phi - phi0, high=high, low=low, axis=axis)


def circstd(phi, phi0=0, high=2 * sp.pi, low=0, axis=None):
    '''
    If two phase vectors age given, calculate for dphi=phi-phi0
    '''
    return stats.circstd(phi - phi0, high=high, low=low, axis=axis)


def circvar(phi, phi0=0, high=2 * sp.pi, low=0, axis=None):
    '''
    If two phase vectors age given, calculate for dphi=phi-phi0
    '''
    return stats.circvar(phi - phi0, high=high, low=low, axis=axis)


def complex_exp_mean(phi, phi0):
    '''
    Mean of complex exponents of array of phases phi, relative to phi0
    '''
    return sp.mean(sp.exp(1j * (phi - phi0)))


def define_circular_mean_phase(phi0):
    '''
    :return: function which computes circular mean phase of phi, relative to phi0.
    '''

    def get_phi_global(phi):  # equivalent up to mod 2pi to sp.angle(sp.mean(sp.exp(1j * (phi - phi0))))
        return circmean(phi - phi0)

    return get_phi_global


def get_order_parameter(phi, phi0):
    return abs(complex_exp_mean(phi, phi0))


### Mean phase
def get_mean_phase(phi):
    return sp.mean(phi)
