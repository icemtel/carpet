import scipy as sp
from scipy.integrate import solve_ivp
import scipy.stats as stats


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

    def solve_cycle(phi_init, tol, phi_global_end=None, max_step=t_max / 20, ncycle=1):
        """
        :param phi_global_end: if not given, the cycle will finish at the initial phase, otherwise at phase_finish up to 2pi
        :param ncycle: solve several cycles; default: 1 cycle
               Works only if global phase increases over 2pi rather than  resetting to zero
               ex: mean phase works; circular mean phase - doesn't
        :return:
        """
        if phi_global_end is None:
            phi_global_end = phi_global_func(phi_init)  # end at starting phase (up to 2pi)

        def end_cycle_event(t, phi):
            '''
            Event triggered when this function returns zero.
            Defined in such a way that it returns zero when global phase increment = 2 * sp.pi * ncycle
            '''
            glob_phase_increment = (phi_global_func(phi) - phi_global_end) * increment_sign
            if glob_phase_increment > 2 * sp.pi * (ncycle - 0.5):
                return sp.sin(glob_phase_increment / 2)
            else:
                return sp.inf

        if backwards:
            increment_sign = -1
            right_side = lambda t, phi: - right_side_of_ODE(t, phi)
        else:
            increment_sign = +1
            right_side = right_side_of_ODE

        end_cycle_event.direction = -1  # event only triggered if return variable passes through zero from positive to negative values
        end_cycle_event.terminal = True  # tell solver to terminate the process in case of the event

        # Local error estimates are kept less than `atol + rtol * abs(y)`
        atol = tol  # absolute tolerance
        rtol = 3 * 10 ** -14  # = 0, but scipy doesn't like zero relative tolerance

        t_span = (0, ncycle * t_max)

        sol = solve_ivp(right_side, t_span, phi_init, rtol=rtol, atol=atol, max_step=max_step,
                        events=end_cycle_event)  # returns a solution class

        # Check that we ended cycle;
        # if not - the list of times when event was triggered will be empty and we raise an error
        if sol.t_events[0].size == 0:
            raise RuntimeError("solve_cycle: end of cycle event was not triggered")

        return sol

    return solve_cycle




## Global phase
### Circular Average
from scipy.stats import circmean, circstd, circvar


def circ_dist(phi, phi0, axis=None):
    '''
    If two phase vectors age given, calculate for dphi=phi-phi0
    '''
    return stats.circstd(phi - phi0, axis=axis)


def complex_exp_mean(phi, phi0=0):
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


def order_parameter(phi, phi0=0):
    return abs(complex_exp_mean(phi, phi0))


### Mean phase
def get_mean_phase(phi): # keep for compatibility
    return sp.mean(phi)

def mean_phase(phi):
    return sp.mean(phi)