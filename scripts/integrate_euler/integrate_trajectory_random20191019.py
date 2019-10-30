'''
Integrate with Euler method with noise

- Random I.Co

INPUT: irun, ncycle, D, dt
'''

# import packages needed below
import  sys, os
import logging
import pickle
import scipy as sp

import carpet
import carpet.lattice_triangular as lattice

carpet.setup_logging('integrate_trajectory.log')
## Parameters
# Physics
set_name = 'machemer_1' # which hydrodynamic coefficients to use
order_g11 = (8,0)
order_g12 = (4,4)
period = 31.25 # [ms] period of cilia beat; freq = 2 * sp.pi / period [rad/ms]

# Geometry
nx = 6
ny = 6 # even number
N = nx * ny
a = 18  # [um] lattice spacing

## Initialize
# Geometry
L1,L2 = lattice.get_domain_sizes(nx,ny ,a)
coords, lattice_ids = lattice.get_nodes_and_ids(nx, ny, a)
N1, T1 = lattice.get_neighbours_list(coords, nx, ny, a)
get_k = lattice.define_get_k(nx, ny, a)
get_mtwist = lattice.define_get_mtwist(coords, nx, ny, a)



# Physics: carpet
gmat_glob, q_glob = lattice.define_gmat_glob_and_q_glob(set_name, a, N1, T1,order_g11, order_g12, period)
right_side_of_ODE = lattice.define_right_side_of_ODE(gmat_glob, q_glob)
solve_cycle = carpet.define_solve_cycle(right_side_of_ODE,2 * period, carpet.get_mean_phase)

def gaussian():  # returns random values, distributed as gaussian
    return sp.randn(N)  # with the same dimension as number of cilia


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
    t = t_span[0]
    t_end = t_span[1]
    y = sp.array(y0)
    noise_coeff = (2 * D * dt) ** (1 / 2)

    ys = [y]
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
    dy = fun(t, y) * noise_coeff * gaussian()
    y = y + dy
    t += dt

    ys.append(y)
    ts.append(t)

    return sp.array(ys), sp.array(ts)


def integrate_cycles(y0, D, dt, T, ncycle, eps):
    '''
    Call integrate Euler N times.
    Use right_side_of_ODE as input function.
    Reset phase to (0, 2 pi) interval after every cycle.
    Save initial state and after every cycle.
    '''
    y = y0
    t = 0
    ys_coarse = [y]
    ts_coarse = [t]
    for icycle in range(ncycle):
        ys, ts = integrate_euler(y, right_side_of_ODE, D, dt, (t, t + T), eps)
        y = ys[-1] % (2 * sp.pi)
        t = ts[-1]
        ys_coarse.append(y)
        ts_coarse.append(t)
    return ys_coarse, ts_coarse



def get_phi_random():
    '''
    Phases from 0 to 2pi & subtract mean phase so that we project initial condition to the Poincare plane.
    '''
    phi = 2 * sp.pi * sp.rand(N)
    return phi

## Prepare input
irun, ncycle, D, dt = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])

output_folder = "traj_random_D={:.3E}_dt={:.3E}/".format(D, dt)
os.makedirs(output_folder, exist_ok=True)
output_name = 'phi_{}.pkl'.format(irun)

## Save rng state
state = sp.random.get_state()
with open(output_folder + 'state_{}.pkl'.format(irun), 'wb') as f:
    pickle.dump(state, f, pickle.HIGHEST_PROTOCOL)

## Run simulation
phi0 = get_phi_random()
phis, ts = integrate_cycles(phi0, D, dt, period, ncycle, eps=10 ** -3 * dt)

with open(output_folder + output_name, 'wb') as f:
    pickle.dump(phis, f, pickle.HIGHEST_PROTOCOL)

logging.info("Finished run {} at random I.Co.; D={:.3E}; dt={:.3E}".format(irun, D, dt))