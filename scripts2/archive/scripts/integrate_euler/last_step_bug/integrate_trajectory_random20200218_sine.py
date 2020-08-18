'''
Integrate with Euler method with noise

- Random I.Co

INPUT: irun, ncycle, D, dt
'''

# import packages needed below
import sys, os
import logging
import pickle
import scipy as sp
import numpy as np

import carpet
import carpet.lattice.triangular as lattice
import carpet.physics.kuramoto as physics

carpet.setup_logging('integrate_trajectory.log')

## Parameters
# Physics
period = 31.25  # [ms] period of cilia beat; freq = 2 * sp.pi / period [rad/ms]
freq = 2 * np.pi / period
# Geometry
nx = 6
ny = 6  # even number
N = nx * ny
a = 18  # [um] lattice spacing

sin_str = 0.003 * freq

## Initialize
# Geometry
L1, L2 = lattice.get_domain_sizes(nx, ny, a)
coords, lattice_ids = lattice.get_nodes_and_ids(nx, ny, a)
N1, T1 = lattice.get_neighbours_list(coords, nx, ny, a)
get_k = lattice.define_get_k(nx, ny, a)
get_mtwist = lattice.define_get_mtwist(coords, nx, ny, a)

# Physics: sine-coupling
coupling = physics.define_sine_coupling(sin_str)
right_side_of_ODE = physics.define_right_side_of_ODE(coupling, freq, N1, T1)
solve_cycle = carpet.define_solve_cycle(right_side_of_ODE, 2 * period, carpet.get_mean_phase)


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
        ys, ts = carpet.integrate_euler(y, right_side_of_ODE, D, dt, (t, t + T), eps)
        y = ys[-1] % (2 * np.pi)
        t = ts[-1]
        ys_coarse.append(y)
        ts_coarse.append(t)
    return ys_coarse, ts_coarse


def get_phi_random():
    '''
    Phases from 0 to 2pi & subtract mean phase so that we project initial condition to the Poincare plane.
    '''
    phi = 2 * np.pi * np.random.rand(N)
    return phi


## Prepare input
irun, ncycle, D, dt = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])

output_folder = "traj_random_D={:.3E}_dt={:.3E}/".format(D, dt)
os.makedirs(output_folder, exist_ok=True)
output_name = 'phi_{}.pkl'.format(irun)

## Save rng state
state = np.random.get_state()
with open(output_folder + 'state_{}.pkl'.format(irun), 'wb') as f:
    pickle.dump(state, f, pickle.HIGHEST_PROTOCOL)

## Run simulation
phi0 = get_phi_random()
phis, ts = integrate_cycles(phi0, D, dt, period, ncycle, eps=10 ** -3 * dt)

with open(output_folder + output_name, 'wb') as f:
    pickle.dump(phis, f, pickle.HIGHEST_PROTOCOL)

# logging.info("Finished run {} at random I.Co.; D={:.3E}; dt={:.3E}".format(irun, D, dt))
