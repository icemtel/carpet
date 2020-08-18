'''
Integrate with Euler method with noise

- Random I.Co

INPUT: irun, ncycle, D, dt, save_every

2020-02-19: upd scipy->numpy; upd carpet
2020-03-23: new argument: save_every -> save every n-th cycle only
'''

# import packages needed below
import  sys, os
import pickle
import numpy as np

import carpet
import carpet.lattice.triangular2 as lattice
import carpet.physics.friction_pairwise as physics

carpet.setup_logging('integrate_trajectory.log')

## Parameters
# Physics
set_name = 'machemer_1' # which hydrodynamic coefficients to use
order_g11 = (8,0)
order_g12 = (4,4)
period = 31.25            # [ms] period of cilia beat
freq = 2 * np.pi / period # [rad/ms] angular frequency

# Geometry
nx = 12
ny = 12  # even number
N = nx * ny
a = 18  # [um] lattice spacing

## Initialize
# Geometry

L1,L2 = lattice.get_domain_sizes(nx,ny ,a)
coords, lattice_ids = lattice.get_nodes_and_ids(nx, ny, a)  # get cilia (nodes) coordinates

distances = [1] #
NN, TT = lattice.get_neighbours_list(coords, nx, ny, a, distances)
e1, e2 = lattice.get_basis()
get_k = lattice.define_get_k(nx, ny, a)
get_mtwist = lattice.define_get_mtwist(coords, nx, ny, a)

# Physics
gmat_glob, q_glob = physics.define_gmat_glob_and_q_glob(set_name, e1, e2, a, NN, TT, order_g11, order_g12, period)
right_side_of_ODE = physics.define_right_side_of_ODE(gmat_glob, q_glob)
solve_cycle = carpet.define_solve_cycle(right_side_of_ODE, 2 * period, phi_global_func=carpet.get_mean_phase)


def integrate_cycles(y0, D, dt, T, ncycle, eps, save_every=1):
    '''
    Call integrate Euler N times.
    Use right_side_of_ODE as input function.
    Reset phase to (0, 2 pi) interval after every cycle.
    Save initial state and after every :save_every:-th cycle
    '''
    y = y0
    t = 0
    ys_coarse = [y]
    ts_coarse = [t]
    save_counter = 1
    for icycle in range(ncycle):
        ys, ts = carpet.integrate_euler(y, right_side_of_ODE, D, dt, (t, t + T), eps)
        y = ys[-1] % (2 * np.pi)
        t = ts[-1]

        if save_counter == save_every:
            save_counter = 1
            ys_coarse.append(y)
            ts_coarse.append(t)
        else:
            save_counter += 1

    return ys_coarse, ts_coarse


## Prepare input
irun, ncycle, D, dt, save_every = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), \
                                  float(sys.argv[4]), int(sys.argv[5])

output_folder = "traj_random_D={:.3E}_dt={:.3E}/".format(D, dt)
os.makedirs(output_folder, exist_ok=True)
output_name = 'phi_{}.pkl'.format(irun)

## Save rng state
state = np.random.get_state()
with open(output_folder + 'state_{}.pkl'.format(irun), 'wb') as f:
    pickle.dump(state, f, pickle.HIGHEST_PROTOCOL)

## Run simulation
phi0 = carpet.get_phase_random(N)
phis, ts = integrate_cycles(phi0, D, dt, period, ncycle, eps=10 ** -3 * dt, save_every=save_every)

with open(output_folder + output_name, 'wb') as f:
    pickle.dump(phis, f, pickle.HIGHEST_PROTOCOL)

# logging.info("Finished run {} at random I.Co.; D={:.3E}; dt={:.3E}".format(irun, D, dt))