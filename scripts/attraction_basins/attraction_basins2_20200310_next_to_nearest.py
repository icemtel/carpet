'''
- 12x12 lattice
- Use m-twist as the fixed point
- ! Don't subtract mean phase from I. Co.
- basic logging

2019-09-11: check that dphi decays
2019-11-25: Don't store phase vectors at intermediate time steps
2020-01-28: Redo analysis for the rotated lattice
2020-02-02: upd numpy and carpet

- To test: new procedure to keep zero mean phase

2020-03-10: Save every 10 cycles, return phis

args: irun, ncycle, tol
'''
import sys
import pickle
import numpy as np
import pandas as pd

import carpet
import carpet.lattice.triangular as lattice
import carpet.physics.friction_pairwise as physics
from carpet import circ_dist
from carpet.various import phases_to_interval

logger = carpet.setup_logging('attraction_basins.log')

## Parameters
# Physics
set_name = 'machemer_1' # hydrodynamic friction coefficients data set
period = 31.25 # [ms] period
freq = 2 * np.pi / period # [rad/ms] angular frequency
order_g11 = (8, 0) # order of Fourier expansion of friction coefficients
order_g12 = (4, 4)
# Geometry
a = 18  # [um]
nx = 6  # number of cilia in x-direction
ny = 6  # must be even
N = nx * ny

## Initialize
# Geometry
L1,L2 = lattice.get_domain_sizes(nx,ny ,a)
coords, lattice_ids = lattice.get_nodes_and_ids(nx, ny, a)  # get cilia (nodes) coordinates
N1, T1 = lattice.get_neighbours_list(coords, nx, ny, a, distances=(1, np.sqrt(3)))  # get list of neighbours and relative positions
e1, e2 = lattice.get_basis()
get_k = lattice.define_get_k(nx, ny, a)
get_mtwist = lattice.define_get_mtwist(coords, nx, ny, a)

# Physics
gmat_glob, q_glob = physics.define_gmat_glob_and_q_glob(set_name, e1, e2, a, N1, T1, order_g11, order_g12, period)
right_side_of_ODE = physics.define_right_side_of_ODE(gmat_glob, q_glob)
solve_cycle = carpet.define_solve_cycle(right_side_of_ODE, 2 * period, phi_global_func=carpet.get_mean_phase)


# Fixpoints: use m-twist instead of fixed point - good enough for basins of attraction
def get_fixpoint(k1,k2):
    fp = get_mtwist(k1,k2)
    fp -= carpet.get_mean_phase(fp)
    fp = phases_to_interval(fp)
    return fp

# Attraction basins definitions

def solve_many_cycles_return_all(phi0, tol, ncycle, save_interval):
    '''
    Returns trajectory: phis: [phi0, phi1, phi2,..]
    Terminate when the distance travelled in one cycle less is than `termination_eps`, or when `ncycle` is reached.
    '''
    dphi_prev = 0
    phis = [phi0]
    ts = [0]
    t = 0
    save_counter = 0
    for icycle in range(ncycle):
        sol = solve_cycle(phi0, tol)
        phi1 = sol.y.T[-1] - 2 * np.pi
        phi1 = phases_to_interval(phi1)
        t += sol.t[-1]

        phi0 = phi1

        save_counter += 1
        if save_counter == save_interval:
            phis.append(phi1)
            ts.append(t)
            save_counter = 0
    return phis,ts

def main(irun, ncycle, tol):
    ## Save rng state
    state = np.random.get_state()
    with open('state_{}.pkl'.format(irun), 'wb') as f:
        pickle.dump(state, f, pickle.HIGHEST_PROTOCOL)

    ## Run simulation
    phi0 = carpet.get_phase_random(N)
    phis,ts  = solve_many_cycles_return_all(phi0, tol, ncycle, save_interval)

    with open('phi_{}.pkl'.format(irun), 'wb') as f:
        pickle.dump(phis, f, pickle.HIGHEST_PROTOCOL)

    with open('ts_{}.pkl'.format(irun), 'wb') as f:
        pickle.dump(ts, f, pickle.HIGHEST_PROTOCOL)

# Run tests if script is run individually

## Prepare input
irun, ncycle, tol = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
save_interval = 10
# Run
main(irun, ncycle, tol)