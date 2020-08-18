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

args: ibatch, nrun, ncycle, tol,  fixpoint_radius, termination_eps
'''
import sys
import pickle
import numpy as np
import pandas as pd

import carpet
import carpet.lattice.triangular as lattice
import carpet.physics.kuramoto as physics
from carpet import circ_dist
from carpet.various import phases_to_interval

logger = carpet.setup_logging('attraction_basins.log')

## Parameters
# Physics
period = 31.25  # [ms] period of cilia beat; freq = 2 * np.pi / period [rad/ms]
freq = 2 * np.pi / period
# Geometry
nx = 6
ny = 6  # even number
N = nx * ny
a = 18  # [um] lattice spacing

sin_str = 0.0016 * freq

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




# Fixpoints: use m-twist instead of fixed point - good enough for basins of attraction
def get_fixpoint(k1,k2):
    fp = get_mtwist(k1,k2)
    fp -= carpet.get_mean_phase(fp)
    fp = phases_to_interval(fp)
    return fp

# Attraction basins definitions

def solve_many_cycles(phi0, tol, ncycle, termination_eps):
    '''
    Returns trajectory: phis: [phi0, phi1, phi2,..]
    Terminate when the distance travelled in one cycle less is than `termination_eps`, or when `ncycle` is reached.
    '''

    def termination_condition(dphi, dphi_prev):
        ddphi = dphi - dphi_prev
        return (ddphi < 0 and dphi < termination_eps)  # termination condition

    dphi_prev = 0
    for icycle in range(ncycle):
        sol = solve_cycle(phi0, tol)
        phi1 = sol.y.T[-1] - 2 * np.pi
        phi1 = phases_to_interval(phi1)

        dphi = circ_dist(phi1, phi0)

        if termination_condition(dphi, dphi_prev):
            break
        else:
            phi0 = phi1
            dphi_prev = dphi

    ncycle_fact = icycle + 1
    return phi1, dphi, ncycle_fact  # length = num of cycles run + 1


def find_closest_fixpoint(phi):
    k1_min, k2_min, dist_min = None, None, np.inf
    for k1 in range(nx):
        for k2 in range(ny):
            fixpoint = get_fixpoint(k1, k2)
            dist = circ_dist(phi, fixpoint)
            if dist < dist_min:
                (k1_min, k2_min) = k1, k2
                dist_min = dist

    return (k1_min, k2_min), dist_min


def initialize_results_df():
    '''
    If dist is too big - k1,k2 -> None
    '''
    return pd.DataFrame(columns=['phi0', 'dist', 'k1', 'k2', 'conv', 'phi1', 'ncycle', 'dphi'])


def integrate_trajectories(phi0s, tol, ncycle, fixpoint_radius, termination_eps, output_name):
    """
    For parallel work: integrates trajectories - compiles dataframe and saves it as a pickle
    """
    res_df = initialize_results_df()

    for phi0 in phi0s:
        phi_end, dphi, ncycle_fact = solve_many_cycles(phi0, tol, ncycle, termination_eps)

        (m1, m2), dist = find_closest_fixpoint(phi_end)

        if dist < fixpoint_radius:
            conv = (m1, m2)
        else:
            m1 = np.nan
            m2 = np.nan
            conv = np.nan

        res_df = res_df.append(dict(phi0=phi0, ncycle=ncycle_fact, phi1=phi_end, dist=dist, k1=m1, k2=m2, conv=conv,
                                    dphi=dphi), ignore_index=True)
    res_df.to_pickle(output_name)

def main(ibatch, nrun, ncycle, tol, fixpoint_radius, termination_eps):
    output_name = 'res_{}.pkl'.format(ibatch)
    ## Save rng state
    state = np.random.get_state()
    with open('state_{}.pkl'.format(ibatch), 'wb') as f:
        pickle.dump(state, f, pickle.HIGHEST_PROTOCOL)

    ## Run simulation
    phi0s = [carpet.get_phase_random(N) for irun in range(nrun)]
    integrate_trajectories(phi0s, tol, ncycle, fixpoint_radius, termination_eps, output_name)


# Run tests if script is run individually

## Prepare input
ibatch, nrun, ncycle, tol, fixpoint_radius, termination_eps = int(sys.argv[1]), int(sys.argv[2]), \
                                                              int(sys.argv[3]), float(sys.argv[4]), \
                                                              float(sys.argv[5]), float(sys.argv[6])
# Run
main(ibatch, nrun, ncycle, tol, fixpoint_radius, termination_eps)