'''
- Load a batch; continue trajectories which have not converged to a vicinity of a fixed point.
- Create a new output file


Prolong trajectories from a specific pickled dataframe file

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


def integrate_trajectories_more(input_df, tol, ncycle, fixpoint_radius, termination_eps, output_name):
    """
    For parallel work: gets trajectories - computes dataframe and saves it as a pickle
    """
    res_df = initialize_results_df()

    for ix in input_df.index:
        row = input_df.loc[ix]
        phi0 = row['phi1']
        # if converged, just resave the result
        if not pd.isnull(row['conv']):
            res_df.loc[ix]= row.copy()
            ## Fix dphi if necessary
            if pd.isnull(row['dphi']):
                sol = solve_cycle(phi0, tol)
                phi1 = sol.y.T[-1]
                dphi = circ_dist(phi1, phi0)
                res_df.loc[ix]['dphi'] = dphi
            continue
        # else solve more cycles to continue the trajectory
        phi_end, dphi, ncycle_extra = solve_many_cycles(phi0, tol, ncycle, termination_eps)

        (m1, m2), dist = find_closest_fixpoint(phi_end)

        if dist < fixpoint_radius:
            conv = (m1, m2)
        else:
            m1 = np.nan
            m2 = np.nan
            conv = np.nan

        res_df = res_df.append(dict(phi0=row['phi0'], ncycle=row['ncycle'] + ncycle_extra,
                                    phi1=phi_end, dist=dist, k1=m1, k2=m2, conv=conv,
                                    dphi=dphi), ignore_index=True)
    res_df.to_pickle(output_name)

## Prepare input
tol = 10 ** -6
ncycle = 3000
fixpoint_radius = 0.2
termination_eps = 1e-4 # terminate cycles if we are close enough to a fixpoint


input_name = 'not_converged.pkl'
output_name = 'not_converged_out.pkl'

## Get input
with open(input_name, 'rb') as f:
    input_df = pickle.load(f)

## Run simulation
integrate_trajectories_more(input_df, tol, ncycle, fixpoint_radius, termination_eps, output_name)