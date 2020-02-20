'''
2019-07-11:
    - Rotated lattice
    - m-twists used instead of fixed points


- Load a batch; continue trajectories which have not converged to a vicinity of a fixed point
  or until `ncycle` is reached (counting starts from 0 again)
- Create a new output batch
'''
import sys
import pickle
import scipy as sp
import pandas as pd

import carpet
import carpet.lattice_triangular2 as lattice
from carpet import circstd # root mean square
carpet.setup_logging('attraction_basins.log')

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

# Fixpoints


def get_fixpoint(k1,k2): # Create a copy of a fixpoint - to avoid editing
    fp = get_mtwist(k1,k2)
    fp -= carpet.get_mean_phase(fp)
    return fp

# Attraction basins definitions
"""
2019-09-11: 
- check that dphi decays

"""

def circ_dist(x, y):
    return circstd(x - y)


def solve_many_cycles(phi0, tol, ncycle, termination_eps):
    '''
    Returns trajectory: phis: [phi0, phi1, phi2,..]
    Terminate when the distance travelled in one cycle less is than `termination_eps`, or when `ncycle` is reached.
    '''

    def termination_condition(dphi, dphi_prev):
        ddphi = dphi - dphi_prev
        return (ddphi < 0 and dphi < termination_eps)  # termination condition

    phis = [phi0]
    dphi_prev = 0
    for icycle in range(ncycle):
        sol = solve_cycle(phi0, tol)
        phi1 = sol.y.T[-1] % (2 * sp.pi)
        phis.append(phi1)

        dphi = circ_dist(phi1, phi0)
        if termination_condition(dphi, dphi_prev):
            break
        else:
            phi0 = phi1
            dphi_prev = dphi
    return phis  # length = num of cycles run + 1


def get_phi_random():
    '''
    Phases from 0 to 2pi & subtract mean phase so that we project initial condition to the Poincare plane.
    '''
    phi = 2 * sp.pi * sp.rand(N)
    phi -= carpet.get_mean_phase(phi)
    return phi


def find_closest_fixpoint(phi):
    k1_min, k2_min, dist_min = None, None, sp.inf
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
        phis = solve_many_cycles(phi0, tol, ncycle, termination_eps)
        phi_end = phis[-1]
        dphi = circ_dist(phi_end, phis[-2])

        (m1, m2), dist = find_closest_fixpoint(phi_end)

        if dist < fixpoint_radius:
            conv = (m1, m2)
        else:
            m1 = sp.nan
            m2 = sp.nan
            conv = sp.nan

        res_df = res_df.append(dict(phi0=row['phi0'], ncycle=row['ncycle'] + len(phis) - 1,
                                    phi1=phi_end, dist=dist, k1=m1, k2=m2, conv=conv,
                                    dphi=dphi), ignore_index=True)
    res_df.to_pickle(output_name)

## Prepare input
ibatch, ncycle, tol,  fixpoint_radius, termination_eps = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]),\
                                                           float(sys.argv[4]), float(sys.argv[5])
input_name = 'res_{}.pkl'.format(ibatch)
output_name = 'res2_{}.pkl'.format(ibatch)

## Get input
with open(input_name, 'rb') as f:
    input_df = pickle.load(f)

## Run simulation
integrate_trajectories_more(input_df, tol, ncycle, fixpoint_radius, termination_eps, output_name)