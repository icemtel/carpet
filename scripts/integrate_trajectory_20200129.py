"""
Integrate for a number of cycles (without phase noise)

Save phases after each cycle in a pickle.


To use in try16c_debug
"""

import sys
import pickle
import scipy as sp
import pandas as pd


import carpet
import carpet.lattice_triangular2 as lattice
from carpet import circ_dist


def phases_to_interval(phi):
    '''
    Take phase vector.
    Return equivalent phase vector with
    phases in the +-2pi interval, centered at the mean of the initial phase vector (mean unchanged)
    '''
    x = sp.array(phi)
    xmin = x.mean() - 2 * sp.pi
    xmax = x.mean() + 2 * sp.pi
    flag = (x.max() > xmax) or (x.min() < xmin)
    while flag:
        imax = x.argmax()
        imin = x.argmin()
        x[imax] -= 2 * sp.pi
        x[imin] += 2 * sp.pi

        flag = (x.max() > xmax) or (x.min() < xmin)

    return x



def integrate_cycles(y0, solve_cycle, tol, ncycle):
    '''
    Call integrate solve_cycle N times.
    Reset phase to a fixed interval after every cycle.
    Save initial state and after every cycle.
    '''
    y = y0
    t = 0
    ys_coarse = [y]
    ts_coarse = [t]
    for icycle in range(ncycle):
        sol = solve_cycle(y, tol)
        y = sol.y.T[-1]
        y = phases_to_interval(y)
        t = sol.t[-1]
        ys_coarse.append(y)
        ts_coarse.append(t)
    return ys_coarse, ts_coarse


# to use in notebook
tol = 1.e-6
ncycle =  200
ntraj = 10
res_path = objfolder + '/traj'

for ix in tqdm(range(ntraj)):
    phi0 = res_df.at[ix, 'phi1']
    phis, ts = integrate_cycles(phi0, solve_cycle, tol, ncycle)

    dump_object(phis, f'phi_{ix}.pkl',path=res_path)
    dump_object(ts, f'T_{ix}.pkl',path=res_path)