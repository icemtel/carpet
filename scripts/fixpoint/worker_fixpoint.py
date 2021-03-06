'''
Search for fixed points in vicinity of a plane wave.

args: k1, k2 - wave numbers

===Optimization method discussion===
Use now - `lm` with non-zero starting phase; If it fails - switch to other methods

`hybr` -  NOT OK: strange jump and breaks my solve_cycle
'lm' - OK if start with non-zero mean phase
'broyden1' - was OK, but need more testing

MAYBE: try scipy.optimize.fixed_point
MAYBE: try `krylov` - it is said to be good for large-scale problems.
'''

import sys
import os
import pickle
import numpy as np
import logging
from sim_physics import carpet, solve_cycle, N, NN, TT, get_mtwist, define_solve_cycle_class
import carpet.classes as cc
import scipy.optimize as opt

carpet.setup_logging('master.log', mode='a')

def get_vec_diff(solve_cycle, tol):
    def vec_diff(phi):
        #        phi_global_end = 0
        #        sol = solve_cycle(phi, tol, phi_global_end=phi_global_end)  # end at global phase = 0
        sol = solve_cycle(phi, tol)  # end at global phase = 0

        phi1 = sol.y.T[-1]

        diff = phi1 - phi - 2 * np.pi - phi.mean()  # force it to prefer zero mean phase
        return diff

    return vec_diff


def find_fixpoint(phi0, tol, mini_tol):
    '''
    v2: optional phi0 as input; otherwise use m-twist
    2019-08-12: - subtract mean phase from the mtwist to work in the same Poincare plane
                - change dynamics for sine coupling
    '''
    # Map to classes
    ix_to_class, class_to_ix = cc.get_classes(phi0)  # cc.get_classes(phi0)
    nclass = len(class_to_ix)
    # Get classes representatives
    # Get one oscillator from each of cilia classes
    unique_cilia_ids = np.array([class_to_ix[iclass][0] for iclass in range(nclass)], dtype=np.int64)
    # Get neighbours
    NN_class, TT_class = cc.get_neighbours_list_class(unique_cilia_ids, ix_to_class, NN, TT)
    solve_cycle_class = define_solve_cycle_class(NN_class,TT_class)

    # Optimize!
    phi0_class = phi0[unique_cilia_ids]  # I.Co.
    vec_diff = get_vec_diff(solve_cycle_class, tol)
    res = opt.root(vec_diff, phi0_class, tol=mini_tol, method='lm')

    if not res.success:
        logging.warning(f'Did not converge, k1,k2=({k1},{k2})')

    fixpoint_class = np.array(res.x)
    fixpoint = fixpoint_class[ix_to_class]  # from classes to cilia ix

    return fixpoint

### Main
k1, k2 = int(sys.argv[1]), int(sys.argv[2])
tol = 10 ** -8

phi0 = get_mtwist(k1, k2)
# phi0 = phi0  - carpet.get_mean_phase(phi0)  # Test without this line - it makes it worse at least in 1 case
phi0 = phi0  - carpet.get_mean_phase(phi0) + 0.01 # Start with SMALL mean phase -> clear direction of minimizing for the solver

if k1 == 0 and k2 == 0:
    # assuming that symmetries of fixpoint are preserved,
    # there is only one possibility up to a constant shift: (0,0,0..,0)
    fixpoint = get_mtwist(0,0)
else:
    fixpoint = find_fixpoint(phi0, tol, tol)

## Test
# sol = solve_cycle(fixpoint, tol)
# fp_image = sol.y.T[-1] - 2 * np.pi
# print("fp mean phase:", fixpoint.mean())
# print("dist to image:", carpet.rms(fp_image - fixpoint))

outfolder = 'out/fixpoint/'
os.makedirs(outfolder, exist_ok=True)
filename = outfolder  + "fixpoint_k1={}_k2={}.npy".format(k1, k2)
np.save(filename, fixpoint)
