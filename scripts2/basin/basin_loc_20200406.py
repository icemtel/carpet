'''
- Small perturbation of fixed point
  - rms(dphi) = eps (on a (n-1)-sphere)
  - mean(dphi) = 0 (on Poincare section)
  - uniformly distributed
- Computed trajectory

2020-03-31: fixed saving interval bug
2020-04-06: !!! perturbation on a hypercube

args: irun, k1,k2, eps, ncycle, tol, save_every
'''
import sys, os
import numpy as np
import pickle
from sim_setup import N, solve_cycle, get_fixpoint
import carpet
from carpet.various import phases_to_interval


carpet.setup_logging('basin_loc.log')

## Define a perturbation - on a certain distance from the fixed point
def gram_schmidt_columns(X):
    Q, R = np.linalg.qr(X)
    return Q


def rand_cube(n):  # get random n-dimensional vector; distributed uniformly on surface of n-1-dim hypercube
    # Cube <=> one component = +-1; others take any value from -1 to 1 (distributed uniformly)
    u = 1 - 2 * np.random.rand(n)  # unifrom from -1 to 1
    ifacet = np.random.randint(0, n) # which index will take value from -1 or +1
    if u[ifacet] > 0: # determine sign
        sign = 1
    else:
        sign = -1
    u[ifacet] = sign
    return u


def get_linear_transformation(n):
    X = np.identity(n)
    X[:, 0] = np.ones(n)
    R = gram_schmidt_columns(X)
    return R


def define_get_pert(n):
    R = get_linear_transformation(n)

    def get_pert(eps):
        v = rand_cube(n - 1)  # get (uniformly) random vector on n-2 dimensional sphere
        v = np.concatenate(([0], v))  # embed the cube in the (n-1)-hyperplane in n-dimensional space
        return eps * R.dot(v)  * N ** (1 / 2)   # sqrt(N) to be consistent: before I used root mean square as a norm

    return get_pert

def solve_many_cycles_return_all(phi0, tol, ncycle, save_every):
    '''
    Returns trajectory: phis: [phi0, phi1, phi2,..]
    Saves only every n-th cycle, where n=:save_every:
    '''
    phis = [phi0]
    ts = [0]
    t = 0
    save_counter = 1
    for icycle in range(ncycle):
        sol = solve_cycle(phi0, tol)
        phi1 = sol.y.T[-1] - 2 * np.pi
        phi1 = phases_to_interval(phi1)
        t += sol.t[-1]

        phi0 = phi1

        if save_counter == save_every:
            phis.append(phi1)
            ts.append(t)
            save_counter = 1
        else:
            save_counter += 1

    return phis, ts


## Prepare input
irun, k1, k2, eps, ncycle, tol, save_every = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), \
                                             float(sys.argv[4]), int(sys.argv[5]), float(sys.argv[6]), int(sys.argv[7])

outfolder = f'basin_loc_k1={k1}_k2={k2}_eps={eps:.3G}/'
os.makedirs(outfolder, exist_ok=True)

get_pert = define_get_pert(N)
dphi0 = get_pert(eps)
phi0 = get_fixpoint(k1, k2) + dphi0
phis, ts = solve_many_cycles_return_all(phi0, tol, ncycle, save_every)

with open(outfolder + 'phi_{}.pkl'.format(irun), 'wb') as f:
    pickle.dump(phis, f, pickle.HIGHEST_PROTOCOL)

with open(outfolder + 'ts_{}.pkl'.format(irun), 'wb') as f:
    pickle.dump(ts, f, pickle.HIGHEST_PROTOCOL)
