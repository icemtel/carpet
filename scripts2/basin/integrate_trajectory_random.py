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

args: irun, ncycle, tol, save_every
'''
import sys
import pickle
import numpy as np

import carpet
from carpet.various import phases_to_interval
from sim_setup import solve_cycle, N


logger = carpet.setup_logging('attraction_basin.log')


def solve_cycles_many(phi0, tol, ncycle, save_every):
    '''
    Returns trajectory: phis: [phi0, phi1, phi2,..]
    Terminate when the distance travelled in one cycle less is than `termination_eps`, or when `ncycle` is reached.
    '''
    phis = [phi0]
    ts = [0]
    t = 0
    save_counter = 0

    for icycle in range(ncycle):
        sol = solve_cycle(phi0, tol, ncycle=1) # MAYBE: set to ncycle=save_every -> need to check precision
        phi1 = sol.y.T[-1] - 2 * np.pi
        phi1 = phases_to_interval(phi1)
        t += sol.t[-1]
        phi0 = phi1.copy()

        save_counter += 1
        if save_counter == save_every:
            phis.append(phi1)
            ts.append(t)
            save_counter = 0

    return np.array(phis), np.array(ts)


# def solve_cycles_many(phi0, tol, ncycle, save_every):
#     '''
#     Returns trajectory: phis: [phi0, phi1, phi2,..]
#     Terminate when the distance travelled in one cycle less is than `termination_eps`, or when `ncycle` is reached.
#     '''
#     phis = [phi0]
#     ts = [0]
#     t = 0
#
#     for icycle in range(ncycle // save_every):
#         sol = solve_cycle(phi0, tol, ncycle=save_every) # MAYBE: set to ncycle=save_every -> need to check precision
#         phi1 = sol.y.T[-1] - 2 * save_every * np.pi
#         phi1 = phases_to_interval(phi1)
#         t += sol.t[-1]
#         phi0 = phi1.copy()
#
#         phis.append(phi1)
#         ts.append(t)
#
#     return np.array(phis), np.array(ts)



# Run tests if script is run individually

## Prepare input
irun, ncycle, tol, save_every = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])

# Save random state
state = np.random.get_state()
with open('obj/basin/state_{}.pkl'.format(irun), 'wb') as f:
    pickle.dump(state, f, pickle.HIGHEST_PROTOCOL)
# Intial condition
phi0 = carpet.get_phase_random(N)
np.save('obj/basin/phi0_{}.npy'.format(irun), phi0)
# Run
phis, ts = solve_cycles_many(phi0, tol, ncycle, save_every)

np.save('out/basin/phi_{}.npy'.format(irun), phis)
np.save('out/basin/ts_{}.npy'.format(irun), ts)
