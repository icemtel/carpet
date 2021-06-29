'''
Integrates trajectory for many cycles
- tries to load previously computed cycles; starts from the last point available
- saved result in a separate file (e.g. phi_0_pt1.npy)
- finishes early if a trajectory almost converged to a fixed point

- IF integrated for many thousands of cycles - may want to uncomment  `phi0 = phases_to_interval(phi0)`
  after finishing each cycle

args: irun, ncycle, tol, save_every, sim_name
'''

# import packages needed below
import carpet
from carpet.various import phases_to_interval
from sim_physics import N, solve_cycle
import sys, os
import numpy as np
import logging

carpet.setup_logging('worker.log')



def solve_cycles_many(phi0, tol, ncycle, save_every, conv_eps):
    '''
    Returns trajectory: phis: [phi0, phi1, phi2,..]
    Terminate when the distance travelled in one cycle less is than `termination_eps`, or when `ncycle` is reached.
    '''
    # save_every must be a positive integer; save_every = 1 => every cycle saved; save_every=2 => every 2nd saved
    if save_every == 0:
        raise NotImplementedError

    phis = [phi0]
    dphis_norm = []
    ts = [0]
    t = 0
    save_counter = 0

    # phidots_mean = []
    # phidots_std = []
    for icycle in range(ncycle):
        sol = solve_cycle(phi0, tol, ncycle=1)
        phi1 = sol.y.T[-1] - 2 * np.pi
        t += sol.t[-1]

        # Save data once every `save_every` cycles
        save_counter += 1
        if save_counter == 1: # Add dphi for the first point & all points which got saved recently
            dphi = (phi1 - phi0)
            dphi_norm = np.sqrt(1 / N) * np.linalg.norm(dphi)
            dphis_norm.append(dphi_norm)
            # END if change in cycle is too small => (therefore close to a fixed point)
            if dphi_norm < conv_eps:
                return np.array(phis), np.array(ts), np.array(dphis_norm)
            # For small dphi; with zero mean phase; the norm above is equivalent to
            # `np.sqrt(1 - carpet.order_parameter(dphi) ** 2)`
        if save_counter == save_every:
            phis.append(phi1)
            ts.append(t)
            save_counter = 0 # reset save counter

        phi0 = phi1.copy()  # set initial condition for the next cycle
        # phi0 = phases_to_interval(phi0)

    return np.array(phis), np.array(ts), np.array(dphis_norm)



def get_traj_filename(irun, ipart, path):
    if ipart == 0:
        return os.path.join(path, f'phi_{irun}.npy')
    else:
        return os.path.join(path, f'phi_{irun}_pt{ipart}.npy')


def get_ts_filename(irun, ipart, path):
    if ipart == 0:
        return os.path.join(path, f'ts_{irun}.npy')
    else:
        return os.path.join(path, f'ts_{irun}_pt{ipart}.npy')


def load_phis(irun, path):
    '''
    Read previous parts of the trajectory
    '''
    phis_list = []
    for ipart in range(64):
        filename = get_traj_filename(irun, ipart, path)
        if os.path.isfile(filename):
            phis_pt = np.load(filename)
            phis_list.append(phis_pt)
        else:
            break
    return np.concatenate(phis_list)  # trajectory glued back from parts

## Prepare input
irun, ncycle_total, tol, save_every, sim_name = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), \
                                                int(sys.argv[4]), str(sys.argv[5])

# Folder names
objfolder = f'obj/{sim_name}/'
outfolder = f'out/{sim_name}/'
conv_eps = 0.99e-4

# Find how many parts the trajectory already has
ipart_last = None
# Find the last existing part of the trajectory
for i in range(64):  # maximum number of parts
    filename = get_traj_filename(irun, i, outfolder)
    if os.path.isfile(filename):
        ipart_last = i
    else:
        break

# If trajectory exists -> load, get initial condition and number of cycles
if ipart_last is not None:
    phis_old = load_phis(irun, outfolder)
    ncycle_old = (len(phis_old) - 1) * save_every  # assume that input save_every is the same as used in prev. sims!
    phi0 = phis_old[-1]
    ipart = ipart_last + 1
    del phis_old  # free up memory
else:
    ipart = 0
    ncycle_old = 0
    # Load input
    input_filename = objfolder + f'phi0_{irun}.npy'
    phi0 = np.load(input_filename)

## Run simulation
ncycle_extra = ncycle_total - ncycle_old
if ncycle_extra > 0:
    phis, ts = solve_cycles_many(phi0, tol, ncycle_extra, save_every, conv_eps)

    if ipart > 0:  # remove the first point because it's the same as the last point in the first part
        phis = phis[1:]
        ts = ts[1:]

    ## Save output
    if len(phis) > 1: # length = 1 if imeediately finished simulation AND part > 0
        os.makedirs(outfolder, exist_ok=True)
        # Mean and std frequency
        # Time points
        filename = get_ts_filename(irun, ipart, outfolder)
        np.save(filename, ts)
        # Phases - saved the last to make sure that everything else is saved as well
        filename = get_traj_filename(irun, ipart, outfolder)
        np.save(filename, phis)

