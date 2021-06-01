'''
Integrate with Euler method with noise

- starts from an initial condition
       if not found raises an error
- If finds previous parts of trajectory
  Prolongs existing trajectories until some total number of cycles ncycle_total.
- `phases_to_interval` is used -> after each cycles phases are mapped to [-2pi,2pi] interval without changing mean phase


INPUT:
 - irun, ncycle_total, D, dt, save_every, sim_name
'''

# import packages needed below
import carpet
from carpet.various import phases_to_interval
from sim_physics import right_side_of_ODE, period, N
import sys, os
import pickle
import numpy as np


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
        y = ys[-1]
        y = phases_to_interval(y)  # phases to -2pi and 2pi; mean phase is preserved up to modulus 2pi
        t = ts[-1]

        if save_counter == save_every:
            save_counter = 1
            ys_coarse.append(y)
            ts_coarse.append(t)
        else:
            save_counter += 1

    return ys_coarse, ts_coarse

def get_traj_filename(irun, ipart, path):
    if ipart == 0:
        return os.path.join(path, f'phi_{irun}.npy')
    else:
        return os.path.join(path, f'phi_{irun}_pt{ipart}.npy')


def get_state_filename(irun, ipart, path):
    if ipart == 0:
        return os.path.join(path, f'state_{irun}.pkl')
    else:
        return os.path.join(path, f'state_{irun}_pt{ipart}.pkl')


def load_phis(irun, path):
    '''
    Read pickled file with phis array
    '''
    phis_list = []
    for ipart in range(64):
        filename = get_traj_filename(irun, ipart, path)
        if os.path.isfile(filename):
            phis_pt = np.load(filename)
            phis_list.append(phis_pt)
        else:
            break
    return np.concatenate(phis_list) # trajectory glued back from parts

## Load arguments
irun, ncycle_total, D, dt, save_every, sim_name = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), \
                                                  float(sys.argv[4]), int(sys.argv[5]), str(sys.argv[6])
# Folder names
objfolder = f'obj/{sim_name}/'
outfolder = f'out/{sim_name}/'


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
    ncycle_old = (len(phis_old) - 1) * save_every # assume that input save_every is the same as used in prev. sims!
    phi0 = phis_old[-1]
    ipart = ipart_last + 1
    del phis_old # free up memory
else:
    ipart = 0
    ncycle_old = 0
    # Load input
    input_filename = objfolder + f'phi0_{irun}.npy'
    phi0 = np.load(input_filename)


## Save rng state
dump_object = carpet.various.define_dump_object(objfolder)
state = np.random.get_state()
with open(get_state_filename(irun, ipart, objfolder), 'wb') as f:
    pickle.dump(state, f, pickle.HIGHEST_PROTOCOL)


## Run simulation
ncycle_extra = ncycle_total - ncycle_old
if ncycle_extra > 0:
    phis, ts = integrate_cycles(phi0, D, dt, period, ncycle_extra, eps=10 ** -3 * dt, save_every=save_every)

    if ipart > 0:  # remove the first point because it's the same as the last point in the first part
        phis = phis[1:]
    ## Save output
    os.makedirs(outfolder, exist_ok=True)
    filename = get_traj_filename(irun, ipart, outfolder)
    np.save(filename, phis)

# logging.info("Finished run {} at random I.Co.; D={:.3E}; dt={:.3E}".format(irun, D, dt))
