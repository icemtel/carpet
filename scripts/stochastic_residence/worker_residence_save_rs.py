'''
- Integrate trajectory - saved order parameters and the phase vectors. Optionally - with different frequencies,
    set by *save_every_r* and  *save_every_phi*:
 - Save order parameter every  *save_every_r* cycles
- Save phase vectors, every *save_every_phi* cycles
'''

# import packages needed below
import carpet
from sim_physics import right_side_of_ODE, period, N
import sys, os
import pickle
import numpy as np


def integrate_cycles(y0, D, dt, T, ncycle, eps, save_every_r=1, save_every_phi=1):
    '''
    Call integrate Euler N times.
    Use right_side_of_ODE as input function.
    Reset phase to (0, 2 pi) interval after every cycle.
    Save initial state and after every :save_every_phi:-th cycle
    Save initial order parameter and after every :save_every_r:-th cycle
    Output: phases modulo 2pi,  winding number increments, order_parameters
    '''
    order_parameter = carpet.order_parameter  # shorthand
    y = y0 % (2 * np.pi)  # initial phase modulo 2pi
    dw = np.zeros_like(y, dtype=np.int)  # will be filled with the winding number increment
    t = 0
    ys_coarse = [y]
    dws_coarse = []
    rs_coarse = [order_parameter(y)]
    save_counter_r = 1
    save_counter_phi = 1
    for icycle in range(ncycle):
        ys, ts = carpet.integrate_euler(y, right_side_of_ODE, D, dt, (t, t + T), eps)
        y = ys[-1]
        dw += np.array(np.floor_divide(y, 2 * np.pi), dtype=np.int)  # .astype(np.int) is slightly slower
        y = y % (2 * np.pi)  # phase vector moduli 2pi
        t = ts[-1]

        if save_counter_r == save_every_r:
            save_counter_r = 0
            r = order_parameter(y)
            rs_coarse.append(r)
        if save_counter_phi == save_every_phi:
            save_counter_phi = 0
            ys_coarse.append(y)
            dws_coarse.append(dw)
            dw = np.zeros_like(y, dtype=np.int)

        save_counter_r += 1
        save_counter_phi += 1

    return ys_coarse, dws_coarse, rs_coarse


def get_rs_filename(irun, ipart, path):
    if ipart == 0:
        return os.path.join(path, f'rs_{irun}.npy')
    else:
        return os.path.join(path, f'rs_{irun}_pt{ipart}.npy')


def get_traj_filename(irun, ipart, path):
    if ipart == 0:
        return os.path.join(path, f'phi_{irun}.npy')
    else:
        return os.path.join(path, f'phi_{irun}_pt{ipart}.npy')


def get_winding_number_increments_filename(irun, ipart, path):
    if ipart == 0:
        return os.path.join(path, f'winding_{irun}.npy')
    else:
        return os.path.join(path, f'winding_{irun}_pt{ipart}.npy')


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
    return np.concatenate(phis_list)  # trajectory glued back from parts


## Load arguments
irun, ncycle_total, D, dt, save_every_r, save_every_phi, sim_name = int(sys.argv[1]), int(sys.argv[2]), float(
    sys.argv[3]), \
                                                                    float(sys.argv[4]), int(sys.argv[5]), int(
    sys.argv[6]), str(sys.argv[7])

# Folder names
objfolder = f'obj/{sim_name}/'
outfolder = f'out/{sim_name}/'

## Logging
carpet.setup_logging(objfolder + f'worker_irun={irun}.log')

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
    ncycle_old = (len(phis_old) - 1) * save_every_phi  # assume that input save_every is the same as used in prev. sims!
    phi0 = phis_old[-1]
    ipart = ipart_last + 1
    del phis_old  # free up memory
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
if ncycle_extra > max(save_every_phi, save_every_r): # would need to handle separately the case of an empty output array
    phis, dws, rs = integrate_cycles(phi0, D, dt, period, ncycle_extra, eps=10 ** -3 * dt, save_every_r=save_every_r,
                                     save_every_phi=save_every_phi)

    if ipart > 0:  # remove the first point because it's the same as the last point in the first part
        phis = phis[1:]
        rs = rs[1:]
    ## Save output
    os.makedirs(outfolder, exist_ok=True)
    # filename = get_traj_filename(irun, ipart, outfolder)
    # np.save(filename, phis)
    filename = get_rs_filename(irun, ipart, outfolder)
    np.save(filename, rs)
    filename = get_traj_filename(irun, ipart, outfolder)
    np.save(filename, phis)

    filename = get_winding_number_increments_filename(irun, ipart, outfolder)
    # dws = np.array(dws, dtype=np.int32) # - > already int32 by default
    np.save(filename, dws)

# logging.info("Finished run {} at random I.Co.; D={:.3E}; dt={:.3E}".format(irun, D, dt))
