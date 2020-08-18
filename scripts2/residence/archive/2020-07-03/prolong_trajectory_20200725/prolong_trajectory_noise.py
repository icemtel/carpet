'''
Integrate with Euler method with noise

- Use objfolder for input and intermediate/non-essential results
- Use outfolder for results (trajectories data)

INPUT: irun, ncycle, D, dt, save_every, sim_name

2020-02-19: upd scipy->numpy; upd carpet
2020-03-23: new argument: save_every -> save every n-th cycle only
2020-07-07: unified framework for residence and escape analysis
- edited from integrate_trajectory_noise
- goal: prolong old trajectories
'''

# import packages needed below
import carpet
from sim_setup import right_side_of_ODE, period, N
import  sys, os
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
        y = ys[-1] % (2 * np.pi)
        t = ts[-1]

        if save_counter == save_every:
            save_counter = 1
            ys_coarse.append(y)
            ts_coarse.append(t)
        else:
            save_counter += 1

    return ys_coarse, ts_coarse


## Load arguments
irun, ncycle, D, dt, save_every, sim_name = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), \
                                            float(sys.argv[4]), int(sys.argv[5]), str(sys.argv[6])
# Folder names
old_folder = f'{sim_name}/'
new_folder = f'{sim_name}_long/'
os.makedirs(new_folder, exist_ok=True)

# Load input
load_object = carpet.various.define_load_object(old_folder)
input_filename = f'phi_{irun}.pkl' # old trajectory
phis = load_object(input_filename)
phi0 = phis[-1]

## Save rng state
dump_object = carpet.various.define_dump_object(new_folder)
state = np.random.get_state()
dump_object(state, f'state_{irun}.pkl')
## Run simulation
phis_extra, ts = integrate_cycles(phi0, D, dt, period, ncycle, eps=10 ** -3 * dt, save_every=save_every)


phis_new = np.concatenate([phis, phis_extra[1:]]) # remove the first point because it's the same as the last point in the first part
## Save output

output_filename = f'phi_{irun}.pkl'
dump_object(phis_new, output_filename)

# logging.info("Finished run {} at random I.Co.; D={:.3E}; dt={:.3E}".format(irun, D, dt))