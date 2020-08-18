'''
Integrate with Euler method with noise

- Use objfolder for input and intermediate/non-essential results
- Use outfolder for results (trajectories data)

INPUT: irun, ncycle, D, dt, save_every, sim_name

2020-02-19: upd scipy->numpy; upd carpet
2020-03-23: new argument: save_every -> save every n-th cycle only
2020-07-07: unified framework for residence and escape analysis
'''

# import packages needed below


from sim_setup import carpet, right_side_of_ODE, period, N
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
objfolder = f'obj/{sim_name}/'
outfolder = f'out/{sim_name}/'

# Load input
input_filename = objfolder + f'phi0_{irun}.npy'
phi0 = np.load(input_filename)

## Save rng state
state = np.random.get_state()
with open(objfolder + 'state_{}.pkl'.format(irun), 'wb') as f:
    pickle.dump(state, f, pickle.HIGHEST_PROTOCOL)

## Run simulation
phis, ts = integrate_cycles(phi0, D, dt, period, ncycle, eps=10 ** -3 * dt, save_every=save_every)

## Save output
os.makedirs(outfolder, exist_ok=True)
output_filename = outfolder + f'phi_{irun}.npy'
np.save(output_filename, phis)

# logging.info("Finished run {} at random I.Co.; D={:.3E}; dt={:.3E}".format(irun, D, dt))