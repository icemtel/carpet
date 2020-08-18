"""
2019-10-08 Master script to run a parallel jobs
2019-11-25: num_of_threads - argument
2019-11-26: removed pickling functions - were unused
"""
import logging
import os
import sys
import subprocess
import carpet.parallel_with_threads as pwt
import carpet
import numpy as np
from sim_setup import period, N

## Script
script_path = '.'
script_name = "prolong_trajectory_noise.py"  # _ROTATED


def integrate_trajectory_script(*args):
    '''
    args: irun, ncycle, D, dt
    '''
    irun = args[0]
    args_str = [str(arg) for arg in args]
    logging.info("Starting: irun={}".format(irun))
    completedProcess = subprocess.run([sys.executable, script_name, *args_str], cwd=script_path)
    if completedProcess.returncode == 1:
        raise RuntimeError("Subprocess finished with an error")
    logging.info(
        "Finished: irun={} Random I.Co. D={:.3e}; dt={:.3E}".format(irun, args[2], args[3]))  # MAYBE: more info


## Setup logging
carpet.setup_logging('master.log')

## Setup simulations
num_threads = int(sys.argv[1])  # number of threads utilized

dt = 0.01 * period
save_every = 50  # must be an integer divisor of ncycle for best efficiency

Ds = 1e-4 * np.array([1.66, 1.83, 2.0]) # noise strength simulated
nstarts =      [0,     0,    0] # read results from nstart to nrun
nruns =        [400,  400,  400]
nrun_ranges = [(nstart, nrun) for nstart,nrun in zip(nstarts, nruns)]
ncycles = [5000, 8000, 6000] # assume that many cycles in trajectory
save_everys =  [100,   80,   60]


# Test arguments
# Ds = [0]
# nrun_ranges = [(0,1)]
# ncycles = [2]
# save_everys = [1]
assert (len(Ds) == len(nrun_ranges) == len(ncycles))


list_of_args = []
for (D, nrun_range, ncycle, save_every) in zip(Ds, nrun_ranges, ncycles, save_everys):
    list_of_args += [[irun, ncycle, D, dt, save_every] for irun in range(*nrun_range)]



# Write initial conditions and compile a list of arguments to pass to the script
list_of_args = []
for (D, nrun_range, ncycle, save_every) in zip(Ds, nrun_ranges, ncycles, save_everys):
    # Input initial condition
    sim_name = "traj_random_D={:.3E}_dt={:.3E}".format(D, dt)

    list_of_args += [[irun, ncycle, D, dt, save_every, sim_name] for irun in range(*nrun_range)]


## Run simulations
pwt.run_parallel(num_threads, integrate_trajectory_script, list_of_args)
logging.info("master_script (integrate_euler; random) finished")