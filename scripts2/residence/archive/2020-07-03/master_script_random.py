"""
2019-10-08 Master script to run a parallel jobs
2019-11-25: num_of_threads - argument
2019-11-26: removed pickling functions - were unused
"""
import logging
import subprocess
import sys
import carpet.parallel_with_threads as pwt
import carpet
import numpy as np

## Script
script_path = '.'
script_name = "integrate_trajectory_random20200323.py"  # _ROTATED


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

period = 31.25
dt = 0.01 * period
save_every = 50  # must be an integer divisor of ncycle for best efficiency

# ncycle = 6000
# nrun_range = (0, 300)  # trajectories to simulate
#
# list_of_args = []
# # Different Ds
# Ds = [1e-5, 1e-4]
# for D in Ds:
#     list_of_args += [[irun, ncycle, D, dt, save_every] for irun in range(*nrun_range)]

Ds =  1e-4 *    np.array([0.1]) #, 1 / 3, 1])
nrun_ranges = [(0, 400), (0, 200), (0,100)]
ncycles = [4500, 8000, 6000] #

# Test arguments
Ds = [0]
nrun_ranges = [(0,1)]
ncycles = [2]
assert (len(Ds) == len(nrun_ranges) == len(ncycles))

list_of_args = []
for (D, nrun_range, ncycle) in zip(Ds, nrun_ranges, ncycles):
    list_of_args += [[irun, ncycle, D, dt, save_every] for irun in range(*nrun_range)]

## Run simulations
pwt.run_parallel(num_threads, integrate_trajectory_script, list_of_args)
logging.info("master_script (integrate_euler; random) finished")
