'''
args: num_threads
'''

import logging
import carpet.parallel_with_threads as pwt  # to find fixpoints faster
import carpet
import subprocess
import sys

## Setup logging
carpet.setup_logging('master.log')

script_path = '.'
script_name = "basin_loc_20200330.py"


def run_script(*args):
    '''
    args: irun, k1, k2, eps, ncycle, tol, save_every
    '''
    irun = args[0]
    k1, k2 = args[1], args[2]
    eps = args[3]
    args_str = [str(arg) for arg in args]
    logging.info("Starting: irun={}, k1={}, k2={}, eps={}".format(irun, k1, k2, eps))
    completedProcess = subprocess.run([sys.executable, script_name, *args_str], cwd=script_path)
    if completedProcess.returncode == 1:
        raise RuntimeError("Subprocess finished with an error")
    logging.info("Finished: irun={}, k1={}, k2={}, eps={}".format(irun, k1, k2, eps))


## Params
num_threads = int(sys.argv[1])
tol = 10 ** -6
k12s = [(1, 3)]
epses = [1 / 3 + 1 / 12, 1 / 3 + 2 / 12, 1 / 3 + 3 / 12]  # [0.01, 0.1, 1/3, 2/3]
nrun = 300  # trajectories per eps
ncycle = 1500  # length of trajectory
save_every = 50

## For testing
# epses = [0.01, 0.1]
# nrun = 2
# ncycle = 1
# save_every = 1

list_of_args = [[irun, *k12, eps, ncycle, tol, save_every] for k12 in k12s for eps in epses for irun in range(nrun)]

## Add also data for another fixed point


pwt.run_parallel(num_threads, run_script, list_of_args)

logging.info("master_script (attraction_basins) finished")
