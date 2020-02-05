'''
Master script to run attraction_basinsXXXXXX.py in parallel

---- (probably for 12x12 carpet), 10 ** -6
- Renate: 8s per cycle -> 2.22 hrs for 1000 cycles
  => guess 3h per 1000 cycles at Taurus
  => 133 traj doable by 16 cores in 25 h
  Be safe: run 100 traj
----
6x6: try 200 in 24h on Renate
----

args: num_threads
Memory? Windows: 60 per thread (12x12 carpet), 30 for main
'''

import logging
import os
import carpet.parallel_with_threads as pwt # to find fixpoints faster
import carpet
import subprocess
import sys

## Setup logging
carpet.setup_logging('master.log')

script_path = '.'
script_name = "attraction_basins_20200128.py"

def run_script(*args):
    '''
    args: ibatch, nrun, ncycle, tol,  fixpoint_radius, termination_eps
    '''
    ibatch = args[0]
    args_str = [str(arg) for arg in args]
    logging.info("Starting: ibatch={}".format(ibatch))
    completedProcess = subprocess.run([sys.executable, script_name, *args_str], cwd=script_path)
    if completedProcess.returncode == 1 :
        raise RuntimeError("Subprocess finished with an error")
    logging.info("Finished: ibatch={}".format(ibatch))

num_threads = int(sys.argv[1])
tol = 10 ** -6
nbatch_range = (0, 20)
nrun = 10 # trajectories per batch
ncycle = 1000
fixpoint_radius = 0.2
termination_eps = 2 * 10 ** - 3 # terminate cycles if we are close enough to a fixpoint


list_of_args = [[ibatch, nrun, ncycle, tol,  fixpoint_radius, termination_eps] for ibatch in range(*nbatch_range)]

pwt.run_parallel(num_threads, run_script, list_of_args)

logging.info("master_script (attraction_basins) finished")