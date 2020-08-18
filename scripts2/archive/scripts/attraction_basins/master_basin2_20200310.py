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
script_name = "attraction_basins2_20200310_next_to_nearest.py"

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
nrun = 2000 # trajectories per batch
ncycle = 1500


list_of_args = [[irun, ncycle, tol] for irun in range(nrun)]

pwt.run_parallel(num_threads, run_script, list_of_args)

logging.info("master_script (attraction_basins) finished")