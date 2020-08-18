'''
Master script to run attraction_basinsXXXXXX.py in parallel


---

args: num_threads
Memory? Windows: 60 per thread (12x12 carpet), 30 for main
'''
import os
import logging
import subprocess
import sys
import carpet.parallel_with_threads as pwt # to find fixpoints faster
import carpet

## Setup logging
carpet.setup_logging('master.log')

script_path = '.'
script_name = "integrate_trajectory_random.py"

def run_script(*args):
    '''
    args: ibatch, nrun, ncycle, tol,  fixpoint_radius, termination_eps
    '''
    irun = args[0]
    args_str = [str(arg) for arg in args]
    logging.info("Starting: irun={}".format(irun))
    completedProcess = subprocess.run([sys.executable, script_name, *args_str], cwd=script_path)
    if completedProcess.returncode == 1 :
        raise RuntimeError("Subprocess finished with an error")
    logging.info("Finished: irun={}".format(irun))


num_threads = int(sys.argv[1])
tol = 10 ** -6
nrun_range = (0,400) # numbers of trajectories
ncycle = 3000
save_every = 30 # ncycle should be divisible by save_every

os.makedirs('out/basin/', exist_ok=True)
os.makedirs('obj/basin/', exist_ok=True)

list_of_args = [[irun, ncycle, tol, save_every] for irun in range(*nrun_range)]
pwt.run_parallel(num_threads, run_script, list_of_args)

logging.info("master_script (attraction_basins) finished")