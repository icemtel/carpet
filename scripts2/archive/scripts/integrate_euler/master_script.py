"""
2019-10-08
Master script to parallelize trajectories integration.
2019-11-25: num_threads as input argument
2019-11-27: remove pickling functions - were unused
"""
import logging
import subprocess
import sys
import carpet.parallel_with_threads as pwt
import carpet

## Script
script_path = '.'
script_name = "integrate_trajectory20190926.py"


def integrate_trajectory_script(*args):
    '''
    args: irun, ncycle, k1,k2, D, dt
    '''
    irun = args[0]
    args_str = [str(arg) for arg in args]
    logging.info("Starting: irun={}".format(irun))
    completedProcess = subprocess.run([sys.executable, script_name, *args_str], cwd=script_path)
    if completedProcess.returncode == 1:
        raise RuntimeError("Subprocess finished with an error")
    logging.info("Finished: irun={} ({},{})-twist D={:.3e}".format(irun, args[2], args[3], args[4]))  # MAYBE: more info


## Setup logging
carpet.setup_logging('master.log')

## Setup simulations
num_threads = int(sys.argv[1])
period = 31.25
dt = 0.01 * period
ncycle = 4000
nrun_range = (30, 100)  # trajectories to simulate

list_of_args = []
k1, k2 = 2, 1
Ds = [10 ** -4]  # Different Ds
for D in Ds:
    list_of_args += [[irun, ncycle, k1, k2, D, dt] for irun in range(*nrun_range)]

## Run simulations
pwt.run_parallel(num_threads, integrate_trajectory_script, list_of_args)
logging.info("master_script (integrate_euler) finished")
