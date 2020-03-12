"""
2020-01-27 Master script to run a parallel jobs; input - number of threads to run


args (all compulsory):
    num_threads
    delta0
    tol
"""
import logging
import subprocess
import sys
import carpet.parallel_with_threads as pwt  # to find fixpoints faster
import carpet

## Script
script_path = '.'
script_name = "compute_evals20200225_next_to_nearest.py" # "compute_evals20200127.py"


def run_script(*args):
    '''
    args: k1, k2, delta0, tol
    '''
    args_str = [str(arg) for arg in args]
    logging.info("Starting: ({},{})".format(*args))
    completedProcess = subprocess.run([sys.executable, script_name, *args_str], cwd=script_path)
    if completedProcess.returncode == 1:
        raise RuntimeError("Subprocess finished with an error")
    logging.info("Finished: ({},{})".format(*args))

## Setup logging
carpet.setup_logging('master.log')

## Setup simulations
num_threads = int(sys.argv[1])  # number of threads utilized
delta0 = float(sys.argv[2])  # number of threads utilized
tol = float(sys.argv[3])  # number of threads utilized

nx, ny = 6,6

list_of_args = [(k1,k2, delta0, tol) for k1 in range(nx) for k2 in range(ny)]


## Run simulations
pwt.run_parallel(num_threads, run_script, list_of_args)
logging.info("master_script (compute_evals) finished")
