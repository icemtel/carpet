"""
2019-10-08 Master script to run a parallel jobs; input - number of threads to run
"""
import logging
import subprocess
import sys
import carpet.parallel_with_threads as pwt
import carpet

## Script
script_path = '.'
script_name = "find_fixpoint20200203.py"


def run_script(*args):
    '''
    args: k1,k2
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
nx, ny = 6,6

list_of_args = [(k1,k2) for k1 in range(nx) for k2 in range(ny)]


## Run simulations
pwt.run_parallel(num_threads, run_script, list_of_args)
logging.info("master_script (find_fixpoint) finished")
