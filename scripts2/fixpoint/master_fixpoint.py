"""
2019-10-08 Master script to run a parallel jobs; input - number of threads to run

Check nx and ny!

args: number of threads
Implicit:
    nx,ny
"""
import logging
import subprocess
import sys
import numpy as np
import carpet.parallel_with_threads as pwt
import carpet
from sim_setup import nx, ny


## Script
script_path = '.'
script_name = "find_fixpoint20200420.py"


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

## Run simulations
list_of_args = [(k1,k2) for k1 in range(nx) for k2 in range(ny)]
pwt.run_parallel(num_threads, run_script, list_of_args)
logging.info("master_script (find_fixpoint) finished")


## Combine results
from carpet.various import define_dump_object

dump_object = define_dump_object()
# Merge
# ## Combine results
fixpoint_dict = {} # Dictionary with eigenvalues, corresponding to each of m-twists

for k1,k2 in [(k1,k2) for k1 in range(nx) for k2 in range(ny)]:
    filename = "out/fixpoint/fixpoint_k1={}_k2={}.npy".format(k1,k2)
    fixpoint = np.load(filename)
    fixpoint_dict[k1,k2] = fixpoint.copy()

filename = 'fixpoint_dict.pkl'
dump_object(fixpoint_dict, filename)