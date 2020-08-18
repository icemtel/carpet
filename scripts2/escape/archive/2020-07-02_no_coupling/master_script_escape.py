"""
2019-10-08 Master script to run a parallel jobs
2019-11-25: num_of_threads - argument
2019-11-26: removed pickling functions - were unused
"""
import logging
import subprocess
import sys, os
import carpet.parallel_with_threads as pwt
import carpet
import numpy as np
from carpet.various import define_load_object

## Script
script_path = '.'
script_name = "integrate_trajectory20200612.py"  # _ROTATED


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

k1, k2 = 3, 1
period = 31.25
dt = 0.01 * period

Ds = 1e-4 * np.array([0.05, 0.1, 0.5, 1.0, 1.5, 2])
nrun_ranges = [(0, 300), (0, 300), (0, 300), (0, 300), (0, 300), (0, 300)]
ncycles = [3000, 2000, 1000, 500, 200, 100]
save_everys = [50, 20, 10, 5, 2, 1]
# Testing arguments
# Ds = [0]
# nrun_ranges = [(0,1)]
# ncycles = [2]
# assert (len(Ds) == len(nrun_ranges) == len(ncycles))

list_of_args = []
for (D, nrun_range, ncycle, save_every) in zip(Ds, nrun_ranges, ncycles, save_everys):
    # Input initial condition
    output_folder = "traj_k1={}_k2={}_D={:.3E}_dt={:.3E}/".format(k1, k2, D, dt)
    os.makedirs(output_folder, exist_ok=True)
    load_object = define_load_object('.')
    fixpoint_dict = load_object('fixpoint_dict.pkl')
    phi0 = fixpoint_dict[k1, k2]
    np.save(output_folder + '/phi0.npy', phi0)

    list_of_args += [[irun, ncycle, D, dt, save_every, output_folder] for irun in range(*nrun_range)]

## Run simulations
pwt.run_parallel(num_threads, integrate_trajectory_script, list_of_args)
logging.info("master_script (integrate_euler; random) finished")
