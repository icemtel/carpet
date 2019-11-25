"""
2019-10-08
Master script to parallelize trajectories integration.
"""
import logging
import time
import os
import pickle
import carpet.parallel_with_threads as pwt
import carpet

## Setup logging
carpet.setup_logging('master.log')

## Shortcuts to save objects
objfolder = 'obj/'
outfolder = 'out/'
os.makedirs(objfolder,exist_ok=True)


def dump_object(obj, filename, path=objfolder):
    filename = os.path.join(path, filename)
    print(filename)
    with open(filename, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_object(filename, path=objfolder):
    filename = os.path.join(path, filename)
    with open(filename, 'rb') as f:
        obj = pickle.load(f)
    return obj


import subprocess
import sys

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
    logging.info("Finished: irun={} ({},{})-twist D={:.3e}".format(irun,args[2],args[3],args[4])) # TODO: more info

period = 31.25

dt = 0.01 * period
ncycle = 4000
nrun_start = 30
nrun_end = 100 # trajectories to simulate

list_of_args = []
# Different Ds
k1,k2 = 2,1
Ds = [10 ** -4]
for D in Ds:
    list_of_args += [[irun, ncycle, k1,k2, D, dt] for irun in range(nrun_start, nrun_end)]


pwt.run_parallel(4, integrate_trajectory_script, list_of_args)

logging.info("master_script (integrate_euler) finished")