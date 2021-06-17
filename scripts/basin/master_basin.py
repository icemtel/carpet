'''
args: num_threads, ncycle
'''
import os
import logging
import subprocess
import sys
import carpet.various.parallel_with_threads as pwt  # to find fixpoints faster
import carpet

## Setup logging
carpet.setup_logging('master.log')

script_path = '.'
script_name = "worker_basin.py"
sim_name = 'basin'


def run_script(*args):
    '''
    args: ibatch, nrun, ncycle, tol,  fixpoint_radius, termination_eps
    '''
    irun = args[0]
    args_str = [str(arg) for arg in args]
    logging.info("Starting: irun={}".format(irun))
    completedProcess = subprocess.run([sys.executable, script_name, *args_str], cwd=script_path)
    if completedProcess.returncode == 1:
        raise RuntimeError("Subprocess finished with an error")
    logging.info("Finished: irun={}".format(irun))


num_threads = int(sys.argv[1])
tol = 10 ** -6
nrun_range = (0, 500)  # numbers of trajectories
ncycle = 1020
save_every = 30  # ncycle should be divisible by save_every

os.makedirs(f'out/{sim_name}', exist_ok=True)
os.makedirs(f'obj/{sim_name}', exist_ok=True)

# Generate initial conditions
# -> in this version they should be created in advance
#
# completedProcess = subprocess.run([sys.executable, 'get_initial_conditions_random.py',
#                                    str(nrun_range[0]), str(nrun_range[1])], cwd=script_path)
# if completedProcess.returncode == 1:
#     raise RuntimeError("Subprocess finished with an error")

# Integrate trajectories
list_of_args = [[irun, ncycle, tol, save_every, sim_name] for irun in range(*nrun_range)]
pwt.run_parallel(num_threads, run_script, list_of_args)

logging.info("master_script (attraction_basins) finished")
