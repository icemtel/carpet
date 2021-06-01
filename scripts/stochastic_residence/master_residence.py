"""
"""
import logging
import os
import sys
import subprocess
import carpet.parallel_with_threads as pwt
import carpet
import numpy as np
from sim_physics import period, N

## Script
script_path = '.'
script_name = "intergrate_trajectory_noise_run_or_prolong_until.py"  # _ROTATED


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

dt = 0.01 * period
D0 = 1e-4
Ds = D0 * np.array([0.5, 1.0, 1.25, 1.5, 1.7, 1.8, 1.9, 2.2, 2.5])  # noise strength simulated
nrun_ranges = [(0, 400), (0, 300), (0, 200), (0, 200), (0, 200), (0, 200), (0, 200), (0, 200), (0, 200)]
ncycles_real = [8000, 14992, 20000, 40000, 70000, 100000, 60000, 14992, 8000]  # assume that many cycles in trajectory
save_everys = [16, 16, 16, 16, 16, 16, 16, 16, 8]

isets = list(range(len(Ds)))

nstarts = [nrun_range[0] for nrun_range in nrun_ranges]
nruns = [nrun_range[1] for nrun_range in nrun_ranges]

ncycles = [ncycle // save_every for ncycle, save_every in zip(ncycles_real, save_everys)]

assert len(Ds) == len(nrun_ranges) == len(ncycles) == len(save_everys)
# Test arguments
# Ds = [0]
# nrun_ranges = [(0,1)]
# ncycles = [7]
# save_everys = [1]


list_of_args = []
for (D, nrun_range, ncycle, save_every) in zip(Ds, nrun_ranges, ncycles_real, save_everys):
    # Input initial condition
    sim_name = "residence/traj_random_D={:.3E}_dt={:.3E}/".format(D, dt)

    # Input initial condition
    sim_name = "residence/traj_random_D={:.3E}_dt={:.3E}/".format(D, dt)
    objfolder = f'obj/{sim_name}/'
    outfolder = f'out/{sim_name}/'
    os.makedirs(objfolder, exist_ok=True)
    # Save input - if does not exist already
    for irun in range(*nrun_range):
        filename = objfolder + f'/phi0_{irun}.npy'
        if not os.path.isfile(filename):
            phi0 = carpet.get_phase_random(N)
            np.save(filename, phi0)

    list_of_args += [[irun, ncycle, D, dt, save_every, sim_name] for irun in range(*nrun_range)]

## Run simulations
pwt.run_parallel(num_threads, integrate_trajectory_script, list_of_args)
logging.info("master_script (integrate_euler; random) finished")
