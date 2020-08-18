"""
2020-01-27 Master script to run a parallel jobs; input - number of threads to run

TODO: implement with .npy

args:
    num_threads
    delta0
    tol
Implicit:
    nx,ny
"""
import logging
import subprocess
import sys
import carpet.parallel_with_threads as pwt  # to find fixpoints faster
import carpet

## Script
script_path = '.'
script_name = "compute_evals20200420.py" # CHECK nx, ny

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

nx, ny = 12,12

list_of_args = [(k1,k2, delta0, tol) for k1 in range(nx) for k2 in range(ny)]


## Load calculated matrices
## Calculate eigenvalues and logarithms
import numpy as np
import scipy.linalg as lin # or use numpy.linalg https://docs.scipy.org/doc/scipy/reference/tutorial/linalg.html
from carpet.various import define_dump_object, define_load_object

dump_object = define_dump_object()
load_object = define_load_object()

def calc_sort_log_evals_evecs(L):
    '''
    :return: logairthm of eigenvalues, eigenvectors
    '''
    # Ljapunov spectrum
    evals, evecs = lin.eig(L)
    evals = np.log(evals)
    # Sort eigenvalues
    idx = evals.argsort()
    evals = evals[idx]
    evecs = evecs[:, idx]
    return evals, evecs


evals_dict = {} # Dictionary with eigenvalues, corresponding to each of m-twists
evecs_dict = {}

output_path =  'L_delta0={:.3E}_tol={:.3E}/'.format(delta0, tol)

for k1,k2, _, _ in list_of_args:
    filename = "L_mtwist_k1={}_k2={}.pkl".format(k1,k2)
    L = load_object(filename, output_path)
    evals, evecs = calc_sort_log_evals_evecs(L)
    evals_dict[k1,k2] = evals
    evecs_dict[k1,k2] = evecs

# Save eigenvalues and eigenvectors

filename = 'evals_dict_nx={}_ny={}_delta0={:.3E}_tol={:.3E}.pkl'.format(nx,ny, delta0, tol)
dump_object(evals_dict, filename)

filename = 'evecs_dict_nx={}_ny={}_delta0={:.3E}_tol={:.3E}.pkl'.format(nx,ny, delta0, tol)
dump_object(evecs_dict, filename)