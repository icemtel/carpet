'''
A script to run master_scripts for eigenvalues
'''
import subprocess
import sys
import carpet.parallel_with_threads as pwt

nxs = [4, 8, 12, 16, 32, 64] # from prepare_computations
nys = nxs
## Make combinations of nx and ny, but without permutations, e.g. if we have (4,8), then don't add (8,4)
# sim_paths = []
# for inx, nx in enumerate(nxs):
#     for ny in nys[inx:]:
#         sim_path = f'nx={nx}_ny={ny}/'
#         sim_paths.append(sim_path)

# All possible combinations of nx,ny
sim_paths = []
for inx, nx in enumerate(nxs):
    for ny in nys:
        sim_path = f'nx={nx}_ny={ny}/'
        sim_paths.append(sim_path)



num_threads = int(sys.argv[1])  # number of threads utilized
num_threads_subprocess = 1
delta0 = 1e-2
tol    = 1e-8


## Script
script_name = "master_evals.py"

def run_script(script_path, *args):
    '''
    args: k1, k2, delta0, tol
    '''
    args_str = [str(arg) for arg in args]
    completedProcess = subprocess.run([sys.executable, script_name, *args_str], cwd=script_path)
    if completedProcess.returncode == 1:
        raise RuntimeError("Subprocess finished with an error")

list_of_args = [(sim_path, num_threads_subprocess, delta0, tol) for sim_path in sim_paths]

pwt.run_parallel(num_threads_subprocess, run_script, list_of_args)