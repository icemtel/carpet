'''
- Create scripts with simulation parameters for various number of oscillators
- Copy master_evals and compute_evals scripts
'''
import numpy as np
import os
from shutil import copyfile


nxs = [4, 8, 12, 16, 32, 64]
nys = nxs
a = 18  # [um] lattice spacing


##
## Parameters
# Physics
period = 31.25            # [ms] period of cilia beat
freq = 2 * np.pi / period # [rad/ms] angular frequency
sin_str = 0.0016 * freq

# # Make combinations of nx and ny, but without permutations, e.g. if we have (4,8), then don't add (8,4)
# nxys = []
# for inx, nx in enumerate(nxs):
#     for ny in nys[inx:]:
#         nxys.append((nx,ny))
## Prepare all different combinations of nx and ny
# Make combinations of nx and ny, but without permutations, e.g. if we have (4,8), then don't add (8,4)
nxys = []
for inx, nx in enumerate(nxs):
    for ny in nys:
        nxys.append((nx,ny))



for nx, ny in nxys:
    # Write a sim_setup
    script_code = f"""import numpy as np
import carpet
import carpet.lattice.triangular as lattice

## Parameters
# Physics
period = {period}        # [ms] period of cilia beat
freq = {freq} # [rad/ms] angular frequency

# Geometry
nx = {nx}
ny = {ny}  # even number
N = nx * ny
a = {a} # [um] lattice spacing

## Initialize
# Geometry

L1,L2 = lattice.get_domain_sizes(nx,ny ,a)
coords, lattice_ids = lattice.get_nodes_and_ids(nx, ny, a)  # get cilia (nodes) coordinates

distances = [1]
NN, TT = lattice.get_neighbours_list(coords, nx, ny, a, distances)
e1, e2 = lattice.get_basis()
get_k = lattice.define_get_k_fbz(nx, ny, a)
get_mtwist = lattice.define_get_mtwist(coords, nx, ny, a)

# Physics
import carpet.physics.kuramoto as physics

sin_str = {sin_str}
coupling = physics.define_sine_coupling(sin_str)
right_side_of_ODE = physics.define_right_side_of_ODE(coupling, freq, NN, TT)
solve_cycle = carpet.define_solve_cycle(right_side_of_ODE, 2 * period, carpet.get_mean_phase)
"""

    sim_path = f'nx={nx}_ny={ny}/'
    os.makedirs(sim_path, exist_ok=True)
    with open(sim_path + 'sim_setup.py', 'w') as f:
        f.write(script_code)


    # Copy master_evals and compute_evals scripts
    evals_script_name = 'compute_evals_no_fp_dict.py'
    master_script_name = 'master_evals.py'

    copyfile(evals_script_name, sim_path + evals_script_name)
    copyfile(master_script_name, sim_path + master_script_name)
