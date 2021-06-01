
Subfolders contain various simulation scripts

Each simulation should contain:
-  Setup: `sim_geometry.py`  and `sim_physics.py` --- contains a setup for geometry and physics.
    A collection of such files is located in `sim_setup` folder
- `worker_foo.py` --- a script to run a single simulation/computation.
Imports sim_geometry.py`  and `sim_physics.py`, such that the script is reusable with different inputs.
- `master_foo.py` 
 calls `worker_foo.py` many times, each time with different parameters, uses multi-threading.

Example: `fixpoint` - the script tries to find a fixed point in vicinity of a plane wave. 
           The master script is iterating through different wave modes.

## sim_setup

`sim_geometry.py`
- get positions of oscillators, their connections
- save data -> next time load it instead of computing again

`sim_physics.py`
- info about coupling between oscillators, e.g., cilia carpet model, Kuramoto model.
- functions to integrate dynamical equation

## Numerical simulations

Scripts use the same folder logic:
- Use `obj/` for input and intermediate/non-essential results
- Use `out/` for results (trajectories data)
     
Finding fixed points - `fixpoint`
- the script tries to find a fixed point in vicinity of a plane wave. 
  The master script is iterating through different wave modes.
  

Linear stability analysis - `linear`

Basins of attraction - `basin`


# Various

- Load trajectories (1 trajectory= several files) -> https://gist.github.com/icemtel/da66f64c8170f2a67bae4a3879a6d396
`slurm_submit.sh` - to submit SLURM jobs on a high-perfomance computer cluster.