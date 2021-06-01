

# sim_setup

`sim_geometry.py`
- get positions of oscillators, their connections
- save data -> next time load it instead of computing again

`sim_physics.py`
- info about coupling between oscillators, e.g., cilia carpet model, Kuramoto model.
- functions to integrate dynamical equation


# Numerical simulations

Scripts use the same folder logic:
- Use `obj/` for input and intermediate/non-essential results
- Use `out/` for results (trajectories data)


Finding fixed points - `fixpoint`

Linear stability analysis - `linear`

Basins of attraction - `basin`



# Various

- Load trajectories (1 trajectory= several files) -> https://gist.github.com/icemtel/da66f64c8170f2a67bae4a3879a6d396
`slurm_submit.sh` - to submit SLURM jobs on a high-perfomance computer cluster.