Subfolders contain various simulation scripts
- Each simulation contains of at least three files:
  - `sim_setup.py`  --- contains a setup for geometry and physics.
    A collection of such files is located in `sim_setup` folder
  - `foo.py` --- a script to run a single simulation/computation.
    Imports `sim_setup.py`, such that the script is reusable with different inputs.
  - `master_foo.py` 
     Calls `foo.py` with different parameters, uses multi-threading.
     
Example: `fixpoint`, the script tries to find a fixed point in vicinity of a plane wave. 
          The master script is iterating through different wave modes.
