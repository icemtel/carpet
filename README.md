## Installation
- Run in terminal `python setup.py develop`. 
- After that one can use the package everywhere.
- Any code changes to the code will be immediately active.

## Structure
Importing **`carpet`** will give direct access to functions inside several modules
- `visualize` (name speaks for itself)
- `dynamics` 
  - Contains functions to solve ODE
  - Define global phase
 
**Other modules** must be imported separately when needed. This is done as following

`import carpet.triangular_lattice as lattice`


- Geometry-specific files, e.g. `lattice_*.py` 
  - Contains functions to build list of cilia positions and neighbours
  - m-twists
  - Functions to load friction matrix, and compile the right side of ODE.
- `parallel_with_threads.py`  - contains code to compute a function for a list of inputs in parallel manner.

- `classes.py` - implementation of cilia classes


