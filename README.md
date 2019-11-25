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
- `various`: 
  - a function to setup logging
  - root mean square
 
**Other modules** must be imported separately when needed. This is done as following

`import carpet.lattice_triangular as lattice`


- Geometry-specific files, e.g. `lattice_triangular.py`, `lattice_1D_chain.py` 
  - Each contains functions to build a list of cilia positions and neighbours
  - m-twists
  - Functions to load friction matrix, and compile the right side of ODE.
- `parallel_with_threads.py`  - contains code to compute a function for a list of inputs in parallel manner.

- `classes.py` - implementation of cilia classes


TODO: update requirements to python 3.7, etc