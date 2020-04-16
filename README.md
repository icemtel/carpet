Anton Solovev 2020

Package to solve ODEs for cilia carpets - hence the name `carpet`. 
Can be applied to lattices and networks of phase oscillators.

## Installation
Run in terminal `python setup.py develop`. 
- After that package is ready to be imported.
- Any changes to the code will be applied immediately.
- Package `numba` is an optional requirement. It speeds up some simulations
  (`carpet.physics.friction_pairwise.py`)

## Structure
- **`carpet`** - importing it immediately gives access to several modules, e.g.
    - `dynamics` 
      - Contains functions to solve ODE
      - Define global phase
    - `various`: 
      - a function to setup logging
      - root mean square
      
**Other modules** must be imported separately when needed. This is done as following
`import carpet.visualize as vis`.

- `visualize` - functions for visualization
- Geometry-specific files, are contained in packages inside `lattice`. 
  Recommended use `import carpet.lattice.triangular as lattice`
  - Each contains functions to build a list of cilia positions and neighbours
  - m-twists - metachronal waves 
- Physics/coupling specific files: in `physics` folder.
  - Recommended use `import carpet.physics.kuramoto as physics` 
  
- `parallel_with_threads.py`  - code to run a function in parallel on a list of inputs.
- `classes.py` - implementation of cilia classes; assume symmetry classes and reduce ODE system.

- Many reusable functions reside outside the main package - they are located in `scripts`

TODO: update requirements to python 3.7, etc

## How to use

- See examples
- General idea: import lattice type as `import carpet.lattice.SOME_TYPE as lattice`,
                import physics/coupling as `import carpet.coupling.SOME_TYPE as coupling`
                => define geometry and physics, define right_side_of_ODE
                Then use `define_solve_cycle` and dynamics can be integrated.