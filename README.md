# Solve dynamics of N>>1 coupled oscillators

Python package to study systems of coupled phase oscillators:
- 2D Kuramoto model
- cilia carpet model (references [1, 2])  - hence the name `carpet`.
- new models can be easily added

Can be applied to lattices and networks of phase oscillators.


## How to use

- See examples
- General idea: import lattice type as `import carpet.lattice.SOME_TYPE as lattice`,
                import physics/coupling as `import carpet.coupling.SOME_TYPE as coupling`
                => define geometry and physics, define right_side_of_ODE
                Then use `define_solve_cycle` and dynamics can be integrated.

## Installation
Run in terminal `python setup.py develop`. 
- After that package is ready to be imported.
- Any changes to the code will be applied immediately.
- Package `numba` is optional, but highly recommended. It speeds up some computations
  (`carpet.physics.friction_pairwise`, `carpet.physics.kuramoto_numpy`).

## Structure
- `carpet` - import it to get access to several modules, e.g.
    - `dynamics` 
      - Functions to solve ODE
    - `various`: 
      - logging for simulations
      - math functions
-  Geometry-specific files (topology: chain, lattice, network?) 
   are contained in packages inside `lattice` folder. 
   Files contain functions to build a list of positions and connetions depending on the geometry of choice.
- Physics-specific files (how are oscillators coupled?): in `physics` folder.
  - Recommended use `import carpet.physics.kuramoto as physics` 
  
Other modules should be imported separately, e.g., `import carpet.visualize as vis`.

- `visualize` - functions for visualization

- `parallel_with_threads.py`  - code to run a function in parallel on a list of inputs.
- `classes.py` - implementation of cilia classes; assume symmetry classes and reduce dimensionality of the ODE system.
- Many reusable functions reside outside the main package - they are located in `scripts`

### Authors

- [Anton Solovev](https://github.com/icemtel)
- [Benjamin M. Friedrich](https://cfaed.tu-dresden.de/friedrich-home) benjamin.m.friedrich@tu-dresden.de

Publication to cite: [2]

- [1]: [Solovev & Friedrich 2020 EPJ E ST](https://link.springer.com/article/10.1140/epje/s10189-021-00016-x);  also available on [arXiv](https://arxiv.org/abs/2010.08111 ) 
- [2]: [Solovev & Friedrich 2020b arXiv](https://arxiv.org/abs/2012.11741)
