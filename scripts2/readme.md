Version 2?

- use a combination of setup and simulations scripts
- Setup: contains everything about geometry and phsyics
- Simulation: contains the main job of the script, e.g. calculate eigenvalues, or fixed points.

Pros:
- This way I don't need to change every script if I change the data set - only the setup.
- Scripts are shorter
- Keep benchmarks in "setup" file
Cons: 
- syntax highlighting is not working properly -> Solution: import explicitly the objects needed