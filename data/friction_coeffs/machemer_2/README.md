Produced with 2020-10-01_LAMAS_SPP_code

- Use machemer beat pattern
- finer cilia mesh (60x8 (+2) nodes)
- fixed FBEM convergence; tolerance=5e-4
- cilia radius -> 0.125 [um], cilia elevation ->0.25[um]


- Computed only g21 from two cilia experiments;
- g12 - a copy of g21
- g11 is taken from single cilium experiment
- g22 is taken as transpose of g11

TODO: check is it mirrored; compared to machemer_4 
[in that case, machemer_4 is correct]