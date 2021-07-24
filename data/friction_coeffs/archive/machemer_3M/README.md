machemer_3M vs machemer_3:

- found out that the beat pattern was mirrored =>
   to account for this, in this data set translations are mirrored around y-axis
   
   
---
machemer_3:

Produced with 2020-10-01_cilia_manuscript_code

- Computed with background cilia
- Fitted Fourier coefficients to data points with different background cilia phase
- IF hydrodynamic computation failed to converge to 1e-3 (normal tol=5e-4), then points were excluded from the fit target.

---
- Use machemer beat pattern
- finer cilia mesh (60x8 (+2) nodes)
- fixed FBEM convergence; tolerance=5e-4
- cilia radius -> 0.125 [um], cilia elevation ->0.25[um]


- Computed only g21 from two cilia experiments;
- g12 - a copy of g21
- g11 is taken from single cilium experiment
- g22 is taken as transpose of g11


UPD: was a bug in g22; use machemer_4 instead