machemer_4 vs machemer_3 and machemer_3M:

2020-11-04: fixed a bug - g22 was constant; now it's the same as g11 (up to order of arguments)
2021-01-25: add more relative positions (for rotated lattice)

---
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