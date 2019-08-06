---
html:
  embed_local_images: true
  embed_svg: true
  offline: true
  toc: true

print_background: false
export_on_save:
  html: true
---

## Notation
- $\mathbf{k}$ - wave vector
- $\mathbf{\Phi}= (\varphi_1, \varphi_2, .., \varphi_{N})^T $ - phase vector
- $ \mathbf{Q(\Phi)}$ - vector of active driving forces
- $ \mathbf{\Gamma}$ - friction coefficients NxN-matrix
- $\mathbf{\Phi_k}$ - m-twist solution
- $\mathcal{L}: H \rightarrow H$ - Poincare map
- $\mathbf{\Phi^*} \in H$ - fix point of the Poincare map
- $\mathbf{\Phi^{*}_k}$ - fix point close to the m-twist solution
- $\mathbf{\Delta_0}$, $\mathbf{\Delta_1}$ - perturbation initial and after one cycle: $\mathcal{L}(\mathbf{\Phi^{*}} + \mathbf{\Delta_0}) = \mathbf{\Phi^{*}} + \mathbf{\Delta_1}$ 
- $\mathbf{L} = \mathrm{D}\mathcal{L}(\mathbf{\Phi^*})$ - linearized Poincare map at a fixed point (matrix). In code `Lmat`
- $\mathbf{L} = e^\mathbf{\Lambda}; \quad \mathbf{\Lambda} = \log \mathbf{L} $ - logarithm of the linearized Poincare map. In code `Lmat_log`
- $\lambda_j$ - eigenvalues of $\Lambda$



## Geometry

*Unit vectors*
$$
\mathbf{e}_1 =\left( \begin{array}{c} 1 \\ 0 \end{array} \right), \quad
\mathbf{e}_2 =\left( \begin{array}{c} \frac{1}{2} \\  \frac{\sqrt{3}}{2} \end{array} \right).
$$

*Position vectors*
$$
\mathbf{x}_{n,m} = n\,a\,\mathbf{e}_1 +  m\,a\,\mathbf{e}_2
$$

*Honeycomb lattice*
$$
\mathcal{H} = \{ \mathbf{x}_{n,m} | 2n+m \equiv 0 \text{ or } 2 \text{ mod } 3\}
$$

*Triangular lattice*
$$
\mathcal{T} = \{ \mathbf{x}_{n,m} | n,m \in \mathrm{Z} \}
$$


## m-twist solutions

*Wave vector of metachronal wave consistent with periodic boundary condition*
$$ \mathbf{k} = \frac{a_1}{L_1}\,k_1\,\mathbf{a}_1^* + \frac{a_2}{L_2}\,k_2\,\mathbf{a}_2^*, $$
where $k_1,k_2\in{Z}$ and $a_1=a$, $a_2 = \sqrt{3}a/2$

*Meta-chronal wave solutions:*

Phase-space vector
$ \Phi (t) = (\varphi_1, \varphi_2, .., \varphi_{N})^T $
with components
$$
\varphi_i(t) = \varphi(t) - \mathbf{k} \cdot \mathbf{x}_i,
$$
for some global phase $\varphi(t)$.
We expect $\dot{\varphi}\approx\omega_0$, but only approximately, since the calibration of active driving forces does not take into account the (weak) hydrodynamic interactions.

### Dynamic equation from force-balance equation

$$ \dot{\Phi} = \mathbf{\Gamma}^{-1}(\Phi)Q(\Phi)  $$
This linear system can be solved directly in Python.


### Calibration from hydrodynamic simulations

*Pair-wise coupling functions:*
$$ \Gamma_{ij}(\Phi) = 0 \quad \text{if not neighbours} $$
$$ \Gamma_{ij}(\Phi) = \Gamma_{12}^{\mathrm{(loc)}}(\varphi_i,\varphi_j; d, \psi) \quad \text{if neighbours} $$
where the translation vector pointing from cilium $i$ to cilium $j$ reads
$$[n(j)-n(i)]\,\mathbf{e}_1 + [m(j)-m(i)]\,\mathbf{e}_2 = d \left( \begin{array}{c} \sin\psi \\ \cos\psi \end{array} \right).$$
*Self-friction:*

$$ \Gamma_{ii}(\Phi) = \Gamma_{11}^{\mathrm{(loc)}}(\varphi_i, \varphi_j; d, \psi) $$

Note, self-friction $\Gamma_{11}^{\mathrm{(loc)}}(\varphi_i, \varphi_j; d, \psi) $ dependence on $\varphi_j$, $d$, $\psi$ is so weak, that we can take
$\Gamma_{11}^{\mathrm{(loc)}}(\varphi_i, \varphi_j; d, \psi)= \Gamma_{11}^{\mathrm{(loc)}}(\varphi_i)$



*Active driving forces (simplest case: calibration for isolated cilium):*
$$ Q_{i}(\Phi) = \Gamma_{11}^{\mathrm{(loc)}}(\varphi_i) \omega_0 $$


### Global Phase

#### Naive phase
For some index $j$, $\varphi = \mathbf{\varphi}_j$

#### Circular average phase
We define order parameters
[[Kuramoto, 1984]](https://link.springer.com/chapter/10.1007/978-3-642-69689-3_7)  [[Strogatz, 2000]](https://www.sciencedirect.com/science/article/pii/S0167278900000944)
for each metachronal wave vector $\mathbf{k}$

$$S_\mathbf{k} = \left| \frac{1}{N}\sum_j \exp\left[ i ( \varphi_j - \varphi_{\mathbf{k},j} ) \right] \right|,$$

where $\varphi_{\mathbf{k},j}$ - j-th component of $\Phi_\mathbf{k}$, m-twist solution defined by vector $\mathbf{k}$.

**Global phase:**
For a fixed $\mathbf{k}$
$$\varphi(t; \mathbf{k}) = \arg \left( \frac{1}{N}\sum_j \exp\left[ i ( \varphi_j(t) - \varphi_{\mathbf{k},j} ) \right] \right) $$

- Only has meaning when $S_\mathbf{k}$ is close to 1.
#### Mean phase
Define the global phase as algebraic mean of all phases
$$
\varphi =  \frac{1}{N}\sum_{j=1}^{N} \varphi_j
$$

Pros
- Additive: if $\mathbf{\Phi} = \mathbf{\Phi_0}+ \mathbf{\Phi_1} $, then $\varphi = \varphi_0 + \varphi_1$
- As an implication, it is independent of which initial condition we consider.

Cons
-$\varphi_j$ can't jump from  $2 \pi$ to $0$, otherwise the global phase will make a jump by $- \frac{2 \pi}{N}$.

## Poicnare map for a m-twist solution
### Poincare map and limit cycle
We consider a Poincare section $H$ defined by

$$
H =\{ \Phi : \varphi(\Phi)= \varphi(\Phi_0) \mod 2\pi \}.
$$

- $\varphi$ denotes the global phase.
- $\Phi_0$ - phase vector at initial time $t_0$.
- $H$ is a $(N-1)$-dimensional hypersurface in $N$-dimensional phase space.
- In case of $\varphi$ - mean phase, $H$ is a hyperplane.

For each of m-twists solutions, we anticipate a corresponding limit cycle $C_\mathbf{k}$, piercing $H$ close to $\Phi_k$ i.e.

$$
C_\mathbf{k} \cap H_0 = \Phi_\mathbf{k} + \mathbf{E}^*,
$$

with a small correction vector $\mathbf{E}^*$.
The reason for the small correction $\mathbf{E}^*$ is that the calibration of active driving forces causes small phase-dependent variations of the instantaneous phase speed *(a cilium will speed up at one part of the beat cycle, and slow down at another; but then the phase difference will vary during the cycle, therefore we don't get a perfect m-twist, which has constant phase difference)*.

Therefore, the first step is to find $\Phi^*_\mathbf{k}=\Phi_\mathbf{k} + \mathbf{E}^*$.

We define **Poincare map** $\mathcal{L}: H \rightarrow H$ as
 $$\mathcal{L}(\Phi_0) = \Phi_1,$$
where
- $\Phi_0 = \Phi(t_0)$ - some initial phase. *Note that for any $\Phi_0$ we can define a Poincare section as defined above.*
- $\Phi_1=\Phi(t_1)$ - where $t_1$ is the next time when our phase trajectory hits the Poincare section.

Then  $\Phi^*_\mathbf{k}$ is a fixpoint of Poincare map (
**TODO:** explain why it cannot be another limit cycle in Poincare plane.
)
#### Procedure to find the fixpoint $\Phi^*$
By definition, the fixpoint is such a point that
$$
\mathcal{L}(\Phi^*) = \Phi^*
$$
We define (**TODO:** this function can be useful in other places - keep this notation?)
$$
D(\Phi) = \mathcal{L}(\Phi) - \Phi
$$
and
$$
d(\Phi) = \lVert D(\Phi) \rVert ^ 2
$$
Function $d(\Phi)$ is the squared distance between a phase vector and its Poincare map image.

Properties of norm imply that
- $d(\Phi) \geq 0$
- $d(\Phi)=0 \iff \Phi = \Phi^*$ - a fixed point

Therefore, to find the fixpoint $\Phi^*_{\mathbf{k}}$, we numerically find a minimum of function $d$ with initial condition in $\Phi_{\mathbf{k}}$. Used library function `scipy.optimize.minimize` with method `BFGS`.


### Linearized Poincare map
- $\mathbf{L} = \mathrm{D}\mathcal{L}(\mathbf{\Phi^*})$ - linearized Poincare map at a fixed point (matrix). In code `Lmat`

$\mathbf{L} = e^\mathbf{\Lambda}; \quad \mathbf{\Lambda} = \log \mathbf{L} $


- $\Lambda$ - is not a symmetrical matrix.

---
# Simulations

### Eigenvectors
- Stored as **columns** of NxN array `evecs`
- Normalized
- Not orthogonal in general, but a lot of them *are* orthogonal
- There are eigenvectors which are complex conjugates of one another (which is what we expect for complex eigenvalues).
- Perturbation made of eigenvector and its complex conjugate will develop after a cycle as predicted to linear theory up to precision of `10 ** -4` - `10 ** -5`.
  - This is not very small, but it remains as small if I increase `delta0` to `10 ** -1`
  - Imaginary part of eigenvalue gives very small contribution.

- Most of eigenvectors have their components lying on a circle - which is often shifted away from origin.

#### Eigenvector decomposition

Decompose eigenvector into basis of complex exponents of m-twists (coefficients are essentially multidimensional discrete Fourier transform output):
$$
\Delta = \sum_{\mathbf{k}} d_k e^{i \Phi_\mathbf{k}}.
$$
Observed (1D carpet,`try09b`, `try12a`) that in fact only two components give major contribution: $\mathbf{k}=\mathbf{0}$ and $\mathbf{k}=\mathbf{k_1} \neq \mathbf{0}$,

$$
\Delta = d_0 +  d_{\mathbf{k_1}} e^{i \Phi_\mathbf{k_1}}  +\mathbf{R_2},
$$

where $\mathbf{R_2} = \sum_{\mathbf{\mathbf{k} \neq \mathbf{0}, \mathbf{k_1}}} d_k e^{i \Phi_\mathbf{k}}$ - residual.

- In most of the cases $\lVert \mathbf{R_2} \rVert < 0.02$. In a few cases $ 0.02 < \lVert \mathbf{R_2} \rVert < 0.1$ and taking just one additional term would lower the residual down to $ 0.02$.

- Representation in the form
$$
\Delta = d_0 +  d_{\mathbf{k_1}} e^{i \Phi_\mathbf{k_1}}
$$
means that eigenvector components are arranged into a circle:
<img alt="carpet_jupyter_summary-2019-07-03-af1e3af7.png" src="assets/carpet_jupyter_summary-2019-07-03-af1e3af7.png" width="" height="" >

- Circle centers $d_0$ are real or almost real.
- Eigenvectors which are complex conjugates of each other have their centers complex conjugated as well, therefore centers aare distributed symmetrical around real axis.

We can construct real perturbation based on complex eigenvectors:
$$
\widetilde \Delta = \frac{1}{2} ( \Delta + \Delta^*) = \operatorname{Re}(\Delta)
= \operatorname{Re}(d_0) + \operatorname{Re}(d_{\mathbf{k_1}}) \cos(\Phi_{\mathbf{k_1}}) - \operatorname{Im}(d_{\mathbf{k_1}}) \sin(\Phi_{\mathbf{k_1}})
$$

$$
\widetilde{\widetilde \Delta} = \frac{1}{2} ( \Delta - \Delta^*) = \operatorname{Im}(\Delta)
= \operatorname{Im}(d_0) + \operatorname{Im}(d_{\mathbf{k_1}}) \cos(\Phi_{\mathbf{k_1}}) + \operatorname{Re}(d_{\mathbf{k_1}}) \sin(\Phi_{\mathbf{k_1}})
$$