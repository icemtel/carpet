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
- $ \mathbf{Q(\Phi)}$ - active driving forces vector
- $ \mathbf{\Gamma}$ - friction coefficients NxN-matrix
- $\mathbf{\Phi_k}$ - m-twist solution
- $\mathcal{L}: H \rightarrow H$ - Poincare map
- $\mathbf{\Phi^*} \in H$ - fix point of the Poincare map
- $\mathbf{\Phi^{*}_k}$ - fix point close to the m-twist solution
- $\mathcal{L}(\mathbf{\Phi^{*}} + \mathbf{\Delta_0}) = \mathbf{\Phi^{*}} + \mathbf{\Delta_1}$ - perturbation initial and after one cycle
- $\mathbf{L} = \mathrm{D}\mathcal{L}(\mathbf{\Phi^*})$ - linearized Poincare map at a fixed point (matrix). In code `Lmat`
- $\mathbf{L} = e^\mathbf{\Lambda}; \quad \mathbf{\Lambda} = \log \mathbf{L} $ - logarithm of the linearized Poincare map. In code `Lmat_log`
- $\lambda_j$ - eigenvalues of $\Lambda$



## Honeycomb/Triangular lattice

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

**Dynamic equation from force-balance equation:**
$$ \dot{\Phi} = \mathbf{\Gamma}^{-1}(\Phi)Q(\Phi)  $$
This linear system can be solved directly in Python.


### Calibration from hydrodynamic simulations

*Self-friction:*

$$ \Gamma_{ii}(\Phi) = \Gamma_{11}^{\mathrm{(loc)}}(\varphi_i, \_, d, \psi) $$
*Pair-wise coupling functions:*
$$ \Gamma_{ij}(\Phi) = 0 \quad \text{if not neighbours} $$
$$ \Gamma_{ij}(\Phi) = \Gamma_{12}^{\mathrm{(loc)}}(\varphi_i,\varphi_j; d, \psi) \quad \text{if neighbours} $$

where the translation vector pointing from cilium $i$ to cilium $j$ reads
$$[n(j)-n(i)]\,\mathbf{e}_1 + [m(j)-m(i)]\,\mathbf{e}_2 = d \left( \begin{array}{c} \sin\psi \\ \cos\psi \end{array} \right).$$

*Active driving forces (simplest case: calibration for isolated cilium):*
$$ Q_{i}(\Phi) = \Gamma_{11}^{\mathrm{(loc)}}(\varphi_i, \_; d, \psi) \omega_0 $$


## Global Phase

#### Naive phase

$$
\varphi = \mathbf{\varphi}_j \text{\quad --- for some index } j
$$

#### Circular average phase
We define order parameters
[[Kuramoto, 1984]](https://link.springer.com/chapter/10.1007/978-3-642-69689-3_7)  [[Strogatz, 2000]](https://www.sciencedirect.com/science/article/pii/S0167278900000944)
for each metachronal wave vector $\mathbf{k}$

$$S_\mathbf{k} = \left| \frac{1}{N}\sum_j \exp\left[ i ( \varphi_j - \varphi_{\mathbf{k},j} ) \right] \right|,$$

where $\varphi_{\mathbf{k},j}$ - j-th component of m-twist solution defined by vector $\mathbf{k}$.

**Global phase:**
For a fixed $\mathbf{k}$
$$\varphi(t; \mathbf{k}) = \arg \left( \frac{1}{N}\sum_j \exp\left[ i ( \varphi_j(t) - \varphi_{\mathbf{k},j} ) \right] \right) $$

#### Mean phase
Use the mean of all phases as a global phase.
$$
\varphi =  \frac{1}{N}\sum_{j=1}^{N} \varphi_j
$$

Pros
- Additive: if $\mathbf{\Phi} = \mathbf{\Phi_0}+ \mathbf{\Phi_1} $, then $\varphi = \varphi_0 + \varphi_1$
- As an implication, it is independent of which initial condition we consider.

Cons
- Must be careful with $2 \pi$-periodicity.

## Poicnare map for a m-twist solution
**(from Ben's notes)**

TODO: notation: E, $\Delta$?

We consider a Poincare plane $H_0$ defined by

$$
H_0 =\{ \Phi : \overline{\varphi}(\Phi)=0 \}.
$$

Here, $\overline{\varphi}_\mathbf{k}$ denotes the global phase.
This is a $(N-1)$-dimensional in $N$-dimensional phase space.
We anticipate a limit cycle $C_\mathbf{k}$ for each m-twist solution with wave vector $\mathbf{k}$,
piercing $H_0$ close to
$\Phi'_\mathbf{k} = \Phi_\mathbf{k}-\overline{\varphi}_\mathbf{k}$, i.e.

$$
C_\mathbf{k} \cap H_0 = \Phi'_\mathbf{k} + \mathbf{E}^\ast,
$$

with a small correction vector $\mathbf{E}^\ast$.
The reason for the small correction $\mathbf{E}$ is that the calibration of active driving forces causes small phase-dependent variations of the instantaneous phase speed. Therefore, a first aim is to find the exact initial conditions $\Phi'_\mathbf{k}+\mathbf{E}^\ast$. We compute $\mathbf{E}^\ast$ as a fixed point of the Poincare map and will proceed using a Newton method. Let $\mathcal{L}$ be the full Poincare map, i.e., the nonlinear map that maps a point $\Phi\in H_0$ with $\overline{\varphi}_\mathbf{k}(\Phi)=0$ to its image $\mathcal{L}(\Phi)\in H_0$ with $\overline{\varphi}_\mathbf{k}(\mathcal{L}(\Phi))=2\pi$ after a full cycle of the global phase.

Let $\Delta(\Phi)$ be change in state after a full beat cycle

$$
\Delta(\Phi) = \mathcal{L}(\Phi) - \Phi - 2\pi .
$$

We are looking for a zero of $\Delta$ close to $\Phi_\mathbf{k}$,
i.e. $\Delta(\Phi_\mathbf{k} + E^\ast )=0$.
This is equivalent to a fixed point of $\mathcal{L}$
$$
\mathcal{L}\left( \Phi_\mathbf{k} + E^\ast \right) = \Phi_\mathbf{k} + E^\ast .
$$

Each iteration step of Newton's method involves two numerical integrations

- Step 1. Compute a deviation:
$$ \Delta_n = \Delta(\Phi_\mathbf{k} + E_n),$$
where $\Phi_\mathbf{k} + E_n$ is the current estimate for the intersection of the m-twist limit cycle for wave vector $\mathbf{k}$ with the plane $H_0$ defined by $\overline{\varphi}=0$. We start with $E_0 = 0$.

- Step 2. Compute the linearized Poincare map:
We compute the linearization
$\mathrm{L}=\nabla \mathcal{L}_{|\Phi_0}$
of the Poincare map $\mathcal{L}(\Phi_0+\delta_0 D)$
for a small deviation $\delta_0 D$ of amplitude $\delta_0$ around the center $\Phi_0$, i.e.,
$$
\mathcal{L}(\Phi_0 + \delta_0 D) = \mathcal{L}(\Phi_0) + \delta_0 \, \mathrm{L}\cdot D
+ \mathcal{O}(\delta_0^2) .
$$
The matrix $\mathrm{L}$ is essentially the matrix *Lmat* computed above.
Special precaution will be needed to project $\mathrm{L}$ on $H_0$.

- Newton update:
$$ E_{n+1} = E_n + \left( \mathbb{1} - \mathrm{L} \right)^{-1} \Delta_n . $$

*Note:*
For a quick computation, we can compute the change in state for the perturbation of a single node
$$ L_{i:}=\delta_0^{-1}\,\Delta(\Phi_\mathbf{k}+E_n+\delta_0 D_i) . $$
Here, $D_i$ is the vector that has all its entries equal to $-1/N$,
except the $i$-th entry, which equals $1-1/N$.
Then, as above, it suffices to compute the row vector $L_{i:}$ for a single index $i$, and get all other rows of the $N\times N$-matrix $\mathrm{L}$ by applying a lattice transformation.
Yet, special precaution will be needed since we expect $\mathrm{L}$ to be degenerate. We should substract the arithmetic mean of $\Phi$.
