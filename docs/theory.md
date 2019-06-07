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
