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
- $\phi(\mathbf{\Phi})$ - global phase, according to one of the definitions. *Currently mean phase is the default choice for global phase.*
- $\mathcal{L}: H \rightarrow H$ - Poincare map
- $\mathbf{\Phi^*} \in H$ - fix point of the Poincare map
- $\mathbf{\Phi^{*}_k}$ - fix point close to the m-twist solution
- $\mathbf{\Delta_0}$, $\mathbf{\Delta_1}$ - perturbation initial and after one cycle: $\mathcal{L}(\mathbf{\Phi^{*}} + \mathbf{\Delta_0}) = \mathbf{\Phi^{*}} + \mathbf{\Delta_1}$
- $\mathbf{L} = \mathrm{D}\mathcal{L}(\mathbf{\Phi^*})$ - linearized Poincare map at a fixed point (matrix). In code `Lmat`
- $\mathbf{L} = e^\mathbf{\Lambda}; \quad \mathbf{\Lambda} = \log \mathbf{L} $ - logarithm of the linearized Poincare map. In code `Lmat_log`
- $\lambda_j$ - eigenvalues of $\Lambda$

**Conflicts:**
- Fixpoint notation and complex conjugation.
- Once $\delta$ and $\Delta$ are taken - how to denote a difference of some values?
- global phase and components of $\mathbf{\Phi}$
- $d(\Phi)$ in procedure to find fixpoint and $d$ - as a distance between cilia. Change $d$ to $a$ - spacing between cilia?
- TODO: $\mathbf{\Delta_0}$- > $\mathbf{\Delta}_0$ -


**Fill the tables to see which letters I can still use?**

|         |   |         |   |
|---------|---|---------|---|
| a       |   | A       |   |
| b       |   | B       |   |
| c       |   | C       |   |
| d       |   | D       |   |
| e       | basis  | E       | Ben's distance from mtwist to fixpoint  |
| f       |   | F       |   |
| g       |   | G       |   |
| h       |   | H       |   |
| i       | idx  | I       |   |
| j       | idx  | J       |   |
| k       | idx, mtwist number, dual space vector  | K       |   |
| l       | idx | L       |   |
| m       |   | M       |   |
| n       | $n_x$, $n_y$ - num cilia in 1 direction  | N       |  num cilia |
| o       |   | O       |   |
| p       |   | P       |   |
| q       |   | Q       |  gen force |
| r       |   | R       |   |
| s       |   | S       |   |
| t       | time, translation vec  | T       | period  |
| u       |   | U       |   |
| v       |   | V       |   |
| w       |   | W       |   |
| x       | coords  | X       |   |
| y       | coords  | Y       |   |
| z</pre> | coords  | Z</pre> |   |
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
where $k_1,k_2\in{Z}$ and $a_1=a$, $a_2 = \sqrt{3} a/2$

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
- Poincare surface will be a plane, defined by normal $(1,1,1,1..,1)$.
  [A language of quotient vector space could be useful; just keep this reference here for now http://mathworld.wolfram.com/QuotientVectorSpace.html]
  (quotient = частное)

Cons
- $\varphi_j$ can't jump from  $2 \pi$ to $0$, otherwise the global phase will make a jump by $- \frac{2 \pi}{N}$.

## Poicnare map for a m-twist solution
### Poincare map and limit cycle
We consider a Poincare section $H$ defined by

$$
H =\{ \Phi : \varphi(\Phi) \equiv \varphi(\Phi_0) \mod 2\pi \}.
$$

- $\varphi$ denotes the global phase.
- $\Phi_0$ - phase vector at initial time $t_0$.
- $H$ is a $(N-1)$-dimensional hypersurface in $N$-dimensional phase space.
  - If $\varphi$ - mean phase, $H$ is a hyperplane.

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

Then  $\Phi^*_\mathbf{k}$ is a fixpoint of Poincare map

**TODO:** explain why it cannot be another limit cycle in Poincare plane.

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

We use linear stability analysis to study stability of limit cycles. Consider a fixpoint of Poincare map $\Phi^*$ (corresponds to a limit cycle) and apply a small perturbation $\mathbf{\Delta_0}$.

$$\mathcal{L}(\mathbf{\Phi^{*}} + \mathbf{\Delta_0}) = \mathbf{\Phi^{*}} + \mathbf{\Delta_1}$$

Power expansion in vicinity of the fixpoint yields

$$\mathcal{L}(\mathbf{\Phi^{*}} + \mathbf{\Delta_0}) = \mathbf{\Phi^{*}} +  \mathrm{D}\mathcal{L}(\mathbf{\Phi^*}) \mathbf{\Delta_0} + \mathcal{O}(\lVert\mathbf{\Delta_0}\rVert^2) $$

where $\mathrm{D}\mathcal{L}(\mathbf{\Phi^*})$ is linear contribution (represented by a matrix), also known as linearized Poincare map [reference].

Further let's use short notation $\mathbf{L}$ = $\mathrm{D}\mathcal{L}(\mathbf{\Phi^*})$. In code `Lmat`

TODO:

$\mathbf{L} = e^\mathbf{\Lambda}; \quad \mathbf{\Lambda} = \log \mathbf{L} $


- $\mathbf{\Lambda}$ - is not a symmetrical matrix.

## Numerical inaccuracy discussion

### Inaccuracy in fixpoint

#### Fixpoint

We find fixpoint by optimization procedure [defined above] up to a tolerance, specified in code. Two tolerances are involved: (i) solver tolerance and (ii) minimizer tolerance. After minimal tests, those were taken to be equal and further in this section both are referred simply as tolerance.

Suppose we found numerically a fixpoint $\widetilde \Phi^*$, but it is not perfect and $ \mathcal{L}(\widetilde \Phi^*) \neq \widetilde \Phi^* $:

$$ \mathcal{L}(\widetilde \Phi^*) = \widetilde \Phi^* + \Delta_*$$

This means, that real fixpoint is somewhere else:

$$\widetilde \Phi^*= \Phi^* + \widetilde \Delta$$

If we consider linear approximation in vicinity of $\Phi^*$


$$ \mathcal{L}(\widetilde \Phi^*) = \Phi^* + \mathbf{L}  \widetilde \Delta$$

$$ \widetilde \Phi^* + \Delta_* = \Phi^* + \mathbf{L}  \widetilde \Delta$$

$$  \Phi^* + \widetilde \Delta + \Delta_* = \Phi^* + \mathbf{L}  \widetilde \Delta$$

Therefore

$$
\Delta_* = ( \mathbf{L} - \mathbf{I} ) \widetilde \Delta + \mathcal{O}(\lVert \widetilde \Delta \rVert^2)
$$

---
Note
- We don't know real $\mathbf{L}$.
- We don't know $\widetilde \Delta$ - distance to real fixpoint.
- We can, however, calculate $\Delta_*$.
---

Let's denote
 $\varepsilon_{F}
 = \lVert \Delta_* \rVert
 =\lVert \mathcal{L}(\widetilde \Phi^*) - \widetilde \Phi^*  \rVert $
and $\delta_F = \lVert \widetilde \Delta \rVert $. Then

$$ \lVert \Delta_* \rVert = \lVert ( \mathbf{L} - \mathbf{I} ) \widetilde \Delta \rVert$$


$$
\min|\lambda| \delta_F
\leq \varepsilon_{F}
\leq \max|\lambda| \delta_F
$$

$$
\varepsilon_{F} \sim |\lambda| \delta_F
$$



#### Perturbed state
Let's consider a perturbed (approximate) fixpoint. Fixpoint error will add up to the perturbation:

$$ \mathcal{L}(\widetilde \Phi^* + \widetilde \Delta_0)
= \mathcal{L}(\Phi^* + \widetilde  \Delta +  \Delta_0 )
= \Phi^* +  \mathbf{L}  ( \widetilde \Delta +  \widetilde \Delta_0 )
$$
The right side denote as initial state plus deviation

$$
\mathcal{L}(\widetilde \Phi^*  + \widetilde \Delta_0)
= \widetilde \Phi^* + \widetilde \Delta_1
= \Phi^* + \widetilde \Delta  + \widetilde \Delta_1
$$

And combining these two equalities

$$
\widetilde \Delta_1
=  \mathbf{L}   \widetilde \Delta_0  + (\mathbf{L} - \mathbf{I}) \widetilde \Delta
$$

$$
\widetilde \Delta_1 - \widetilde \Delta_0
=  (\mathbf{L} - \mathbf{I})  \widetilde \Delta_0  + (\mathbf{L} - \mathbf{I}) \widetilde \Delta
$$

$$
\widetilde \Delta_1 - \widetilde \Delta_0
=  (\mathbf{L} - \mathbf{I})  \widetilde \Delta_0  + \Delta_*
$$

Now let's estimate norms

$$
\lVert \widetilde \Delta_1 - \widetilde \Delta_0 \rVert
\sim
|\lambda|(\delta_0 + \delta_F)
\sim
|\lambda|\delta_0 + \varepsilon_F
$$

$$
\frac{ \lVert \Delta_* \rVert }  {\lVert \widetilde \Delta_1 - \widetilde \Delta_0 \rVert }
\sim
\frac{ \varepsilon_F } { |\lambda| \delta_0 + \varepsilon_F}
$$

If $|\lambda| \delta_0 >>  \varepsilon_{F}:$

$$
\frac{ \lVert \Delta_* \rVert }  {\lVert \widetilde \Delta_1 - \widetilde \Delta_0 \rVert }
\sim
\frac{ \varepsilon_F } { |\lambda| \delta_0}
\sim
\frac{ \delta_F } {  \delta_0 }
$$


##### Conclusions
- We want to make sure that $\frac{ \lVert \Delta_* \rVert }  {\lVert \widetilde \Delta_1 - \widetilde \Delta_0 \rVert } $ is small. Otherwise, $\widetilde \Delta_1 - \widetilde \Delta_0$ is significantly affected by $\Delta_*$ contribution, and we will have big error when we calculate matrix $\mathbf{L}$. TODO: quantify

- $\widetilde \Delta_1 - \widetilde \Delta_0$ has a constant contribution equal to $\Delta_*$.
- Note that $\widetilde \Delta_0$, $\widetilde \Delta_1$ - simulation input and outputs, so that's something we can directly obtain.
- $ \Delta_* $ - is something that we can measure, and can control indirectly, by tuning
  algorithm tolerance.
- Therefore, we must make sure that
  - Test at $N=6$ showed that this ratio is around or below $10^{-2}$ if fixpoint tolerance is $10^{-8}$.

- Also, just for fun, let's note
$$
 \mathcal{L}(\widetilde \Phi^*  + \widetilde \Delta_0) - \mathcal{L}( \Phi^*  + \widetilde \Delta_0) = D\mathcal{L} \vert_{\widetilde \Phi^*  + \widetilde \Delta_0} \widetilde \Delta + \mathcal{O}(\delta_F^2)
$$


#### Linearized map
**TODO**: double-check with Ben; is there a better estimation?

Expand Poincare map $\mathcal{L}$in $\Phi^*$

$$ \mathcal{L}(\widetilde \Phi^* + \widetilde \Delta + \widetilde \Delta_0)
= \Phi^* +  \mathbf{L}  ( \widetilde \Delta +  \widetilde \Delta_0 ) + \mathcal{O}(( \delta_F + \delta_0)^2) = ... + \mathcal{O}(\delta_F^2, \delta_F \delta_0,  \delta_0^2)
$$
and in $\widetilde \Phi^*$
$$
\mathcal{L}(\widetilde \Phi^* + \widetilde \Delta + \widetilde \Delta_0)
= \mathcal{L}(\widetilde \Phi^*) +  \widetilde \mathbf{L} \widetilde \Delta_0 + \mathcal{O}(\delta_0^2) =
$$
$$
= \widetilde \Phi^* + \Delta_*   +  \widetilde \mathbf{L} \widetilde \Delta_0 + \mathcal{O}(\delta_0^2)
= \Phi^* + \widetilde \Delta + \Delta_* +  \widetilde \mathbf{L} \widetilde \Delta_0 + \mathcal{O}(\delta_0^2)
$$

Combining those two
$$
\Phi^* +  \mathbf{L}  ( \widetilde \Delta +  \widetilde \Delta_0 )
= \Phi^* + \widetilde \Delta + \Delta_* +  \widetilde \mathbf{L} \widetilde \Delta_0 + \mathcal{O}(\delta_F^2, \delta_F \delta_0, \delta_0^2)
$$

$$
 \mathbf{L} \widetilde \Delta_0
=  \widetilde \mathbf{L} \widetilde \Delta_0  + \Delta^* - (\mathbf{L} - I)\widetilde \Delta + \mathcal{O}(\delta_F^2, \delta_F \delta_0, \delta_F \delta_0, \delta_0^2)
$$
Since $ \Delta^* - (\mathbf{L} - I)\widetilde \Delta = \mathcal{O}(\delta_F^2)$,

$$
 \mathbf{L} \widetilde \Delta_0
=  \widetilde \mathbf{L} \widetilde \Delta_0  + \mathcal{O}(\delta_F^2, \delta_F \delta_0, \delta_0^2)
$$

$$
( \mathbf{L} -  \widetilde \mathbf{L} ) \widetilde \Delta_0  = \mathcal{O}(\delta_F^2, \delta_F \delta_0, \delta_0^2)
$$


$$
 \mathbf{L} -  \widetilde \mathbf{L}
 = \mathcal{O}(\delta_F^2, \delta_F \delta_0, \delta_0^2)  = \mathcal{O}(\delta_F)
$$
#### Numerical estimations
If $N = 6$, fixpoint with `tol = 10 ** -8`:

$|\lambda| = 10 ^ {-2}$ - $10 ^ {-3}$

$\varepsilon_F = 10 ^ {-6}$ - $10^{-7}$

Therefore $\delta_F \sim 10 ^{ -4}$

If $\delta_0 = 10 ^{-3}$, $\frac{ \delta_F } {  \delta_0 } = 10^{-1}$, and indeed,
$ \frac{ \lVert \Delta_* \rVert }  {\lVert \widetilde \Delta_1 - \widetilde \Delta_0 \rVert }$ lie in range $2 * 10^{-2} - 10 ^{-3} $

Goals/questions:
- What can we do with the accuracy that we have? Which parameters matter?
- Estimate distance to the real fixpoint: YES
- What order of errors we get in $\widetilde \mathbf{L}$ compared to $\mathbf{L}$
  - Done: due to error in fixpoint
  - finite perturbation length

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


[todo why  there are values on a line?]



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
