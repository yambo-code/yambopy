# Wannier representation
We refer to $|u_{n\mathbf{k}}>$ as the Bloch function computed by the DFT code (i.e. QE).
We distinguish between Bloch-like functions in the `Wannier basis`, **Wannier gauge** (which for our purpose can be any localized basis set), and the ones expressed in a Bloch-like band basis **Hamiltonian gauge**.

- **Wannier gauge**: 
  
  $$
  |u^{W}_{n\mathbf{k}}> = \sum_{\mathbf{R}} e^{-i \mathbf{k} \cdot(\hat{\mathbf{r}}-\mathbf{R})}|\mathbf{R} n\rangle
  $$ (eq-wannier-bloch-function)
  
- **Hamiltonian gauge**:
  
  $$
    \left|u_{n \mathbf{k}}^{(\mathrm{H})}\right\rangle=\sum_m\left|u_{m \mathbf{k}}^{(\mathrm{W})}\right\rangle U_{m n}(\mathbf{k})
  $$ (eq-hamiltonian-bloch-function)

where $U$ is the unitary matrix that diagonalizes {eq}`eq-wannier-bloch-function` 

# Representations in band theory
In 1962 Blount {cite}`blount1962formalisms` reviewed the formalism of band theory, pointing out the existence of different representations for electronic states: the crystal momentum representation (CMR) developed by Adams, the Kohn-Luttinger CMR or Modified CMR (MCMR), and the Wannier one.
In this work, he focused on the Schroedinger, Pauli and Dirac Hamiltonians but his results are easily generalizable for applications with DFT Hamiltonians:

 1) Schrödinger: $H=\frac{p^2}{2 m}+U$
 2) Pauli: $H=\frac{p^2}{2 m}+\frac{e^2}{4 m^2 c^2} \mathbf{p} \cdot \mathbf{\sigma} \times \nabla \mathrm{U}+U$
 3) Dirac: $H=c \boldsymbol{\alpha} \cdot \mathrm{p}+U$

The velocity operator is defined as:

$$
\mathfrak{B}=-\frac{i}{\hbar}[\mathbf{x}, H]=\nabla_{\mathbf{p}}H
$$

When U is a periodic potential, Bloch's theorem applies and we classify eigenfunctions $\psi_{n\mathbf{k}}(\mathbf{x})$ by ($\mathbf{k}, n$) quantum numbers

$$
\psi_{n \mathbf{k}}(\mathbf{x})=e^{(i \mathbf{k} \cdot \mathbf{x})} u_{n \mathbf{k}}(\mathbf{x})
$$ (eq:blochwf)
with $u_{n\mathbf{k}}$ the periodic part of the Bloch function.

We write an EOM for $u$ in the form
$$
H(\mathbf{k}) u_{n \mathbf{k}}(\mathbf{x})=E_n(\mathbf{k}) u_{n \mathbf{k}}(\mathbf{x})
$$ (eq:EOM-u)

with

$$
H(\mathbf{k}) \equiv e^{(-i \mathbf{k} \mathbf{x})} H e^{(i \mathbf{k} \cdot \mathbf{x})}
$$

The solutions of {eq}`eq:EOM-u` are peridic. Hence, $\psi_{n\mathbf{k}}$ and the functions $\psi_{n\mathbf{k+K}}$ span the same space and it is enough
to restrict ourself to the first unit cell in $\mathbf{k}$ k space.

## Wavefunction representations
Any wavefunction $f(\mathbf{x})$ can be expressed as a suporposition of Bloch functions with $f_{n}(\mathbf{k})$ the wavefunction in the CMR.

$$
f(\mathbf{x})=\sum_n \int d^3 k f_n(\mathbf{k}) \psi_{n \mathbf{k}}(\mathbf{x})
$$

### Crystal momentum
The crystal momentum operator is $\mathbf{p_e} = \hbar \mathbf{k}$, with matrix elements

$$
\mathbf{p}_{e,nn^\prime}(\mathbf{k},\mathbf{k^\prime})
$$ (eq:crystalmomentumoperator)

while the true momentum has the form

$$
\mathbf{p}_{n n^{\prime}}\left(\mathbf{k}, \mathbf{k}^{\prime}\right)=\delta\left(\mathbf{k}-\mathbf{k}^{\prime}\right)\left(\hbar \mathbf{k} \delta_{n n^{\prime}}-i \hbar \int_{uc} u_n{ }^* \frac{\partial u_{n^{\prime}}}{\partial \mathbf{x}} d \tau\right)
$$

the velocity operator has matrix elements $\mathfrak{B}_{nn\prime}(\mathbf{k})$ ($n=n\prime$ are intrabands, $n\neq n^\prime$ are interband)$

### Position representation
Representation of $\mathbf{x}$ involevs evaluating the following integral

$$
\begin{align}
I_{n^{\prime} \mathbf{k}^{\prime} n \mathbf{k}} & =\int \psi_{n^{\prime} \mathbf{k}^{\prime}}^* \mathbf{x} \psi_{n \mathbf{k}} d^3 \mathbf{x} \nonumber\\
&=\delta_{n n^{\prime}} \sum_R e^{\left[i\left(\mathbf{k}-\mathbf{k}^{\prime}\right) \cdot \mathbf{R}\right]} \mathbf{R}+ \nonumber \\
&\sum e^{\left[i\left(\mathbf{k}-\mathbf{k}^{\prime}\right) \cdot \mathbf{R}\right]} \xi_{n^{\prime} n}(\mathbf{k})
\end{align}
$$ (eq:posrepr1)

where we used the usual trick of replacing the integration over the whole cell with the integration over the unit cell.

$$
\mathbf{\xi}_{n^{\prime} n}(k)=\int_{uc} u_{n^{\prime} \mathbf{k}}^* \mathbf{x} u_{n \mathbf{k}} d \tau
$$

{eq}`eq:posrep1` has two problems. The first term is not well defined and the second term depends on the choice of the unit cell.
A different approach could be to write

$$
\begin{align}
I_{n^{\prime} \mathbf{k}^{\prime} n \mathbf{k}} & =\int \psi_{n^{\prime} \mathbf{k}^{\prime}}^* x^\mu \psi_{n \mathbf{k}} d^3 \mathbf{x} \nonumber\\
&=-i \frac{\partial}{\partial k^\mu} \int \psi_{n^{\prime} \mathbf{k}^{\prime}}^* \psi_{n \mathbf{k}} d^3 x \\
& \quad+\int u_{n^{\prime} \mathbf{k}^{\prime}}^* e^{\left[i\left(\mathbf{k}-\mathbf{k}^{\prime}\right) \cdot \mathbf{x}\right]} \frac{i \partial u_{n \mathbf{k}}}{\partial k^\mu} d^3 x \\
& =-i \frac{\partial}{\partial k^\mu} \Delta_{n^{\prime} n}\left(\mathbf{k}^{\prime}, \mathbf{k}\right)+\delta\left(\mathbf{k}-\mathbf{k}^{\prime}\right) \mathfrak{X}_{n^{\prime} n}^\mu(\mathbf{k})
\end{align}
$$ (eq:posrepr2)

where $\Delta_{n^{\prime} n}\left(\mathbf{k}^{\prime}, \mathbf{k}\right)=\int \psi_{x^{\prime} \mathbf{k}^{\prime}}^* \psi_{n \mathbf{k}} d^8 x$ (cannot be assumed to be $\delta$ because we need to differentiate) and $\mathfrak{X}_{n^{\prime} n}(\mathbf{k})=\int u_{n^{\prime} \mathbf{k}}^* i \frac{\partial u_{n \mathbf{k}}}{\partial \mathbf{k}} d \tau$.

$\mathfrak{X}$ is not sensitive to the choice of the unit cell but $\psi_{n\mathbf{k}}$ is not differential w.r.t. $\mathbf{k}$

For all pratical calculation representation {eq}`eq:posrepr2` is used.
The action of $\mathbf{x}$ on a generic wave function $f_n(\mathbf{k})$ is given by:

$$
\mathbf{x} = i\frac{\partial}{f_n}{\mathbf{k}} + \sum \mathbf{\mathfrak{X}}_{nn^\prime}f_{n^\prime}
$$ (eq:xoff)

The arbitrariness in $\mathfrak{X}$ arises from the arbitrary phase phactor $e^{-i\phi(\mathbf{k})}$in the $u's$: $i\frac{\partial u}{\partial\mathbf{k}}\rightarrow e^{-i\phi}(i\frac{\partial u}{\partial{\mathbf{k}}} + u\frac{\partial \phi}{\partial \mathbf{k}})$
Under this phase transformation $\mathfrak{X}_{nn^\prime}$ does not transform as an operator for $n = n^\prime$.

$$
\begin{align}
& \mathfrak{X}_{n n^{\prime}}^{\prime}=\exp \left(i \varphi_n\right) \mathfrak{X}_{n n^{\prime}} \exp \left(-i \varphi_{n^{\prime}}\right) \quad n \neq n^{\prime} \\
& \mathfrak{X}_{n n}^{\prime}=\mathfrak{X}_{n n}+\frac{\partial \varphi_n}{\partial \mathbf{k}}
\end{align}
$$

However, a compensantory term in the first term of {eq}`eq:xoff` makes it in such a way that $\mathbf{x_c} = i\partial/\partial\mathbf{k} + \mathfrak{X}_{nn}$ transforms like an operator. So that $\mathbf{x} = $\mathbf{x}_c$ + X$ is a well defined position operator.
$X$ has only interband matrix elements equal to $\mathfrak{X}_{nn^\prime}$ ($n\neq n^\prime$)

# References

```{bibliography}
