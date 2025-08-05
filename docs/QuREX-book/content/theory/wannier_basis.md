# Wannier representation
We refer to {math}`|u_{n\mathbf{k}}>` as the Bloch function computed by the DFT code (i.e. QE).
We distinguish between Bloch-like functions in the `Wannier basis`, **Wannier gauge** (which for our purpose can be any localized basis set), and the ones expressed in a Bloch-like band basis **Hamiltonian gauge**.

- **Wannier gauge**: 
  
  {math}`
  |u^{W}_{n\mathbf{k}}> = \sum_{\mathbf{R}} e^{-i \mathbf{k} \cdot(\hat{\mathbf{r}}-\mathbf{R})}|\mathbf{R} n\rangle
  {math}` (eq-wannier-bloch-function)
  
- **Hamiltonian gauge**:
  
  {math}`
    \left|u_{n \mathbf{k}}^{(\mathrm{H})}\right\rangle=\sum_m\left|u_{m \mathbf{k}}^{(\mathrm{W})}\right\rangle U_{m n}(\mathbf{k})
  {math}` (eq-hamiltonian-bloch-function)

where {math}`U` is the unitary matrix that diagonalizes {eq}`eq-wannier-bloch-function` 

# Representations in band theory
In 1962 Blount {cite}`blount1962formalisms` reviewed the formalism of band theory, pointing out the existence of different representations for electronic states: the crystal momentum representation (CMR) developed by Adams, the Kohn-Luttinger CMR or Modified CMR (MCMR), and the Wannier one.
In this work, he focused on the Schroedinger, Pauli and Dirac Hamiltonians but his results are easily generalizable for applications with DFT Hamiltonians:

 1) Schrödinger: {math}`H=\frac{p^2}{2 m}+U`
 2) Pauli: {math}`H=\frac{p^2}{2 m}+\frac{e^2}{4 m^2 c^2} \mathbf{p} \cdot \mathbf{\sigma} \times \nabla \mathrm{U}+U`
 3) Dirac: {math}`H=c \boldsymbol{\alpha} \cdot \mathrm{p}+U`

The velocity operator is defined as:

{math}`
\mathfrak{B}=-\frac{i}{\hbar}[\mathbf{x}, H]=\nabla_{\mathbf{p}}H
{math}`

When U is a periodic potential, Bloch's theorem applies and we classify eigenfunctions {math}`\psi_{n\mathbf{k}}(\mathbf{x})` by ({math}`\mathbf{k}, n`) quantum numbers

{math}`
\psi_{n \mathbf{k}}(\mathbf{x})=e^{(i \mathbf{k} \cdot \mathbf{x})} u_{n \mathbf{k}}(\mathbf{x})
{math}` (eq:blochwf)
with {math}`u_{n\mathbf{k}}` the periodic part of the Bloch function.

We write an EOM for {math}`u` in the form
{math}`
H(\mathbf{k}) u_{n \mathbf{k}}(\mathbf{x})=E_n(\mathbf{k}) u_{n \mathbf{k}}(\mathbf{x})
{math}` (eq:EOM-u)

with

{math}`
H(\mathbf{k}) \equiv e^{(-i \mathbf{k} \mathbf{x})} H e^{(i \mathbf{k} \cdot \mathbf{x})}
{math}`

The solutions of {eq}`eq:EOM-u` are peridic. Hence, {math}`\psi_{n\mathbf{k}}` and the functions {math}`\psi_{n\mathbf{k+K}}` span the same space and it is enough
to restrict ourself to the first unit cell in {math}`\mathbf{k}` k space.

## Wavefunction representations
Any wavefunction {math}`f(\mathbf{x})` can be expressed as a suporposition of Bloch functions with {math}`f_{n}(\mathbf{k})` the wavefunction in the CMR.

{math}`
f(\mathbf{x})=\sum_n \int d^3 k f_n(\mathbf{k}) \psi_{n \mathbf{k}}(\mathbf{x})
{math}`

### Crystal momentum
The crystal momentum operator is {math}`\mathbf{p_e} = \hbar \mathbf{k}`, with matrix elements

{math}`
\mathbf{p}_{e,nn^\prime}(\mathbf{k},\mathbf{k^\prime})
{math}` (eq:crystalmomentumoperator)

while the true momentum has the form

{math}`
\mathbf{p}_{n n^{\prime}}\left(\mathbf{k}, \mathbf{k}^{\prime}\right)=\delta\left(\mathbf{k}-\mathbf{k}^{\prime}\right)\left(\hbar \mathbf{k} \delta_{n n^{\prime}}-i \hbar \int_{uc} u_n{ }^* \frac{\partial u_{n^{\prime}}}{\partial \mathbf{x}} d \tau\right)
{math}`

the velocity operator has matrix elements {math}`\mathfrak{B}_{nn\prime}(\mathbf{k})` ({math}`n=n\prime` are intrabands, {math}`n\neq n^\prime` are interband)$

### Position representation
Representation of {math}`\mathbf{x}` involevs evaluating the following integral

{math}`
\begin{align}
I_{n^{\prime} \mathbf{k}^{\prime} n \mathbf{k}} & =\int \psi_{n^{\prime} \mathbf{k}^{\prime}}^* \mathbf{x} \psi_{n \mathbf{k}} d^3 \mathbf{x} \nonumber\\
&=\delta_{n n^{\prime}} \sum_R e^{\left[i\left(\mathbf{k}-\mathbf{k}^{\prime}\right) \cdot \mathbf{R}\right]} \mathbf{R}+ \nonumber \\
&\sum e^{\left[i\left(\mathbf{k}-\mathbf{k}^{\prime}\right) \cdot \mathbf{R}\right]} \xi_{n^{\prime} n}(\mathbf{k})
\end{align}
{math}` (eq:posrepr1)

where we used the usual trick of replacing the integration over the whole cell with the integration over the unit cell.

{math}`
\mathbf{\xi}_{n^{\prime} n}(k)=\int_{uc} u_{n^{\prime} \mathbf{k}}^* \mathbf{x} u_{n \mathbf{k}} d \tau
{math}`

{eq}`eq:posrep1` has two problems. The first term is not well defined and the second term depends on the choice of the unit cell.
A different approach could be to write

{math}`
\begin{align}
I_{n^{\prime} \mathbf{k}^{\prime} n \mathbf{k}} & =\int \psi_{n^{\prime} \mathbf{k}^{\prime}}^* x^\mu \psi_{n \mathbf{k}} d^3 \mathbf{x} \nonumber\\
&=-i \frac{\partial}{\partial k^\mu} \int \psi_{n^{\prime} \mathbf{k}^{\prime}}^* \psi_{n \mathbf{k}} d^3 x \\
& \quad+\int u_{n^{\prime} \mathbf{k}^{\prime}}^* e^{\left[i\left(\mathbf{k}-\mathbf{k}^{\prime}\right) \cdot \mathbf{x}\right]} \frac{i \partial u_{n \mathbf{k}}}{\partial k^\mu} d^3 x \\
& =-i \frac{\partial}{\partial k^\mu} \Delta_{n^{\prime} n}\left(\mathbf{k}^{\prime}, \mathbf{k}\right)+\delta\left(\mathbf{k}-\mathbf{k}^{\prime}\right) \mathfrak{X}_{n^{\prime} n}^\mu(\mathbf{k})
\end{align}
{math}` (eq:posrepr2)

where {math}`\Delta_{n^{\prime} n}\left(\mathbf{k}^{\prime}, \mathbf{k}\right)=\int \psi_{x^{\prime} \mathbf{k}^{\prime}}^* \psi_{n \mathbf{k}} d^8 x` (cannot be assumed to be {math}`\delta` because we need to differentiate) and {math}`\mathfrak{X}_{n^{\prime} n}(\mathbf{k})=\int u_{n^{\prime} \mathbf{k}}^* i \frac{\partial u_{n \mathbf{k}}}{\partial \mathbf{k}} d \tau`.

{math}`\mathfrak{X}` is not sensitive to the choice of the unit cell but {math}`\psi_{n\mathbf{k}}` is not differential w.r.t. {math}`\mathbf{k}`

For all pratical calculation representation {eq}`eq:posrepr2` is used.
The action of {math}`\mathbf{x}` on a generic wave function {math}`f_n(\mathbf{k})` is given by:

{math}`
\mathbf{x} = i\frac{\partial}{f_n}{\mathbf{k}} + \sum \mathbf{\mathfrak{X}}_{nn^\prime}f_{n^\prime}
{math}` (eq:xoff)

The arbitrariness in {math}`\mathfrak{X}` arises from the arbitrary phase phactor {math}`e^{-i\phi(\mathbf{k})}`in the {math}`u's`: {math}`i\frac{\partial u}{\partial\mathbf{k}}\rightarrow e^{-i\phi}(i\frac{\partial u}{\partial{\mathbf{k}}} + u\frac{\partial \phi}{\partial \mathbf{k}})`
Under this phase transformation {math}`\mathfrak{X}_{nn^\prime}` does not transform as an operator for {math}`n = n^\prime`.

{math}`
\begin{align}
& \mathfrak{X}_{n n^{\prime}}^{\prime}=\exp \left(i \varphi_n\right) \mathfrak{X}_{n n^{\prime}} \exp \left(-i \varphi_{n^{\prime}}\right) \quad n \neq n^{\prime} \\
& \mathfrak{X}_{n n}^{\prime}=\mathfrak{X}_{n n}+\frac{\partial \varphi_n}{\partial \mathbf{k}}
\end{align}
{math}`

However, a compensantory term in the first term of {eq}`eq:xoff` makes it in such a way that {math}`\mathbf{x_c} = i\partial/\partial\mathbf{k} + \mathfrak{X}_{nn}` transforms like an operator. So that {math}`\mathbf{x} = `\mathbf{x}_c{math}` + X` is a well defined position operator.
{math}`X` has only interband matrix elements equal to {math}`\mathfrak{X}_{nn^\prime}` ({math}`n\neq n^\prime`)

# References

```{bibliography}
