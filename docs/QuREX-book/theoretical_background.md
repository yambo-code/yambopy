# Theoretical Background

This section includes the theoretical formalism needed to understand and contribute to the QuREX library:

- [model Hamiltonian](model_hamiltonian)
- [excitonic Hamiltonian](h2p)
- [Wannier basis](wannier_basis)
- [BSE equation](bse_equation)
- [Wannier Exciton](wannier_exciton)

# Introduction
Density Functional Theory (DFT){cite}`kohn1965self` is the workhorse for materials simulation in nanoscale systems. It allows for a quantitative microscopic description of the ground state properties of a material from first principles, where only the crystal structure and the atomic species are provided as an input to the code.
The excited state cannot be described by DFT, and the Ab Initio Many Body Perturbation theory (Ai-MBPT){cite}`onida2002electronic` has been developed, over the years, to predict the physical properties of excited systems under an external perturbing potential.
The Ai-MBPT exploits perturbation techniques based on a Green's function approach, for the description of the propagation of particles. Theoretical spectroscopy studies the interaction of matter with light. State of the art calculations employs the GW approximation and Bethe-Salpeter equation (BSE) to describe the absorption spectra of materials, and these methods have been developed in open-source codes such as `yambo` {cite}`sangalli2019many,marini2009yambo`.
Simulating complex systems (such as van der Waals heterostructures) requires large simulation cells (hundred atoms), hence becoming extremely demanding in terms of computational time and memory. These challenges can be overcome either with the advent of the exascale era (1018 double precision floating point operations per second) or through the development of more cost-effective theoretical models that maintain high accuracy and precision, such as Wannierized models ~{cite}`marzari2012maximally`. The Wannierization technique offers a refined description of quantum Hilbert space using a subset of localized basis functions, significantly reducing computational expenses. While many ground state DFT codes currently integrate the Wannierization process, a similar integration is desirable for first principles ab initio codes that treat excited systems and are capable of including correlations and excitonic effects~{cite}`haber2023maximally`.
The description of excitonic effect requires constructing and solving the two-particle Hamiltonian ($H_{2p}$, composed of the difference in energy between the valence and conduction electronic states and the many body kernel $K$.
The many body kernel account for the direct and exchange Coulomb interaction, screened by the dielectric environment. The computation of the dielectric screening function from first-principles is computational expensive in terms of memory and computing time. Approaches, based on model dielectric functions and Coulomb potentials, have been employed in the literature but we find out to fail for the description of complex systems such as 2D van der Waals (vDW) transition metal dichalcogenides heterostructures (TMDs HS).
In this manuscript, we review the theoretical framework of maximally localized exciton Wannier functions (MLXWF) and describe the implementation done in the [yambopy Python library](https://github.com/rreho/yambopy). Our method allows to build and solve the ${H_{2p}}$ via model Coulomb potentials or extracting the many body kernel ${K}$ from first-principles. 

# Gauge issues in Non Linear Optical Responses
In the study of NLOR the long-range wavelength limit is assumed. This implies that the spatial dependence of the radiation electric field is neglected.
Hence, one has to deal with two representation of the radiation field {cite}`ventura2017gauge`:

1) Velocity gauge or vector potential approach
$$
\mathbf{E}(t) = -\frac{\partial \mathbf{A}}{\partial t} 
$$ (eq:vel-gauge)
  - ✅ keeps crystal translational symmetry
  - ❌ plenty of numerical divergencies that have been shown to be 0

2) Length gauge
$$
V(r) = e\mathbf{E}(t)\cdot\mathbf{r}
$$ (eq:length-gauge)
  - ❌ breaks crystal translational symmetry
  - ✅ no numerical divergencies and used in actual calculation

The two gauges are related by unitary transformation $\mathcal{U}(t)$:

$$
\mathcal{U}(t)=\exp \left[i \frac{e}{\hbar} \int d^d \mathbf{r} \mathbf{A}(t) \cdot \mathbf{r} \rho(\mathbf{r})\right],
$$ (eq:Uvel-lengthgauge)

whcih can be used to go from vector/velocity $A$ to length/dipole $E$ gauges
$$
\mathcal{U}(t) H_A(t) \mathcal{U}^{\dagger}(t)+i \hbar \frac{d \mathcal{U}(t)}{d t} \mathcal{U}^{\dagger}(t)=H_E(t)
$$

$$
O_E:=\mathcal{U}(t) O_A(t) \mathcal{U}^{\dagger}(t)
$$

# References

```{bibliography}

