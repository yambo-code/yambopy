# Theoretical Background

## Overview

This section provides the comprehensive theoretical foundation needed to understand and contribute to the QuREX library and beyond. The formalism presented here bridges fundamental many-body theory with practical computational implementations for excitonic systems.

## Chapter Organization

The theoretical framework is organized into interconnected topics:

- **[Model Hamiltonian](model_hamiltonian)**: Construction of effective Hamiltonians for excitonic systems
- **[Excitonic Hamiltonian](h2p)**: Two-particle formalism and electron-hole interactions  
- **[Wannier Basis](wannier_basis)**: Localized basis functions for real-space analysis
- **[BSE Equation](bse_equation)**: Bethe-Salpeter equation for optical excitations
- **[Wannier Excitons](wannier_exciton)**: Real-space representation of excitonic states
- **[Quantum Wells](quantum_well)**: Confined systems and dimensional effects
- **[Exciton-Phonon Coupling](exciton_phonon_coupling)**: Vibrational interactions with excitons
- **[Group Theory Analysis](exciton_group_theory)**: Symmetry properties of excitonic states

## Introduction to Many-Body Theory for Excitons

### From Ground State to Excited States

Density Functional Theory (DFT){cite}`kohn1965self` serves as the foundation for materials simulation in nanoscale systems, providing a quantitative microscopic description of ground state properties from first principles. However, the description of excited states—particularly optical excitations and excitons—requires going beyond the single-particle picture of DFT.

The **Ab Initio Many-Body Perturbation Theory (AI-MBPT)**{cite}`onida2002electronic` framework has been developed to predict the physical properties of excited systems under external perturbations. This approach exploits perturbation techniques based on Green's function methods to describe particle propagation and interactions in the excited state manifold.

### Theoretical Spectroscopy Framework

Theoretical spectroscopy investigates the interaction between matter and electromagnetic radiation. State-of-the-art calculations employ:

1. **GW Approximation**: For accurate quasiparticle energies and band gaps
2. **Bethe-Salpeter Equation (BSE)**: For optical absorption spectra including excitonic effects
3. **Advanced Kernels**: For electron-hole interactions and screening effects

These methods have been implemented in sophisticated open-source codes such as **Yambo**{cite}`sangalli2019many,marini2009yambo`, providing the computational foundation for modern excitonic calculations.

### Computational Challenges and Solutions

#### The Scale Problem

Simulating complex systems such as van der Waals heterostructures requires large simulation cells containing hundreds of atoms, making calculations extremely demanding in terms of:
- **Computational Time**: Scaling with system size and number of states
- **Memory Requirements**: Storage of wavefunctions and interaction matrices
- **Algorithmic Complexity**: Efficient handling of many-body interactions

#### Wannierization Strategy

These challenges can be addressed through the development of cost-effective theoretical models that maintain high accuracy, particularly **Wannierized models**{cite}`marzari2012maximally`. The Wannierization technique provides:

- **Localized Basis Functions**: Refined description of quantum Hilbert space
- **Computational Efficiency**: Significant reduction in computational expenses
- **Physical Insight**: Real-space understanding of electronic properties
- **Transferability**: Models applicable across different conditions

While ground-state DFT codes commonly integrate Wannierization, similar integration is highly desirable for first-principles codes treating excited systems with correlations and excitonic effects{cite}`haber2023maximally`.

### The Two-Particle Problem

#### Excitonic Hamiltonian Construction

The description of excitonic effects requires constructing and solving the **two-particle Hamiltonian** $H_{2p}$, which consists of:

1. **Single-Particle Terms**: Energy differences between valence and conduction states
2. **Many-Body Kernel** $K$: Electron-hole interactions including:
   - Direct Coulomb interaction
   - Exchange interactions  
   - Screening effects from the dielectric environment

#### Computational Strategies

The many-body kernel accounts for Coulomb interactions screened by the dielectric environment. Computing the dielectric screening function from first principles is computationally expensive, leading to two main approaches:

1. **First-Principles Approach**: 
   - Direct calculation of screening functions
   - High accuracy but computationally demanding
   - Suitable for small to medium systems

2. **Model-Based Approach**:
   - Model dielectric functions and Coulomb potentials
   - Computationally efficient
   - Requires careful validation for complex systems

#### Challenges in Complex Systems

Traditional model approaches often fail for complex systems such as 2D van der Waals transition metal dichalcogenide heterostructures (TMD-HS) due to:
- **Anisotropic Screening**: Different screening in different directions
- **Interface Effects**: Modified dielectric properties at interfaces
- **Quantum Confinement**: Dimensional effects on Coulomb interactions

### The QuREX Approach

#### Maximally Localized Exciton Wannier Functions

This documentation reviews the theoretical framework of **Maximally Localized Exciton Wannier Functions (MLXWF)** and describes the implementation in the [yambopy Python library](https://github.com/rreho/yambopy). Our method enables:

- **Flexible Hamiltonian Construction**: Build and solve $H_{2p}$ using various approaches
- **Model Potentials**: Efficient calculations with validated model Coulomb potentials  
- **First-Principles Integration**: Extract many-body kernel $K$ from ab initio calculations
- **Real-Space Analysis**: Understand exciton localization and spatial properties

#### Methodological Advantages

The QuREX framework provides:
- **Computational Efficiency**: Reduced computational cost through Wannierization
- **Physical Insight**: Real-space understanding of excitonic properties
- **Flexibility**: Adaptable to different materials and conditions
- **Accuracy**: Maintains first-principles accuracy where needed 

## Gauge Issues in Nonlinear Optical Responses

### The Long-Wavelength Approximation

In the study of nonlinear optical responses (NLOR), the **long-wavelength limit** is typically assumed. This approximation implies that the spatial dependence of the radiation electric field is neglected, leading to the dipole approximation where the electromagnetic field is treated as spatially uniform over the extent of the system.

This approximation necessitates careful consideration of gauge choices, as different representations of the electromagnetic field can lead to different computational approaches{cite}`ventura2017gauge`.

### Gauge Representations

#### 1. Velocity Gauge (Vector Potential Approach)

In the velocity gauge, the electric field is expressed in terms of the vector potential:

$$
\mathbf{E}(t) = -\frac{\partial \mathbf{A}(t)}{\partial t} 
$$ (eq:vel-gauge)

**Advantages:**
- ✅ **Preserves Crystal Symmetry**: Maintains translational symmetry of the crystal
- ✅ **Gauge Invariance**: Natural gauge-invariant formulation
- ✅ **Physical Interpretation**: Direct connection to current operators

**Disadvantages:**
- ❌ **Numerical Divergences**: Contains apparent divergences that must be shown to cancel
- ❌ **Computational Complexity**: More complex numerical implementation
- ❌ **Convergence Issues**: Slower convergence with respect to k-point sampling

#### 2. Length Gauge (Dipole Approach)

In the length gauge, the interaction is expressed through the position operator:

$$
V(\mathbf{r}) = e\mathbf{E}(t) \cdot \mathbf{r}
$$ (eq:length-gauge)

**Advantages:**
- ✅ **Numerical Stability**: No spurious divergences in practical calculations
- ✅ **Computational Efficiency**: Simpler numerical implementation
- ✅ **Faster Convergence**: Better convergence properties
- ✅ **Standard Implementation**: Commonly used in actual calculations

**Disadvantages:**
- ❌ **Broken Translational Symmetry**: Violates crystal translational symmetry
- ❌ **Surface Effects**: Can introduce artificial surface contributions
- ❌ **Gauge Dependence**: Requires careful treatment of gauge issues

### Gauge Transformation

#### Unitary Transformation Between Gauges

The two gauge representations are related by a unitary transformation $\mathcal{U}(t)$:

$$
\mathcal{U}(t) = \exp \left[i \frac{e}{\hbar} \int d^3\mathbf{r} \, \mathbf{A}(t) \cdot \mathbf{r} \, \rho(\mathbf{r})\right]
$$ (eq:Uvel-lengthgauge)

where $\rho(\mathbf{r})$ is the charge density operator.

#### Hamiltonian Transformation

This unitary operator transforms the Hamiltonian from velocity gauge to length gauge:

$$
H_E(t) = \mathcal{U}(t) H_A(t) \mathcal{U}^{\dagger}(t) + i\hbar \frac{d\mathcal{U}(t)}{dt} \mathcal{U}^{\dagger}(t)
$$ (eq:hamiltonian-transform)

#### Observable Transformation

Similarly, observables transform according to:

$$
O_E(t) = \mathcal{U}(t) O_A(t) \mathcal{U}^{\dagger}(t)
$$ (eq:observable-transform)

### Practical Considerations

#### Computational Strategy

In practice, most calculations use the **length gauge** due to its numerical advantages, while ensuring that:

1. **Gauge Invariance**: Physical results are independent of gauge choice
2. **Proper Limits**: Correct behavior in appropriate limits
3. **Symmetry Restoration**: Careful treatment of broken symmetries

#### Validation Approaches

To ensure reliability, calculations should:
- Compare results between different gauges when possible
- Check gauge invariance of physical observables
- Validate against experimental data
- Use appropriate convergence criteria

### Applications in QuREX

Within the QuREX framework, gauge considerations are particularly important for:
- **Optical Matrix Elements**: Calculation of transition dipole moments
- **Nonlinear Responses**: Higher-order optical processes
- **Real-Space Analysis**: Spatial distribution of optical transitions
- **Symmetry Analysis**: Group theory applications to optical properties

The choice of gauge affects both the computational implementation and the physical interpretation of results, making it a crucial consideration in excitonic calculations.

# References

```{bibliography}

