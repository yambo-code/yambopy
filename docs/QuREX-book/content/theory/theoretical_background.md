# Theoretical Background

## Overview

This section provides the theoretical foundation for understanding excitonic systems within the QuREX framework. The formalism bridges many-body theory with computational implementations for excitonic calculations.

## Chapter Organization

The theoretical framework covers:

- **[GW Approximation](GW)**: Quasiparticle energies and self-energy corrections
- **[BSE Equation](bse_equation)**: Bethe-Salpeter equation for optical excitations
- **[Model Hamiltonian](model_hamiltonian)**: Effective Hamiltonians for excitonic systems
- **[Coulomb Potential](coulomb_potential)**: Electron-hole interaction modeling
- **[Wannier Basis](wannier_basis)**: Localized basis functions for real-space analysis
- **[Wannier Excitons](wannier_exciton)**: Real-space representation of excitonic states
- **[Two-Particle Hamiltonian](h2p)**: Electron-hole interaction formalism
- **[Quantum Wells](quantum_well)**: Confined systems and dimensional effects
- **[Exciton-Phonon Coupling](exciton_phonon_coupling)**: Vibrational interactions
- **[Group Theory Analysis](exciton_group_theory)**: Symmetry properties
- **[Wannier-Chern](wannier_chern)**: Topological properties

## Introduction to Many-Body Theory for Excitons

### From Ground State to Excited States

Density Functional Theory (DFT){cite}`kohn1965self` provides ground state properties from first principles. However, excited states—particularly optical excitations and excitons—require many-body methods beyond the single-particle DFT picture.

**Ab Initio Many-Body Perturbation Theory (AI-MBPT)**{cite}`onida2002electronic` describes excited systems using Green's function methods and perturbation theory to account for particle interactions in the excited state manifold.

### Theoretical Spectroscopy Framework

Theoretical spectroscopy calculations employ:

1. **GW Approximation**: Quasiparticle energies and band gaps
2. **Bethe-Salpeter Equation (BSE)**: Optical absorption spectra with excitonic effects
3. **Interaction Kernels**: Electron-hole interactions and screening effects

These methods are implemented in codes such as **Yambo**{cite}`sangalli2019many,marini2009yambo`.

### Computational Challenges and Solutions

#### Computational Challenges

Large systems like van der Waals heterostructures present challenges:
- **Computational Time**: Scaling with system size and states
- **Memory Requirements**: Storage of wavefunctions and matrices
- **Algorithmic Complexity**: Many-body interaction handling

#### Wannierization Approach

**Wannierized models**{cite}`marzari2012maximally` address these challenges by providing:

- **Localized Basis Functions**: Efficient Hilbert space description
- **Computational Efficiency**: Reduced computational cost
- **Physical Insight**: Real-space electronic properties
- **Transferability**: Applicable across different conditions

Wannierization is standard in DFT codes and beneficial for excited-state calculations{cite}`haber2023maximally`.

### The Two-Particle Problem

#### Two-Particle Hamiltonian

Excitonic effects require the **two-particle Hamiltonian** {math}`H_{2p}`:

1. **Single-Particle Terms**: Valence-conduction energy differences
2. **Many-Body Kernel** {math}`K`: Electron-hole interactions including:
   - Direct Coulomb interaction
   - Exchange interactions  
   - Dielectric screening effects

#### Computational Approaches

Two strategies exist for the many-body kernel:

1. **First-Principles**: Direct screening calculation (accurate but expensive)
2. **Model-Based**: Model dielectric functions (efficient but requires validation)

#### Complex System Challenges

2D van der Waals heterostructures present difficulties:
- **Anisotropic Screening**: Direction-dependent screening
- **Interface Effects**: Modified dielectric properties
- **Quantum Confinement**: Dimensional effects on interactions

### The QuREX Approach

#### Maximally Localized Exciton Wannier Functions

The **Maximally Localized Exciton Wannier Functions (MLXWF)** framework implemented in yambopy enables:

- **Flexible Hamiltonian Construction**: Multiple approaches for building {math}`H_{2p}`
- **Model Potentials**: Efficient calculations with validated Coulomb potentials  
- **First-Principles Integration**: Extract kernel {math}`K` from ab initio calculations
- **Real-Space Analysis**: Exciton localization and spatial properties

#### Framework Advantages

QuREX provides:
- **Computational Efficiency**: Reduced cost through Wannierization
- **Physical Insight**: Real-space excitonic properties
- **Flexibility**: Adaptable to different materials
- **Accuracy**: Maintains first-principles accuracy 

## Gauge Issues in Nonlinear Optical Responses

### Long-Wavelength Approximation

In nonlinear optical responses, the **long-wavelength limit** assumes spatially uniform electromagnetic fields (dipole approximation). This requires careful gauge choice considerations{cite}`ventura2017gauge`.

### Gauge Representations

#### Velocity Gauge

Electric field via vector potential:
```{math}
:label: eq:vel-gauge
\mathbf{E}(t) = -\frac{\partial \mathbf{A}(t)}{\partial t}
```

**Advantages:** Preserves crystal symmetry, gauge invariant, physical interpretation  
**Disadvantages:** Numerical divergences, computational complexity, slower convergence

#### Length Gauge  

Interaction via position operator:
```{math}
:label: eq:length-gauge
V(\mathbf{r}) = e\mathbf{E}(t) \cdot \mathbf{r}
```

**Advantages:** Numerical stability, computational efficiency, faster convergence  
**Disadvantages:** Broken translational symmetry, surface effects, gauge dependence

### Gauge Transformation

#### Gauge Transformation

Gauges are related by unitary transformation:
```{math}
:label: eq:gauge-transform
\mathcal{U}(t) = \exp \left[i \frac{e}{\hbar} \int d^3\mathbf{r} \, \mathbf{A}(t) \cdot \mathbf{r} \, \rho(\mathbf{r})\right]
```

Hamiltonian transformation:
```{math}
:label: eq:hamiltonian-transform
H_E(t) = \mathcal{U}(t) H_A(t) \mathcal{U}^{\dagger}(t) + i\hbar \frac{d\mathcal{U}(t)}{dt} \mathcal{U}^{\dagger}(t)
```

Observable transformation:
```{math}
:label: eq:observable-transform
O_E(t) = \mathcal{U}(t) O_A(t) \mathcal{U}^{\dagger}(t)
```

### Practical Considerations

Most calculations use **length gauge** for numerical advantages while ensuring:
1. **Gauge Invariance**: Physical results independent of gauge choice
2. **Proper Limits**: Correct behavior in appropriate limits  
3. **Symmetry Restoration**: Careful treatment of broken symmetries

Validation requires comparing gauges, checking invariance, and validating against experiments.

### Applications in QuREX

Gauge considerations affect:
- **Optical Matrix Elements**: Transition dipole moments
- **Nonlinear Responses**: Higher-order optical processes
- **Real-Space Analysis**: Spatial distribution of transitions
- **Symmetry Analysis**: Group theory applications

Gauge choice affects both implementation and physical interpretation in excitonic calculations.

# References

```{bibliography}

