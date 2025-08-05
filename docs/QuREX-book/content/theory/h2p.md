# Two-Particle Hamiltonian (H2P)

## Introduction

The **two-particle Hamiltonian** {math}`H_{2p}` forms the core of excitonic calculations within the QuREX framework. This Hamiltonian describes the coupled electron-hole system and incorporates both single-particle energies and many-body interactions that give rise to excitonic effects.

## Theoretical Framework

### General Form

The two-particle Hamiltonian in the electron-hole basis can be written as:

```{math}
:label: eq:h2p-general
H_{2p} = H_0 + V_{eh}
```

where:
- {math}`H_0` represents the non-interacting electron-hole energies
- {math}`V_{eh}` describes the electron-hole interactions (Coulomb kernel)

### Matrix Representation

In the basis of electron-hole pairs {math}`|vc\mathbf{k}\rangle`, the Hamiltonian matrix elements are:

```{math}
:label: eq:h2p-matrix
H^{2p}_{vc\mathbf{k},v'c'\mathbf{k}'} = (\epsilon_{c\mathbf{k}} - \epsilon_{v\mathbf{k}}) \delta_{cc'} \delta_{vv'} \delta_{\mathbf{k}\mathbf{k}'} + K_{vc\mathbf{k},v'c'\mathbf{k}'}
```

where {math}`K` is the electron-hole interaction kernel from the [BSE equation](bse_equation).

## Implementation in QuREX

### Construction Methods

The BSE Hamiltonian in `wann_H2P.py` can be constructed using three approaches:

1. **Model Coulomb Potential**: Using analytical or fitted [Coulomb potentials](coulomb_potential)
2. **Yambo Kernel**: Direct import from `YamboBSEKernelDB` 
3. **Reconstruction**: From `YamboExcitonDB` eigenvalues and eigenvectors

### Wannier Representation

In the Wannier basis, the Hamiltonian becomes:

```{math}
:label: eq:h2p-wannier
H^{2p}_{WW'} = \langle W | H_{2p} | W' \rangle
```

where {math}`|W\rangle` represents Wannier exciton states as described in [Wannier Excitons](wannier_exciton).

## Computational Aspects

### Matrix Diagonalization

The eigenvalue problem:

```{math}
:label: eq:h2p-eigen
H^{2p} |\psi_\lambda\rangle = E_\lambda |\psi_\lambda\rangle
```

yields exciton energies {math}`E_\lambda` and wavefunctions {math}`|\psi_\lambda\rangle`.

### Convergence Considerations

Key parameters affecting convergence:
- **k-point sampling**: Density of the Brillouin zone mesh
- **Band range**: Number of valence and conduction bands included
- **Cutoff energies**: For Coulomb interaction truncation

## Applications

The two-particle Hamiltonian enables calculation of:
- Exciton binding energies
- Optical absorption spectra  
- Exciton spatial distributions
- Temperature-dependent properties
- Phonon coupling matrix elements

## Connection to Other Methods

The H2P formalism connects to:
- **GW calculations**: Through quasiparticle energies {math}`\epsilon_{c/v\mathbf{k}}`
- **BSE theory**: Via the interaction kernel {math}`K`
- **Model Hamiltonians**: Through effective parameters
- **Real-space analysis**: Via Wannier function transformations
