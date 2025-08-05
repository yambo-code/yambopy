(model_ham_intro)=
# Model Hamiltonian Construction

## Introduction to Tight-Binding Parametrization

The construction of effective model Hamiltonians represents a crucial bridge between first-principles calculations and computationally efficient descriptions of electronic and excitonic properties. The **maximally localized Wannier function tight-binding (MLWF-TB)** parametrization{cite}`marzari2012maximally` serves as the workhorse method for this transformation, enabling the extraction of physically meaningful tight-binding parameters from density functional theory (DFT) calculations.

## Theoretical Framework

### From DFT to Tight-Binding

The MLWF-TB approach extracts the electronic tight-binding Hamiltonian ({math}`H`) from converged DFT calculations using established open-source codes such as:

- **Quantum ESPRESSO**{cite}`giannozzi2009quantum,giannozzi2017advanced`: Ground-state DFT calculations
- **Wannier90**{cite}`mostofi2008wannier90`: Wannier function construction and tight-binding extraction

This procedure yields a **real-space representation** of the electronic Hamiltonian:

{math}`H_{nm}(\mathbf{R})`

where:
- {math}`\mathbf{R}` are lattice vectors within a supercell conjugate to the DFT k-mesh
- {math}`n, m` label the electronic band indices (or Wannier function indices)
- The matrix elements represent hopping integrals between Wannier orbitals

### Slater-Koster Interpolation

The real-space Hamiltonian enables **k-space interpolation** via the Slater-Koster scheme{cite}`yates2007spectral`, allowing calculation of the reciprocal space Hamiltonian on arbitrarily fine k-meshes:

```{math}
:label: eq:HR
H_{nm}(\mathbf{k}) = \sum_{\mathbf{R}} e^{i\mathbf{k} \cdot \mathbf{R}} H_{nm}(\mathbf{R})
```

This interpolation provides several key advantages:

1. **Computational Efficiency**: Avoid expensive DFT calculations on fine k-meshes
2. **Physical Insight**: Real-space hopping parameters reveal bonding characteristics
3. **Transferability**: Model parameters applicable across different conditions
4. **Scalability**: Efficient calculations for large systems and complex geometries

### Hamilton Gauge Convention

After diagonalization of the tight-binding Hamiltonian, the resulting eigenvectors are expressed in the **Hamilton gauge**, which provides a natural basis for subsequent many-body calculations and ensures proper transformation properties under symmetry operations.

## Applications in Excitonic Systems

### Two-Particle Hamiltonian Construction

The tight-binding framework naturally extends to excitonic systems through the construction of **two-particle Hamiltonians** that describe electron-hole interactions. The model Hamiltonian approach enables:

- **Efficient BSE Calculations**: Reduced computational cost for large systems
- **Parameter Studies**: Systematic investigation of material properties
- **Real-Space Analysis**: Understanding of exciton localization and binding
- **Heterostructure Modeling**: Treatment of complex multi-layer systems

### Advantages for Complex Systems

Model Hamiltonians prove particularly valuable for:

1. **Large-Scale Systems**: Hundreds to thousands of atoms
2. **Heterostructures**: Multiple materials and interfaces
3. **Defect Studies**: Point defects, grain boundaries, edges
4. **Temperature Effects**: Thermal expansion and phonon coupling
5. **External Fields**: Electric and magnetic field responses

## Implementation Strategy

### Workflow Overview

The typical workflow for model Hamiltonian construction involves:

1. **DFT Ground State**: Converged electronic structure calculation
2. **Wannier Construction**: Maximally localized Wannier functions
3. **Tight-Binding Extraction**: Real-space Hamiltonian matrix elements
4. **Validation**: Comparison with original DFT band structure
5. **Application**: Use in many-body calculations (GW, BSE)

### Quality Control

Ensuring reliable model Hamiltonians requires:

- **Convergence Testing**: k-mesh, energy windows, localization spreads
- **Band Structure Reproduction**: Accurate interpolation of DFT results
- **Wannier Function Analysis**: Proper localization and chemical bonding
- **Symmetry Verification**: Preservation of crystal symmetries

## Physical Interpretation

### Hopping Integrals

The real-space matrix elements {math}`H_{nm}(\mathbf{R})` provide direct physical insight:

- **On-site Terms** ({math}`\mathbf{R} = 0`): Atomic energy levels and local environment
- **Nearest-Neighbor Hopping**: Primary bonding interactions
- **Long-Range Terms**: Extended interactions and screening effects
- **Orbital Character**: s, p, d orbital contributions and hybridization

### Chemical Bonding

The tight-binding parameters reveal:

- **Bond Strengths**: Magnitude of hopping integrals
- **Directional Bonding**: Anisotropy in hopping patterns
- **Orbital Overlap**: Spatial extent of Wannier functions
- **Electronic Correlations**: Deviations from simple tight-binding behavior

## Advanced Topics

### Beyond Standard Tight-Binding

Modern implementations extend the basic framework:

- **Non-Orthogonal Basis**: Overlap matrix inclusion for improved accuracy
- **Spin-Orbit Coupling**: Relativistic effects in heavy-element systems
- **Many-Body Corrections**: GW quasiparticle energy incorporation
- **Dynamic Effects**: Frequency-dependent interactions

### Integration with Many-Body Theory

The model Hamiltonian serves as the foundation for:

- **GW Calculations**: Quasiparticle energy corrections
- **BSE Solutions**: Excitonic effects and optical spectra
- **DMFT Applications**: Strong correlation effects
- **Transport Properties**: Conductivity and mobility calculations

---

*The model Hamiltonian approach represents a powerful methodology that combines the accuracy of first-principles calculations with the efficiency and insight of tight-binding models, enabling comprehensive studies of excitonic materials across multiple length and energy scales.*