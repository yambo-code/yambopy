# Exciton Group Theory Analysis

## Introduction

The `ExcitonGroupTheory` class in Yambopy provides a comprehensive framework for analyzing the symmetry properties of exciton states using group theory. This analysis is crucial for understanding the optical selection rules, degeneracies, and symmetry-allowed transitions in excitonic systems.

## Theoretical Background

### Exciton States and Symmetry

An exciton is a bound state of an electron and a hole, described by the wavefunction:

{math}`|\psi^{\lambda}(\mathbf{Q})\rangle = \sum_{\mathbf{k},v,c} A^{\lambda}_{vc}(\mathbf{k},\mathbf{Q}) |v\mathbf{k-Q}, c(\mathbf{Q})\rangle`

where:
- {math}`\lambda` is the exciton state index
- {math}`\mathbf{Q}` is the exciton center-of-mass momentum
- {math}`A^{\lambda}_{vc}(\mathbf{k},\mathbf{Q})` are the exciton amplitudes
- {math}`|v\mathbf{k-Q}, c(\mathbf{k})\rangle` represents an electron-hole pair state

### Group Theory Analysis

The symmetry properties of exciton states are determined by the little group of the exciton momentum {math}`\mathbf{Q}`. The little group {math}`G_{\mathbf{Q}}` consists of all symmetry operations {math}`g` of the crystal point group that leave {math}`\mathbf{Q}` invariant:

{math}`G_{\mathbf{Q}} = \{g \in G : g\mathbf{Q} = \mathbf{Q} + \mathbf{G}\}`

where {math}`\mathbf{G}` is a reciprocal lattice vector.

### Symmetry Operations on Exciton States

Under a symmetry operation {math}`g`, the exciton wavefunction transforms as:

{math}`g|\psi^{\lambda}(\mathbf{Q})\rangle = \sum_{\mu} D_{\mu\lambda}^{(g)}|\psi_{\mu}(\mathbf{Q})\rangle`

where {math}`D^{(g)}` is the representation matrix of the symmetry operation {math}`g`.

The representation matrix elements are calculated as:

{math}`D_{\mu\lambda}^{(g)} = \langle\psi_{\mu}(\mathbf{Q})|g|\psi_{\lambda}(\mathbf{Q})\rangle`

### Character Analysis

The character of a representation for symmetry operation {math}`g` is:

{math}`\chi^{(g)} = \text{Tr}[D^{(g)}] = \sum_{\lambda} D_{\lambda\lambda}^{(g)}`

For degenerate states with the same energy, the character is computed over the degenerate subspace.

### Irreducible Representation Decomposition

The reducible representation {math}`\Gamma` can be decomposed into irreducible representations using the reduction formula:

{math}`a_i = \frac{1}{|G|} \sum_{g \in G} \chi^{(g)} \chi_i^{(g)*}`

where:
- {math}`a_i` is the multiplicity of irreducible representation {math}`\Gamma_i`
- {math}`|G|` is the order of the group
- {math}`\chi_i^{(g)}` is the character of irreducible representation {math}`\Gamma_i` for operation {math}`g`

## Implementation Details

### Class Structure

The `ExcitonGroupTheory` class implements the following key components:

1. **Database Reading**: Reads Yambo databases including lattice, wavefunction, BSE, and electron-phonon data
2. **Symmetry Analysis**: Determines the little group and applies symmetry operations
3. **Wavefunction Rotation**: Rotates exciton wavefunctions using D-matrices
4. **Character Calculation**: Computes representation characters
5. **Irrep Decomposition**: Decomposes reducible representations

### Mathematical Implementation

#### Wavefunction Rotation

The rotation of exciton wavefunctions under symmetry operations involves:

1. **K-point mapping**: For each k-point {math}`\mathbf{k}`, find {math}`\mathbf{k}' = g\mathbf{k}`
2. **D-matrix application**: Apply the D-matrix to rotate the Bloch functions:
   {math}`\psi_{n\mathbf{k}'}(\mathbf{r}) = \sum_{m} D_{mn}^{(g)}(\mathbf{k}) \psi_{m\mathbf{k}}(\mathbf{r})`
3. **Phase factor**: Include the phase factor from fractional translations:
   {math}`e^{i\mathbf{Q} \cdot \boldsymbol{\tau}_g}`

#### Representation Matrix Calculation

The representation matrix is computed as:

{math}`D_{\mu\lambda}^{(g)} = \sum_{\mathbf{k},v,c} A^{\mu*}_{vc}(\mathbf{k}',\mathbf{Q}) \sum_{v',c'} D_{v'v}^{(g)}(\mathbf{k}-\mathbf{Q}) D_{c'c}^{(g)}(\mathbf{k}+\mathbf{Q}) A^{\lambda}_{v'c'}(\mathbf{k},\mathbf{Q}) e^{i\mathbf{Q} \cdot \boldsymbol{\tau}_g}`

where {math}`\mathbf{k}' = g\mathbf{k}`.

### Degeneracy Analysis

States are considered degenerate if their energy difference is below a threshold:

{math}`|E_{\lambda} - E_{\mu}| < \epsilon_{\text{deg}}`

The analysis groups degenerate states and computes the representation for each degenerate subspace.

## Usage Examples

### Basic Usage

```python
from yambopy.optical_properties import ExcitonGroupTheory

# Initialize the class
egt = ExcitonGroupTheory(
    path='.',
    save='SAVE',
    BSE_dir='bse',
    LELPH_dir='lelph',
    bands_range=[1, 20]
)

# Perform group theory analysis
results = egt.analyze_exciton_symmetry(
    iQ=1,           # Q-point index
    nstates=10,     # Number of states
    degen_thres=0.001  # Degeneracy threshold in eV
)

# Save results
egt.save_analysis_results(results, 'exciton_symmetry.txt')
```

### Advanced Analysis

```python
# Access detailed results
print(f"Point group: {results['point_group_label']}")
print(f"Little group operations: {results['little_group']}")

# Analyze each energy level
for i, (energy, degen, irrep) in enumerate(zip(
    results['unique_energies'],
    results['degeneracies'],
    results['irrep_decomposition'])):
    print(f"Level {i+1}: {energy:.4f} eV (deg={degen}) -> {irrep}")
```

## Required Input Files

The analysis requires the following Yambo database files:

1. **`SAVE/ns.db1`**: Lattice and symmetry information
2. **`SAVE/ns.wf`**: Wavefunction database
3. **`BSE_dir/ndb.BS_diago_Q{iQ}`**: BSE eigenvalues and eigenvectors
4. **`LELPH_dir/ndb.elph`**: Electron-phonon database (for symmetries)
5. **`LELPH_dir/ndb.Dmats`**: D-matrices for wavefunction rotation

## Output and Results

The analysis provides:

- **Point group identification**: Determines the little group of the Q-point
- **Energy level classification**: Groups states by energy with degeneracy analysis
- **Irreducible representations**: Decomposes each energy level into irreps
- **Symmetry characters**: Provides character tables and representation matrices

### Result Dictionary Structure

```python
results = {
    'q_point': array,                    # Q-point coordinates
    'little_group': array,               # Little group operations
    'point_group_label': str,            # Point group symbol
    'unique_energies': array,            # Unique energy levels
    'degeneracies': array,               # Degeneracy of each level
    'irrep_decomposition': list,         # Irrep decomposition strings
    'exciton_energies': array,           # All exciton energies
    'classes': list,                     # Symmetry classes
    'class_dict': dict                   # Class to operation mapping
}
```

## Applications

### Optical Selection Rules

The irreducible representation analysis helps determine:
- Which transitions are optically allowed
- Polarization dependence of optical transitions
- Dark vs. bright exciton classification

### Symmetry-Breaking Effects

The analysis can reveal:
- Effects of strain or external fields
- Symmetry lowering in heterostructures
- Interface-induced symmetry breaking

### Material Design

Group theory analysis aids in:
- Predicting optical properties of new materials
- Understanding exciton fine structure
- Designing materials with specific symmetries

## Limitations and Considerations

1. **Approximations**: The analysis assumes the validity of the BSE within the chosen approximations
2. **Numerical precision**: Results depend on the degeneracy threshold and numerical accuracy
3. **Database quality**: Requires high-quality Yambo calculations with sufficient k-point sampling
4. **Point group coverage**: The current implementation covers common point groups but may need extension for exotic symmetries

## References

1. Tinkham, M. "Group Theory and Quantum Mechanics" (1964)
2. Dresselhaus, M. S., Dresselhaus, G., & Jorio, A. "Group Theory: Application to the Physics of Condensed Matter" (2007)
3. Rohlfing, M. & Louie, S. G. "Electron-hole excitations and optical spectra from first principles" Phys. Rev. B 62, 4927 (2000)
4. Onida, G., Reining, L. & Rubio, A. "Electronic excitations: density-functional versus many-body Green's-function approaches" Rev. Mod. Phys. 74, 601 (2002)