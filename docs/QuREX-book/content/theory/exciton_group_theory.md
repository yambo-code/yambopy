# Exciton Group Theory Analysis

## Introduction

The `ExcitonGroupTheory` class in Yambopy provides a comprehensive framework for analyzing the symmetry properties of exciton states using group theory. This analysis is crucial for understanding the optical selection rules, degeneracies, and symmetry-allowed transitions in excitonic systems.

**Key Features:**
- **Universal symmetry classification** using spglib for all 230 space groups
- **General crystallographic analysis** supporting all 7 crystal systems
- **Non-symmorphic operations** including screw rotations and glide reflections
- **Automatic space group detection** with International Tables notation
- **Comprehensive operation classification** with proper crystallographic symbols
- **Publication-ready output** with standardized mathematical notation

**Supported Crystal Systems:**
- ✅ **Triclinic** (P1, P-1) - Space groups 1-2
- ✅ **Monoclinic** (P2, P2/m, C2/m) - Space groups 3-15
- ✅ **Orthorhombic** (Pmmm, Cmcm, Fddd) - Space groups 16-74
- ✅ **Tetragonal** (P4, P4/mmm, I4/mcm) - Space groups 75-142
- ✅ **Trigonal** (P3, R3m, P3m1) - Space groups 143-167
- ✅ **Hexagonal** (P6, P6/mmm, P6₃/mmc) - Space groups 168-194 - **Validated with hBN**
- ✅ **Cubic** (Pm3m, Fd3m, Im3m) - Space groups 195-230

**Operation Types Supported:**
- **Symmorphic**: Identity (E), rotations (Cₙ), reflections (σ), inversion (i), rotoinversions (Sₙ)
- **Non-symmorphic**: Screw rotations (2₁, 3₁, 6₁), glide reflections (a, b, c, n, d)

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

## General Symmetry Classification

### Spglib Integration

The implementation leverages **spglib** for universal crystallographic analysis:

1. **Automatic Space Group Detection**: Identifies space group from crystal structure
2. **Operation Matching**: Maps Yambo symmetry matrices to spglib operations
3. **Translation Vector Analysis**: Handles non-symmorphic operations properly
4. **International Tables Compliance**: Uses standard crystallographic notation

### Classification Algorithm

The general classification method `classify_symmetry_operations()` works as follows:

```python
def classify_symmetry_operations(self):
    """
    Classify symmetry operations using spglib for general space group support.
    Works for all 230 space groups with comprehensive operation analysis.
    """
    # 1. Get crystal structure for spglib
    cell = (lattice, positions, numbers)
    
    # 2. Obtain spglib symmetry information
    symmetry = spglib.get_symmetry(cell)
    dataset = spglib.get_symmetry_dataset(cell)
    
    # 3. Match Yambo operations with spglib operations
    # 4. Classify each operation by type and properties
    # 5. Return comprehensive analysis with crystal system info
```

### Operation Classification Types

The method classifies operations into these categories:

- **`identity`**: Identity operation (E)
- **`rotation`**: Proper rotations (C₂, C₃, C₄, C₆)
- **`reflection`**: Mirror planes (σₕ, σᵥ, σₐ)
- **`inversion`**: Inversion center (i)
- **`rotoinversion`**: Improper rotations (S₃, S₄, S₆)
- **`screw`**: Screw rotations (2₁, 3₁, 4₁, 6₁, etc.)
- **`glide`**: Glide reflections (a, b, c, n, d)

## Implementation Details

### Class Structure

The `ExcitonGroupTheory` class implements the following key components:

1. **Database Reading**: Reads Yambo databases including lattice, wavefunction, BSE, and electron-phonon data
2. **Universal Symmetry Analysis**: Uses spglib for general space group identification and operation classification
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

### General Symmetry Classification

```python
# NEW: Universal symmetry operation classification
operations = egt.classify_symmetry_operations()
summary = operations.get('_summary', {})

print(f"Space Group: {summary.get('space_group')} (#{summary.get('space_group_number')})")
print(f"Point Group: {summary.get('point_group')}")
print(f"Crystal System: {summary.get('crystal_system')}")

# Show operation breakdown
operation_types = ['identity', 'rotation', 'reflection', 'inversion', 
                  'rotoinversion', 'screw', 'glide']

for op_type in operation_types:
    op_list = operations.get(op_type, [])
    if op_list:
        print(f"{op_type.title()}: {len(op_list)} operations")
        # Show first operation as example
        if len(op_list[0]) >= 4:
            idx, mat, desc, symbol, spglib_info = op_list[0]
            print(f"  Example: {desc} ({symbol})")

# Display comprehensive analysis
egt.display_symmetry_operations()
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

# NEW: Crystal system specific analysis
crystal_system = summary.get('crystal_system', '').lower()
if crystal_system == 'hexagonal':
    print("Hexagonal system detected - analyzing D6h operations")
elif crystal_system == 'cubic':
    print("Cubic system detected - analyzing high-symmetry operations")
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