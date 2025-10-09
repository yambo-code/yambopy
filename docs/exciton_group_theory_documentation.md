# Exciton Group Theory Analysis

## Theory

The symmetry analysis of exciton states provides crucial insight into their optical selection rules and physical properties. The `ExcitonGroupTheory` class implements a comprehensive framework for analyzing the irreducible representations of exciton states under crystallographic point groups.

### Mathematical Framework

For an exciton state at momentum **Q**, the little group G(**Q**) consists of all symmetry operations **R** that leave **Q** invariant:

```
G(Q) = {R ∈ G | RQ = Q + G}
```

where **G** is a reciprocal lattice vector. The exciton wavefunction transforms under symmetry operations as:

```
Ψ_λ^R(k) = e^{i2π Q·τ_R} D_R(k) Ψ_λ(R^(-1)k)
```

where:
- `D_R(k)` is the rotation matrix in the Bloch basis
- `τ_R` is the fractional translation associated with operation **R**
- `λ` labels the exciton state

### Implementation Architecture

The implementation uses a dual-symmetry approach:

1. **Yambo symmetries**: Used for wavefunction rotations and representation matrix construction
2. **spglib/spgrep symmetries**: Used for crystallographic point group identification and irrep analysis

This separation ensures compatibility with Yambo's internal symmetry handling while leveraging the comprehensive crystallographic databases in spglib/spgrep.

### Point Group Identification

The algorithm automatically identifies the crystallographic point group through:

1. **Crystal structure analysis**: Reads atomic positions and lattice vectors from Yambo SAVE database
2. **Space group determination**: Uses spglib to identify the space group from crystal structure
3. **Point group extraction**: Extracts point group operations from space group symmetries
4. **spgrep classification**: Uses spgrep's internal database for point group identification

### Irreducible Representation Analysis

The irrep decomposition follows standard group theory procedures:

1. **Character calculation**: Computes traces of representation matrices for each symmetry class
2. **Projection operators**: Uses projection formulas to decompose reducible representations
3. **Degeneracy analysis**: Groups states by energy degeneracy within specified threshold
4. **Label assignment**: Assigns systematic Γ_n labels to irreducible representations

## API Reference

### ExcitonGroupTheory Class

```python
class ExcitonGroupTheory(BaseOpticalProperties):
    """
    Group theory analysis of exciton states under crystallographic symmetries.
    
    Automatically reads crystal structure from Yambo databases and performs
    comprehensive symmetry analysis using spglib/spgrep libraries.
    """
```

#### Constructor

```python
def __init__(self, path=None, save='SAVE', lelph_db=None, latdb=None, wfdb=None, 
             bands_range=None, BSE_dir='bse', LELPH_dir='lelph', 
             read_symm_from_ns_db_file=True):
```

**Parameters:**
- `path` (str): Calculation directory path
- `save` (str): SAVE directory name (default: 'SAVE')
- `lelph_db` (LetzElphElectronPhononDB): Pre-loaded electron-phonon database
- `latdb` (YamboLatticeDB): Pre-loaded lattice database  
- `wfdb` (YamboWFDB): Pre-loaded wavefunction database
- `bands_range` (list): Band range for analysis
- `BSE_dir` (str): BSE calculation directory (default: 'bse')
- `LELPH_dir` (str): Electron-phonon directory (default: 'lelph')
- `read_symm_from_ns_db_file` (bool): Read symmetries from ns.db1 (default: True)

#### Main Analysis Method

```python
def analyze_exciton_symmetry(self, iQ, nstates, degen_thres=0.001):
    """
    Perform comprehensive group theory analysis for exciton states.
    
    Parameters
    ----------
    iQ : int
        Q-point index (1-based, following Yambo convention)
    nstates : int
        Number of exciton states to analyze
    degen_thres : float
        Energy degeneracy threshold in eV (default: 0.001)
        
    Returns
    -------
    results : dict
        Analysis results containing:
        - 'little_group': Little group symmetry operations
        - 'point_group_label': Crystallographic point group symbol
        - 'unique_energies': Unique energy levels
        - 'degeneracies': Degeneracy of each level
        - 'irrep_decomposition': Irreducible representation labels
        - 'character_table': Character table for the point group
        - 'symmetry_classes': Conjugacy classes of symmetry operations
    """
```

#### Crystal Structure Reading

```python
def _get_crystal_structure(self):
    """
    Extract crystal structure from Yambo SAVE database.
    
    Automatically reads:
    - Lattice vectors from lattice database
    - Atomic positions (fractional coordinates)
    - Atomic numbers for all atoms in unit cell
    
    Returns
    -------
    tuple
        (lattice_vectors, atomic_positions, atomic_numbers)
    """
```

### Integration with spglib/spgrep

The implementation leverages two key external libraries:

#### spglib Integration
- **Space group identification**: `spglib.get_spacegroup()`
- **Symmetry operations**: `spglib.get_symmetry()`
- **Standardization**: Ensures crystallographic conventions

#### spgrep Integration  
- **Point group classification**: `spgrep.pointgroup.get_pointgroup()`
- **Irrep generation**: `spgrep.get_crystallographic_pointgroup_irreps_from_symmetry()`
- **Character tables**: Automatic generation from symmetry operations

### Example Usage

```python
from yambopy.optical_properties.exciton_group_theory import ExcitonGroupTheory

# Initialize analysis
egt = ExcitonGroupTheory(
    path='/path/to/calculation',
    BSE_dir='bse',
    bands_range=[6, 10]
)

# Analyze exciton symmetry at Γ point
results = egt.analyze_exciton_symmetry(
    iQ=1,           # Γ point (Q=0)
    nstates=10,     # First 10 exciton states
    degen_thres=0.001  # 1 meV degeneracy threshold
)

# Access results
print(f"Point group: {results['point_group_label']}")
print(f"Little group size: {len(results['little_group'])}")

for i, (E, deg, irrep) in enumerate(zip(
    results['unique_energies'], 
    results['degeneracies'],
    results['irrep_decomposition']
)):
    print(f"State {i+1}: {E:.3f} eV (deg={deg}) → {irrep}")
```

### Output Format

The analysis provides systematic irrep labels using the universal Γ_n notation:

```
Energy (eV)    Degeneracy    Representation
7.364          2             Γ₅
7.891          1             Γ₁  
8.123          2             Γ₆
```

This notation is independent of specific point group conventions and universally understood in the solid-state physics community.

## References

1. **Group Theory**: Tinkham, M. "Group Theory and Quantum Mechanics" (Dover, 2003)
2. **Crystallographic Point Groups**: Bradley, C. J. & Cracknell, A. P. "The Mathematical Theory of Symmetry in Solids" (Oxford, 1972)
3. **spglib**: Togo, A. & Tanaka, I. "Spglib: a software library for crystal symmetry search" arXiv:1808.01590 (2018)
4. **spgrep**: Watanabe, H. et al. "spgrep: On-the-fly generator of space-group irreducible representations" J. Open Source Softw. 8, 5269 (2023)
5. **Exciton Symmetry**: Cho, K. "Optical Response of Nanostructures" (Springer, 2003)