# ExcitonGroupTheory API Reference

## Class: ExcitonGroupTheory

### Overview

The `ExcitonGroupTheory` class provides comprehensive group theory analysis of exciton states in crystalline materials. It determines the irreducible representations of exciton states under the little group of the exciton momentum. **This class now inherits from `BaseOpticalProperties` for improved code organization and consistency.**

### Class Definition

```python
class ExcitonGroupTheory(BaseOpticalProperties):
    """
    Group theory analysis of exciton states.
    
    This class analyzes the irreducible representations of exciton states under the 
    little group of the exciton momentum, providing insight into the symmetry
    properties of excitonic states.
    
    Inherits from BaseOpticalProperties for common database handling and utilities.
    """
```

### Constructor

```python
def __init__(self, path=None, save='SAVE', lelph_db=None, latdb=None, wfdb=None, 
             bands_range=None, BSE_dir='bse', LELPH_dir='lelph', 
             read_symm_from_ns_db_file=False, save_files=True):
```

#### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `path` | `str` or `None` | `None` | Path to the calculation directory. If `None`, uses current working directory. |
| `save` | `str` | `'SAVE'` | Name of the folder containing Yambo save files. |
| `lelph_db` | `LetzElphElectronPhononDB` or `None` | `None` | Pre-loaded electron-phonon database object. If `None`, reads from file. |
| `latdb` | `YamboLatticeDB` or `None` | `None` | Pre-loaded lattice database object. If `None`, reads from file. |
| `wfdb` | `YamboWFDB` or `None` | `None` | Pre-loaded wavefunction database object. If `None`, reads from file. |
| `bands_range` | `list` or `None` | `None` | Range of bands for analysis `[min_band, max_band]`. If `None`, uses all bands. |
| `BSE_dir` | `str` | `'bse'` | Directory containing BSE calculation files. |
| `LELPH_dir` | `str` | `'lelph'` | Directory containing electron-phonon matrix elements. |
| `read_symm_from_ns_db_file` | `bool` | `False` | If `True`, reads symmetries from `ns.db1`; if `False`, from `ndb.elph`. |
| `save_files` | `bool` | `True` | Whether to save intermediate results to files. |

#### Raises

| Exception | Condition |
|-----------|-----------|
| `IOError` | When required database files cannot be read |
| `FileNotFoundError` | When database files are missing |

#### Example

```python
# Basic initialization
egt = ExcitonGroupTheory(
    path='/path/to/calculation',
    save='SAVE',
    BSE_dir='bse',
    LELPH_dir='lelph'
)

# Advanced initialization with custom parameters
egt = ExcitonGroupTheory(
    path='.',
    save='SAVE',
    bands_range=[1, 20],
    BSE_dir='GW_BSE',
    LELPH_dir='elph_data',
    read_symm_from_ns_db_file=True
)
```

### Attributes

#### Core Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `path` | `str` | Calculation directory path |
| `SAVE_dir` | `str` | Full path to SAVE directory |
| `BSE_dir` | `str` | Full path to BSE directory |
| `LELPH_dir` | `str` | Full path to electron-phonon directory |

#### Database Objects

| Attribute | Type | Description |
|-----------|------|-------------|
| `ydb` | `YamboLatticeDB` | Lattice database object |
| `wfdb` | `YamboWFDB` | Wavefunction database object |
| `lelph_db` | `LetzElphElectronPhononDB` | Electron-phonon database object |

#### Structural Information

| Attribute | Type | Description |
|-----------|------|-------------|
| `lat_vecs` | `numpy.ndarray` | Lattice vectors (3×3 array) |
| `blat_vecs` | `numpy.ndarray` | Reciprocal lattice vectors (3×3 array) |
| `nibz` | `int` | Number of IBZ k-points |
| `nkpoints` | `int` | Total number of k-points |
| `nspin` | `int` | Number of spin channels |
| `nspinor` | `int` | Number of spinor components |
| `nbands` | `int` | Number of bands |

#### Symmetry Information

| Attribute | Type | Description |
|-----------|------|-------------|
| `symm_mats` | `numpy.ndarray` | Symmetry matrices in Cartesian coordinates |
| `sym_red` | `numpy.ndarray` | Symmetry matrices in reduced coordinates |
| `time_rev` | `bool` | Time reversal symmetry flag |
| `frac_trans` | `numpy.ndarray` | Fractional translations for each symmetry |

#### K-point Information

| Attribute | Type | Description |
|-----------|------|-------------|
| `kpts` | `numpy.ndarray` | All k-points |
| `kpts_ibz` | `numpy.ndarray` | IBZ k-points |
| `kpt_tree` | `KDTree` | K-point search tree |
| `qpts` | `numpy.ndarray` | Q-points from electron-phonon calculation |

#### Matrix Elements

| Attribute | Type | Description |
|-----------|------|-------------|
| `Dmats` | `numpy.ndarray` | D-matrices for wavefunction rotation |
| `elph_bnds_range` | `list` | Band range from electron-phonon calculation |

### Methods

#### compute()

```python
def compute(self):
```

Main computation method - placeholder for group theory analysis.

**Returns:**
- `dict`: Empty dictionary (placeholder implementation)

**Note:** This is a placeholder method. Use `analyze_exciton_symmetry()` for actual analysis.

**Example:**
```python
# Main compute method (placeholder)
results = egt.compute()
```

#### read()

```python
def read(self, lelph_db=None, latdb=None, wfdb=None, bands_range=None):
```

Read and initialize all required database objects. **Now uses base class functionality for common operations.**

**Parameters:**
- `lelph_db` (`LetzElphElectronPhononDB`, optional): Pre-loaded electron-phonon database
- `latdb` (`YamboLatticeDB`, optional): Pre-loaded lattice database  
- `wfdb` (`YamboWFDB`, optional): Pre-loaded wavefunction database
- `bands_range` (`list`, optional): Band range for analysis

**Raises:**
- `IOError`: If database files cannot be read

**Example:**
```python
# Re-read databases with different parameters
egt.read(bands_range=[1, 15])
```

#### read_excdb_single()

```python
def read_excdb_single(self, BSE_dir, iQ, nstates):
```

Read Yambo exciton database for a specific Q-point. **Renamed from `read_excdb()` for clarity.**

**Parameters:**
- `BSE_dir` (`str`): Directory containing BSE calculation data
- `iQ` (`int`): Q-point index (1-based indexing as in Yambo)
- `nstates` (`int`): Number of exciton states to read

**Returns:**
- `tuple`: (bands_range, BS_eigs, BS_wfcs) for the specific Q-point

**Raises:**
- `IOError`: If BSE database file cannot be read

**Example:**
```python
bands, energies, wavefunctions = egt.read_excdb_single('bse', iQ=1, nstates=10)
```

#### analyze_exciton_symmetry()

```python
def analyze_exciton_symmetry(self, iQ, nstates, degen_thres=0.001):
```

Perform comprehensive group theory analysis for exciton states.

**Parameters:**
- `iQ` (`int`): Q-point index (1-based indexing as in Yambo)
- `nstates` (`int`): Number of exciton states to analyze
- `degen_thres` (`float`, optional): Degeneracy threshold in eV (default: 0.001)

**Returns:**
- `results` (`dict`): Analysis results dictionary

**Result Dictionary Structure:**
```python
{
    'q_point': numpy.ndarray,           # Q-point coordinates
    'little_group': numpy.ndarray,      # Little group operation indices
    'point_group_label': str,           # Point group symbol
    'unique_energies': numpy.ndarray,   # Unique energy levels (eV)
    'degeneracies': numpy.ndarray,      # Degeneracy of each level
    'irrep_decomposition': list,        # Irrep decomposition strings
    'exciton_energies': numpy.ndarray,  # All exciton energies (eV)
    'classes': list,                    # Symmetry classes
    'class_dict': dict                  # Class to operation mapping
}
```

**Mathematical Details:**

The method performs the following steps:

1. **Little Group Determination**: Finds all symmetry operations {math}`g` such that {math}`g\mathbf{Q} = \mathbf{Q} + \mathbf{G}`
2. **Wavefunction Rotation**: Applies symmetry operations to exciton wavefunctions
3. **Representation Matrix**: Computes {math}`D_{\mu\lambda}^{(g)} = \langle\psi_{\mu}|g|\psi_{\lambda}\rangle`
4. **Character Calculation**: Computes traces {math}`\chi^{(g)} = \text{Tr}[D^{(g)}]`
5. **Irrep Decomposition**: Uses reduction formula to decompose representations

**Example:**
```python
results = egt.analyze_exciton_symmetry(
    iQ=1,
    nstates=10,
    degen_thres=0.001
)

print(f"Point group: {results['point_group_label']}")
for i, (E, deg, irrep) in enumerate(zip(
    results['unique_energies'],
    results['degeneracies'],
    results['irrep_decomposition'])):
    print(f"Level {i+1}: {E:.4f} eV (deg={deg}) -> {irrep}")
```

#### save_analysis_results()

```python
def save_analysis_results(self, results, filename=None):
```

Save group theory analysis results to a text file.

**Parameters:**
- `results` (`dict`): Results dictionary from `analyze_exciton_symmetry()`
- `filename` (`str`, optional): Output filename. If `None`, generates automatic name.

**Output File Format:**
```
Exciton Group Theory Analysis
========================================
Q-point: [0.0000, 0.0000, 0.0000]
Little group: C2v
Little group symmetries: [1 2 3 4]

Classes:
  E
  C2
  σv
  σv'

Exciton representations:
Energy (eV)    Degeneracy    Representation
--------------------------------------------------
  2.1500            1         A1
  2.3200            2         E
  2.4100            1         B1
```

**Example:**
```python
# Save with automatic filename
egt.save_analysis_results(results)

# Save with custom filename
egt.save_analysis_results(results, 'my_analysis.txt')
```

## Point Group Operations Module

### Overview

The point group operations module now uses the **spgrep library** for state-of-the-art crystallographic point group analysis. This modern implementation provides:

- **Automatic point group identification** using International Tables standards
- **Comprehensive character tables** from the spgrep database  
- **Robust irreducible representation matrices**
- **Automatic fallback** to original implementation if spgrep unavailable

### spgrep Integration Benefits

- **Enhanced Accuracy**: Uses crystallographic standards from International Tables
- **Comprehensive Coverage**: Supports all crystallographic point groups
- **Maintained Compatibility**: Automatic fallback ensures existing workflows continue
- **Performance**: Optimized algorithms from the spgrep library

### Functions

#### get_pg_info()

```python
def get_pg_info(symm_mats):
```

Identify point group and return character table information. **Now follows the original algorithm exactly.**

**Parameters:**
- `symm_mats` (`numpy.ndarray`): Array of symmetry matrices

**Returns:**
- `pg_label` (`str`): Point group label
- `classes` (`list`): List of symmetry classes
- `class_dict` (`dict`): Mapping of class indices to operations
- `char_tab` (`numpy.ndarray`): Character table
- `irreps` (`list`): List of irreducible representation labels

**Algorithm Details:**
1. **Point group identification** using `get_point_grp()`
2. **Symmetry element generation** via `pg_to_symels()`
3. **Matrix transformation** to standard orientation
4. **Character table lookup** from comprehensive database
5. **Class mapping** using original classification logic

#### decompose_rep2irrep()

```python
def decompose_rep2irrep(red_rep, char_table, pg_order, class_order, irreps):
```

Decompose reducible representation into irreducible components. **Uses the exact reduction formula from the original implementation.**

**Parameters:**
- `red_rep` (`numpy.ndarray`): Characters of reducible representation
- `char_table` (`numpy.ndarray`): Character table
- `pg_order` (`int`): Order of point group
- `class_order` (`numpy.ndarray`): Order of each class
- `irreps` (`list`): Irreducible representation labels

**Returns:**
- `decomposition` (`str`): String representation of decomposition

**Mathematical Formula:**
```python
irrep_coeff = np.einsum('j,j,rj->r', class_order, red_rep, char_table, optimize=True) / pg_order
```

#### Core Algorithm Functions

| Function | Description |
|----------|-------------|
| `find_symm_axis(sym_mats)` | Find symmetry axes and n-fold values | 
| `get_point_grp(symm_mats)` | Classify point group using flowchart | 
| `find_axis_angle(Rmat)` | Extract rotation axis and angle | 
| `fix_axis_angle_gauge(axis, nfold)` | Fix axis gauge convention |

#### Utility Functions

| Function | Description | Returns | Improvements |
|----------|-------------|---------|--------------|
| `normalize(a)` | Normalize a vector | `numpy.ndarray` | **Optimized for zero vectors** |
| `rotation_matrix(axis, theta)` | Create rotation matrix | `numpy.ndarray` | **Rodrigues' formula implementation** |
| `reflection_matrix(axis)` | Create reflection matrix | `numpy.ndarray` | **Optimized matrix construction** |
| `inversion_matrix()` | Create inversion matrix | `numpy.ndarray` | **Simple -I implementation** |
| `transform_matrix(old, new)` | Find transformation between groups | `numpy.ndarray` | **Enhanced axis alignment** |

#### Character Table Generation

```python
def pg_to_chartab(PG):
```

Generate character tables for common point groups. **Expanded database with accurate character tables.**

**Supported Point Groups:**
- **C groups**: C1, C2, Cs, Ci, C2v, C3v, C4v, C6v, Cnh series
- **D groups**: D2, D3, D4, D6, D2h, D3h, D4h, D6h, Dnd series  
- **Cubic groups**: T, Td, Th, O, Oh
- **Special cases**: S4, S6, and other improper rotation groups

#### Symmetry Element Generation

```python
def pg_to_symels(PG):
```

Generate symmetry elements for point groups. **Follows original MolSym structure.**

**Features:**
- **Automatic element generation** for all supported point groups
- **Proper matrix representations** for all symmetry operations
- **Consistent axis conventions** matching crystallographic standards

## Error Handling

### Common Exceptions

| Exception | Cause | Solution |
|-----------|-------|----------|
| `IOError` | Missing database files | Check file paths and ensure calculations completed |
| `FileNotFoundError` | Incorrect directory structure | Verify directory names and file locations |
| `ValueError` | Invalid parameters | Check parameter ranges and types |
| `MemoryError` | Large system size | Reduce number of states or bands |
| `ImportError` | Missing dependencies | Install required packages (netCDF4, scipy, etc.) |

### Debugging Tips

1. **Check File Existence:**
   ```python
   import os
   required_files = ['SAVE/ns.db1', 'SAVE/ns.wf', 'bse/ndb.BS_diago_Q1']
   for f in required_files:
       print(f"{f}: {'✓' if os.path.exists(f) else '✗'}")
   ```

2. **Verify Database Integrity:**
   ```python
   try:
       from netCDF4 import Dataset
       with Dataset('lelph/ndb.elph', 'r') as f:
           print("Available variables:", list(f.variables.keys()))
   except Exception as e:
       print(f"Database error: {e}")
   ```

3. **Memory Usage Monitoring:**
   ```python
   import psutil
   process = psutil.Process()
   print(f"Memory usage: {process.memory_info().rss / 1024**2:.1f} MB")
   ```

## Performance Considerations

### Optimization Strategies

1. **Reduce System Size:**
   - Limit `nstates` to essential states
   - Use smaller `bands_range`
   - Analyze fewer Q-points

2. **Memory Management:**
   - Process Q-points sequentially
   - Clear intermediate arrays
   - Use appropriate data types

3. **Computational Efficiency:**
   - Pre-load database objects
   - Use optimized linear algebra libraries
   - Enable parallel processing where available

### Scaling Guidelines

| System Size | Recommended Settings | Expected Performance |
|-------------|---------------------|---------------------|
| Small (< 100 atoms) | `nstates=20`, `bands_range=[1,30]` | < 1 minute |
| Medium (100-500 atoms) | `nstates=10`, `bands_range=[1,20]` | 1-10 minutes |
| Large (> 500 atoms) | `nstates=5`, `bands_range=[1,15]` | 10+ minutes |


### Compatibility Notes

- Compatible with Yambo 5.x database formats
- Requires LetzElPhC for D-matrix generation
- Works with both serial and parallel Yambo calculations
- Supports both spin-polarized and non-spin-polarized systems