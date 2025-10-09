# Optical Properties API Reference

## Overview

The optical properties module provides a unified framework for computing various optical and excitonic properties in crystalline materials. The module is organized around a base class `BaseOpticalProperties` with specialized derived classes for specific calculations.

## Class Hierarchy

```
BaseOpticalProperties
├── ExcitonGroupTheory
├── ExcitonDipole  
├── ExcitonPhonon
└── Luminescence
```

## BaseOpticalProperties

### Overview

The `BaseOpticalProperties` class provides common functionality for all optical properties calculations, including database handling, file I/O, and utility methods.

### Constructor

```python
def __init__(self, path=None, save='SAVE', latdb=None, wfdb=None, 
             bands_range=None, BSE_dir='bse', save_files=True):
```

#### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `path` | `str` or `None` | `None` | Path to calculation directory |
| `save` | `str` | `'SAVE'` | SAVE directory name |
| `latdb` | `YamboLatticeDB` or `None` | `None` | Pre-loaded lattice database |
| `wfdb` | `YamboWFDB` or `None` | `None` | Pre-loaded wavefunction database |
| `bands_range` | `list` or `None` | `None` | Range of bands to load |
| `BSE_dir` | `str` | `'bse'` | BSE directory name |
| `save_files` | `bool` | `True` | Whether to save intermediate files |

### Common Methods

#### read_common_databases()

```python
def read_common_databases(self, latdb=None, wfdb=None, bands_range=None):
```

Read lattice and wavefunction databases common to all optical calculations.

#### read_excdb()

```python
def read_excdb(self, BSE_dir):
```

Read exciton databases for all Q-points.

#### compute()

```python
def compute(self):
```

Abstract method to be implemented by derived classes.

## ExcitonDipole

### Overview

Computes exciton-photon coupling matrix elements for optical absorption and emission processes.

### Constructor

```python
def __init__(self, path=None, save='SAVE', latdb=None, wfdb=None, 
             ydipdb=None, bands_range=None, BSE_dir='bse', 
             DIP_dir='gw', save_files=True):
```

### Key Methods

#### compute_Exdipole()

```python
def compute_Exdipole(self):
```

Compute exciton-dipole matrix elements.

**Returns:**
- `numpy.ndarray`: Exciton-dipole matrix elements

**Example:**
```python
ex_dip = ExcitonDipole(path='./calculation')
dipole_matrix = ex_dip.compute()
```

## ExcitonPhonon

### Overview

Computes exciton-phonon coupling matrix elements for studying vibronic effects and luminescence.

### Constructor

```python
def __init__(self, path=None, save='SAVE', lelph_db=None, latdb=None, wfdb=None, 
             ydipdb=None, bands_range=None, BSE_dir='bse', LELPH_dir='lelph', 
             DIP_dir='gw', save_files=True):
```

### Key Methods

#### compute_ExPhonon()

```python
def compute_ExPhonon(self):
```

Compute exciton-phonon coupling matrix elements.

**Returns:**
- `numpy.ndarray`: Placeholder return (full implementation needed)

**Example:**
```python
ex_ph = ExcitonPhonon(path='./calculation', LELPH_dir='lelph')
phonon_matrix = ex_ph.compute()
```

## Luminescence

### Overview

Combines exciton-dipole and exciton-phonon coupling to compute luminescence properties.

### Constructor

```python
def __init__(self, path=None, save='SAVE', lelph_db=None, latdb=None, wfdb=None, 
             ydipdb=None, bands_range=None, BSE_dir='bse', LELPH_dir='lelph', 
             DIP_dir='gw', save_files=True):
```

### Key Methods

#### compute_luminescence()

```python
def compute_luminescence(self):
```

Compute luminescence properties using combined dipole and phonon matrix elements.

#### compute_luminescence_spectrum()

```python
def compute_luminescence_spectrum(self, ome_range, temp=20, broadening=0.00124, 
                                 npol=2, ph_thr=1e-9):
```

Compute luminescence spectrum intensities.

**Parameters:**
- `ome_range` (`tuple`): Energy range (start, end, num_points)
- `temp` (`float`): Temperature in Kelvin
- `broadening` (`float`): Peak broadening in Hartree
- `npol` (`int`): Number of polarizations
- `ph_thr` (`float`): Phonon frequency threshold

**Example:**
```python
lum = Luminescence(path='./calculation')
results = lum.compute_luminescence()
spectrum = lum.compute_luminescence_spectrum((1.0, 3.0, 1000), temp=300)
```

## Point Group Operations with spgrep

### Overview

The point group operations now use the **spgrep library** for enhanced accuracy and comprehensive analysis. This provides:

- **Automatic point group identification** using crystallographic standards
- **Complete character tables** from comprehensive databases
- **Robust irreducible representation analysis**
- **Fallback to original implementation** if spgrep is not available

### Key Functions

#### get_pg_info()

```python
def get_pg_info(symm_mats):
```

Analyze point group with automatic spgrep/fallback selection.

**Parameters:**
- `symm_mats` (`numpy.ndarray`): Symmetry matrices (nsym, 3, 3)

**Returns:**
- `pg_label` (`str`): Point group label (e.g., 'C2v', 'D3h')
- `classes` (`list`): Symmetry class labels
- `class_dict` (`dict`): Class to operation mapping
- `char_tab` (`numpy.ndarray`): Character table
- `irreps` (`list`): Irreducible representation labels

#### decompose_rep2irrep()

```python
def decompose_rep2irrep(red_rep, char_table, pg_order, class_order, irreps):
```

Decompose reducible representation using spgrep-enhanced analysis.

### Installation

To use spgrep features:
```bash
pip install spgrep
```

If spgrep is not available, the system automatically falls back to the original implementation.

## Utility Functions

The `utils.py` module provides common utility functions:

### Database Utilities

- `read_lelph_database()`: Read electron-phonon database
- `compute_symmetry_matrices()`: Compute symmetry matrices in reduced coordinates
- `create_kpoint_mapping()`: Create k-point mapping arrays

### File Operations

- `validate_path()`: Validate and normalize file paths
- `setup_directories()`: Setup multiple directories
- `safe_file_operation()`: Safe file operations with error handling

### Data Processing

- `process_bands_range()`: Process and validate band ranges
- `convert_energy_units()`: Convert between energy units
- `process_dipoles_by_spin()`: Process dipole elements by spin

### Progress Tracking

- `create_progress_bar()`: Create consistent progress bars

## Best Practices

### Initialization

```python
# Basic usage
egt = ExcitonGroupTheory(path='./calc', BSE_dir='bse', LELPH_dir='lelph')

# Advanced usage with pre-loaded databases
from yambopy.dbs.latticedb import YamboLatticeDB
latdb = YamboLatticeDB.from_db_file('SAVE/ns.db1')
egt = ExcitonGroupTheory(path='./calc', latdb=latdb, bands_range=[1, 20])
```

### Error Handling

```python
try:
    egt = ExcitonGroupTheory(path='./calc')
    results = egt.analyze_exciton_symmetry(iQ=1, nstates=10)
except IOError as e:
    print(f"Database error: {e}")
except ValueError as e:
    print(f"Parameter error: {e}")
```

### Memory Management

```python
# For large systems, limit bands and states
egt = ExcitonGroupTheory(
    path='./calc',
    bands_range=[5, 15],  # Limit band range
    save_files=False      # Disable file saving if not needed
)
```