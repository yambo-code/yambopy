# ExcitonGroupTheory API Reference (Auto-Generated)

This page contains automatically generated API documentation from the source code docstrings.

```{eval-rst}
.. currentmodule:: yambopy.optical_properties

.. autoclass:: ExcitonGroupTheory
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__
```

## Class Methods

### Public Methods

```{eval-rst}
.. currentmodule:: yambopy.optical_properties

.. automethod:: ExcitonGroupTheory.analyze_exciton_symmetry

.. automethod:: ExcitonGroupTheory.classify_symmetry_operations

.. automethod:: ExcitonGroupTheory.display_symmetry_operations

.. automethod:: ExcitonGroupTheory.get_latex_labels

.. automethod:: ExcitonGroupTheory.read

```

### Internal Methods

These methods are called internally by the main analysis methods:

```{eval-rst}
.. currentmodule:: yambopy.optical_properties

.. automethod:: ExcitonGroupTheory._setup_symmetry

.. automethod:: ExcitonGroupTheory._compute_spglib_dmats

.. automethod:: ExcitonGroupTheory._analyze_activity

.. automethod:: ExcitonGroupTheory._read_bse_data

.. automethod:: ExcitonGroupTheory._group_degenerate_states

.. automethod:: ExcitonGroupTheory._compute_representation_matrices

.. automethod:: ExcitonGroupTheory._get_crystal_system

```

## Usage Examples

### Basic Usage

```python
from yambopy.optical_properties import ExcitonGroupTheory

# Initialize the class
egt = ExcitonGroupTheory(
    path='./',
    save='SAVE',
    BSE_dir='./GW_BSE/bse',
    LELPH_dir='./lelph',
    bands_range=[6, 10]
)

# Perform symmetry analysis
results = egt.analyze_exciton_symmetry(iQ=1, nstates=10)

# Get LaTeX labels for plotting
latex_labels = egt.get_latex_labels(['A1g', 'E2u', 'B1u'])
```

### General Symmetry Classification

```python
# Universal symmetry operation classification (works for all 230 space groups)
operations = egt.classify_symmetry_operations()
summary = operations.get('_summary', {})

print(f"Space Group: {summary.get('space_group')} (#{summary.get('space_group_number')})")
print(f"Crystal System: {summary.get('crystal_system')}")

# Show operation breakdown
for op_type, op_list in operations.items():
    if op_type != '_summary' and op_list:
        print(f"{op_type.title()}: {len(op_list)} operations")

# Display comprehensive analysis
egt.display_symmetry_operations()
```

### Crystal System Examples

```python
# Works with any crystal system
# Hexagonal (hBN): P6₃/mmc
# Cubic (diamond): Fd3m  
# Tetragonal (TiO₂): P4₂/mnm
# Orthorhombic (Pnma): Pnma
# Monoclinic (β-Ga₂O₃): C2/m
# Triclinic (CuSO₄·5H₂O): P-1

# The same code works for all systems!
operations = egt.classify_symmetry_operations()
crystal_system = operations['_summary']['crystal_system']
print(f"Detected {crystal_system} crystal system")
```

## Key Methods

### `classify_symmetry_operations()`
- **Universal**: Works with all 230 space groups
- **Comprehensive**: Includes symmorphic and non-symmorphic operations
- **Accurate**: Uses spglib for crystallographic validation
- **Detailed**: Provides operation matrices, symbols, and translations

### `display_symmetry_operations()`
- **Publication-ready**: Formatted output with proper notation
- **Educational**: Includes crystal system characteristics
- **Complete**: Shows all operation details and classifications
