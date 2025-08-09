# ExcitonGroupTheory API Reference (Auto-Generated)

This page contains automatically generated API documentation from the source code docstrings.

```{eval-rst}
.. currentmodule:: yambopy.optical_properties.exciton_group_theory

.. autoclass:: ExcitonGroupTheory
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__
```

## Class Methods

### Public Methods

```{eval-rst}
.. currentmodule:: yambopy.optical_properties.exciton_group_theory


.. automethod:: ExcitonGroupTheory.analyze_exciton_symmetry

.. automethod:: ExcitonGroupTheory.compute

.. automethod:: ExcitonGroupTheory.read

.. automethod:: ExcitonGroupTheory.read_common_databases

.. automethod:: ExcitonGroupTheory.read_excdb

.. automethod:: ExcitonGroupTheory.read_excdb_single

.. automethod:: ExcitonGroupTheory.save_analysis_results

```

### Internal Methods

These methods are called internally by the main analysis methods:

```{eval-rst}
.. currentmodule:: yambopy.optical_properties.exciton_group_theory


.. automethod:: ExcitonGroupTheory._analyze_with_spglib_symmetries

.. automethod:: ExcitonGroupTheory._build_kpoint_tree

.. automethod:: ExcitonGroupTheory._get_crystal_structure

.. automethod:: ExcitonGroupTheory._get_crystallographic_point_group

.. automethod:: ExcitonGroupTheory._get_point_group_from_spglib

.. automethod:: ExcitonGroupTheory._match_operations_to_point_group

.. automethod:: ExcitonGroupTheory._read_dipoles_db

.. automethod:: ExcitonGroupTheory._read_lattice_db

.. automethod:: ExcitonGroupTheory._read_wavefunction_db

.. automethod:: ExcitonGroupTheory._setup_directories

.. automethod:: ExcitonGroupTheory._setup_kpoint_mapping

.. automethod:: ExcitonGroupTheory._setup_symmetry_data

```

## Usage Examples

### Basic Usage

```python
from yambopy.optical_properties.exciton_group_theory import ExcitonGroupTheory

# Initialize the class
instance = ExcitonGroupTheory(
    path='.',
    save='SAVE'
)

# Use the main methods
results = instance.analyze_exciton_symmetry(iQ=1, nstates=10)
```

## Notes

- This documentation is automatically generated from the source code docstrings
- For detailed theoretical background, see the theory section
- For practical examples, see the tutorials and example notebooks
