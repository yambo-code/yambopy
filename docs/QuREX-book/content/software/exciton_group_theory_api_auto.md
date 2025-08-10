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

## Notes

- This documentation is automatically generated from the source code docstrings
- For detailed theoretical background, see the theory section
- For practical examples, see the tutorials and example notebooks
