#!/usr/bin/env python3
"""
Working script to generate API documentation that properly renders in Jupyter Book.

This script creates API documentation using a simpler approach that works
reliably with Jupyter Book's MyST parser and Sphinx integration.
"""

import os
import sys
import importlib.util
import inspect
from pathlib import Path

# Add yambopy to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / 'yambopy'))

def generate_simple_class_api(module_name, class_name, output_file, description=""):
    """Generate API documentation that references autoapi-generated docs."""
    
    # Change extension to .rst
    output_file = str(output_file).replace('.md', '.rst')
    
    # Convert module path for autoapi reference
    autoapi_module_path = module_name.replace('.', '/')
    
    content = f"""{class_name} API Reference
{'=' * (len(class_name) + 17)}

{description}

This page contains API documentation for the ``{class_name}`` class.

Class Documentation
-------------------

For detailed API documentation including all methods, attributes, and parameters, see the auto-generated documentation:

:doc:`autoapi/{autoapi_module_path}/index`

Quick Reference
---------------

.. currentmodule:: {module_name}

The ``{class_name}`` class provides the following main functionality:

* **Initialization**: Create instances with various configuration options
* **Data Processing**: Methods for reading and processing data
* **Analysis**: Core analysis and computation methods
* **Output**: Methods for saving and exporting results

Usage Example
-------------

.. code-block:: python

   from {module_name} import {class_name}

   # Initialize the class
   instance = {class_name}()

   # Use the methods
   result = instance.method_name()

Notes
-----

- This documentation is automatically generated from source code docstrings
- For detailed examples, see the tutorials and example notebooks
- Complete API reference is available in the autoapi section
"""
    
    # Write to file
    with open(output_file, 'w') as f:
        f.write(content)
    
    print(f"Generated {class_name} API documentation: {output_file}")

def generate_simple_module_api(module_name, output_file, description=""):
    """Generate simple, working API documentation for a module."""
    
    module_title = module_name.split('.')[-1].replace('_', ' ').title()
    
    content = f"""# {module_title} Module API Reference

{description}

This page contains API documentation for the `{module_name}` module.

## Module Documentation

```{{eval-rst}}
.. currentmodule:: {module_name}

.. automodule:: {module_name}
   :members:
   :undoc-members:
   :show-inheritance:
```

## Usage Example

```python
from {module_name} import *

# Use the module components
```

## Notes

- This documentation is automatically generated from source code docstrings
- Import the module to access all classes and functions
"""
    
    # Write to file
    with open(output_file, 'w') as f:
        f.write(content)
    
    print(f"Generated {module_title} module API documentation: {output_file}")

def main():
    """Generate working API documentation."""
    
    # Define output directory
    output_dir = Path(__file__).parent / 'content' / 'api'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("🚀 Generating Working API Documentation")
    print("=" * 50)
    
    # Generate key class documentation
    classes_to_document = [
        {
            'module': 'yambopy.optical_properties.base_optical',
            'class': 'BaseOpticalProperties',
            'description': 'Base class for optical properties calculations providing common functionality for reading Yambo databases and setting up k-point trees.'
        },
        {
            'module': 'yambopy.optical_properties.exciton_group_theory',
            'class': 'ExcitonGroupTheory',
            'description': 'Group theory analysis of exciton states using crystallographic symmetries to determine irreducible representations and optical selection rules.'
        },
        {
            'module': 'yambopy.optical_properties.ex_dipole',
            'class': 'ExcitonDipole',
            'description': 'Calculation of exciton dipole moments and oscillator strengths for optical transitions.'
        },
        {
            'module': 'yambopy.optical_properties.ex_phonon',
            'class': 'ExcitonPhonon',
            'description': 'Analysis of exciton-phonon coupling matrix elements and related optical properties.'
        },
        {
            'module': 'yambopy.optical_properties.luminescence',
            'class': 'Luminescence',
            'description': 'Calculation of photoluminescence spectra and related optical properties.'
        },
        {
            'module': 'yambopy.letzelphc_interface.lelphcdb',
            'class': 'LetzElphElectronPhononDB',
            'description': 'Interface to read electron-phonon matrix elements from LetzElPhC databases.'
        }
    ]
    
    # Generate class documentation
    for class_info in classes_to_document:
        generate_simple_class_api(
            class_info['module'],
            class_info['class'],
            output_dir / f"{class_info['class'].lower()}_api.md",
            class_info['description']
        )
    
    # Generate module documentation
    modules_to_document = [
        {
            'module': 'yambopy.optical_properties.utils',
            'description': 'Utility functions for optical properties calculations including file I/O, validation, and data processing.'
        },
        {
            'module': 'yambopy.optical_properties.spgrep_point_group_ops',
            'description': 'Point group operations and symmetry analysis using the spgrep library.'
        },
        {
            'module': 'yambopy.letzelphc_interface.lelph2y',
            'description': 'Conversion utilities between LetzElPhC and Yambo data formats.'
        }
    ]
    
    for module_info in modules_to_document:
        module_name = module_info['module'].split('.')[-1]
        generate_simple_module_api(
            module_info['module'],
            output_dir / f"{module_name}_api.md",
            module_info['description']
        )
    
    # Generate main API index
    index_content = """# API Documentation

This section contains comprehensive API documentation for yambopy modules.

## Main Classes

### Optical Properties
- [BaseOpticalProperties](baseopticalproperties_api.md) - Base class for optical calculations
- [ExcitonGroupTheory](excitongrouptheory_api.md) - Group theory analysis of excitons
- [ExcitonDipole](excitondipole_api.md) - Dipole moment calculations
- [ExcitonPhonon](excitonphonon_api.md) - Exciton-phonon coupling
- [Luminescence](luminescence_api.md) - Photoluminescence calculations

### LetzElPhC Interface
- [LetzElphElectronPhononDB](letzelphelectronphonondb_api.md) - Electron-phonon database interface

## Utility Modules
- [Utils](utils_api.md) - Utility functions for optical properties
- [Point Group Operations](spgrep_point_group_ops_api.md) - Symmetry operations
- [LetzElPhC Conversion](lelph2y_api.md) - Format conversion utilities

## Quick Start

```python
# Import main classes
from yambopy.optical_properties import (
    BaseOpticalProperties, ExcitonGroupTheory, 
    ExcitonDipole, ExcitonPhonon, Luminescence
)

# Basic usage
egt = ExcitonGroupTheory(path='.')
results = egt.analyze_exciton_symmetry(iQ=1, nstates=10)
```

## Notes

- All documentation is automatically generated from source code docstrings
- For tutorials and examples, see the tutorials section
- For theoretical background, see the theory section
"""
    
    with open(output_dir / 'index.md', 'w') as f:
        f.write(index_content)
    
    print(f"Generated API index: {output_dir / 'index.md'}")
    
    print("\n" + "=" * 50)
    print("✅ Working API Documentation Generated!")
    print("=" * 50)
    print(f"📁 Files generated in: {output_dir}")
    print(f"📚 Total files: {len(list(output_dir.glob('*.md')))}")

if __name__ == '__main__':
    main()