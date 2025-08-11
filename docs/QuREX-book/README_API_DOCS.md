# Automatic API Documentation System

This directory contains an automatic API documentation system that generates documentation directly from Python docstrings in the source code.

## Overview

The API documentation is now **automatically generated** from the docstrings in the source code, ensuring that:

- ✅ Documentation stays synchronized with code changes
- ✅ No manual maintenance of API files required
- ✅ Comprehensive coverage of all methods and classes
- ✅ Consistent formatting and structure

## Files

### Auto-Generation Scripts

- **`generate_api_docs.py`**: Main script to generate API documentation from docstrings
- **`build_docs.py`**: Complete build script that generates API docs and builds the book
- **`README_API_DOCS.md`**: This documentation file

### Generated Documentation

- **`content/software/exciton_group_theory_api_auto.md`**: Auto-generated ExcitonGroupTheory API (**UPDATED 2024**)
- **`content/software/point_group_operations_api_auto.md`**: Auto-generated point group operations API

### Manual Documentation (Updated 2024)

- **`content/theory/exciton_group_theory.md`**: Theoretical background (**UPDATED** - universal space group support)
- **`content/software/exciton_group_theory_summary.md`**: Comprehensive summary (**NEW** - complete feature overview)
- **`content/software/yambopy_improvements_2024.md`**: 2024 improvements documentation (**NEW**)
- **`notebooks/exciton_group_theory_*.ipynb`**: Example notebooks (**UPDATED** - new universal features)

## Usage

### Quick Build

```bash
# Generate API docs and build the complete book
python build_docs.py

# Only generate API documentation
python build_docs.py --api-only

# Clean build and regenerate everything
python build_docs.py --clean
```

### Manual API Generation

```bash
# Generate only the API documentation
python generate_api_docs.py
```

## How It Works

### 1. Docstring Extraction

The `generate_api_docs.py` script:

1. **Imports the Python modules** using `importlib`
2. **Extracts class and method information** using `inspect`
3. **Generates Markdown files** with Sphinx `autodoc` directives
4. **Organizes methods** into public and internal categories

### 2. Sphinx Integration

The generated Markdown files use Sphinx `eval-rst` blocks with:

- `.. autoclass::` for complete class documentation
- `.. automethod::` for individual method documentation
- `.. automodule::` for module-level documentation

### 3. Jupyter Book Integration

The Jupyter Book configuration (`conf.py`) includes:

- `sphinx.ext.autodoc`: Automatic documentation generation
- `sphinx.ext.autosummary`: Summary tables
- `sys_path = ['../../yambopy']`: Python path configuration

## Benefits

### For Developers

- **No manual API maintenance**: Just write good docstrings
- **Automatic synchronization**: Documentation updates with code changes
- **Consistent formatting**: Standardized API documentation structure
- **Easy updates**: Run one script to regenerate everything

### For Users

- **Always up-to-date**: API docs reflect current code state
- **Comprehensive coverage**: All methods and classes documented
- **Rich formatting**: Full Sphinx features (math, cross-references, etc.)
- **Integrated experience**: API docs seamlessly integrated with tutorials

## Docstring Guidelines

To ensure high-quality auto-generated documentation, follow these guidelines:

### Class Docstrings

```python
class ExcitonGroupTheory(BaseOpticalProperties):
    """
    Universal group theory analysis of exciton states for all 230 space groups.
    
    **NEW 2024**: Complete rewrite with universal space group support using spglib.
    This class performs comprehensive symmetry analysis including non-symmorphic 
    operations (screw rotations, glide reflections) for any crystal system.
    
    **Key Features**
    - Universal space group support (all 230 space groups)
    - Non-symmorphic operations (screw rotations, glide reflections)
    - Automatic crystal system detection
    - Professional crystallographic accuracy via spglib
    - Clean implementation with no duplicate methods
    
    **Theoretical Background**
    
    The symmetry of exciton states is determined by the little group G_k:
    
        D^(n)_R = ⟨ψ_n(Rk)| U(R) |ψ_n(k)⟩
        χ^(n)(R) = Tr[D^(n)_R]
    
    Parameters
    ----------
    path : str, optional
        Path to the calculation directory.
    BSE_dir : str, optional
        Name of the BSE directory. Default is 'bse'.
    LELPH_dir : str, optional
        Name of the electron-phonon directory. Default is 'lelph'.
    bands_range : list, optional
        Range of bands to include in the analysis.
    
    Attributes
    ----------
    point_group_label : str
        Identified crystallographic point group.
    spacegroup_label : str
        Identified space group.
    
    Examples
    --------
    >>> from yambopy.optical_properties import ExcitonGroupTheory
    >>> egt = ExcitonGroupTheory(path='./', BSE_dir='bse', LELPH_dir='lelph')
    
    # NEW: Universal symmetry classification
    >>> operations = egt.classify_symmetry_operations()
    >>> summary = operations['_summary']
    >>> print(f"Space Group: {summary['space_group']} (#{summary['space_group_number']})")
    
    # Legacy: Exciton group theory analysis
    >>> results = egt.analyze_exciton_symmetry(iQ=1, nstates=10)
    >>> latex_labels = egt.get_latex_labels(['A1g', 'E2u'])
    """
```

### Method Docstrings

```python
def analyze_exciton_symmetry(self, iQ, nstates, degen_thres=0.001):
    """
    Brief description of what the method does.
    
    Longer description with implementation details and mathematical
    background if relevant.
    
    Parameters
    ----------
    iQ : int
        Q-point index (1-based indexing as in Yambo).
    nstates : int
        Number of exciton states to analyze.
    degen_thres : float, optional
        Degeneracy threshold in eV. Default is 0.001.
    
    Returns
    -------
    results : dict
        Dictionary containing analysis results with keys:
        - 'point_group_label': Point group symbol
        - 'irrep_decomposition': List of irrep labels
        - 'optical_activity': Optical selection rules
    
    Raises
    ------
    IOError
        If required database files cannot be read.
    ValueError
        If parameters are invalid.
    
    Examples
    --------
    >>> results = egt.analyze_exciton_symmetry(iQ=1, nstates=5)
    >>> print(results['point_group_label'])
    'D3h'
    
    Notes
    -----
    This method implements the complete group theory analysis workflow
    following the algorithm described in [1]_.
    
    References
    ----------
    .. [1] Author, "Title", Journal, Year.
    """
```

## Maintenance

### Adding New Classes/Modules

To add documentation for new classes or modules:

1. **Add to `generate_api_docs.py`**:
   ```python
   # In main() function
   generate_class_api_doc(
       'yambopy.new_module.new_class',
       'NewClass',
       output_dir / 'new_class_api_auto.md'
   )
   ```

2. **Add to `_toc.yml`**:
   ```yaml
   - file: content/software/new_class_api_auto
   ```

3. **Regenerate documentation**:
   ```bash
   python build_docs.py
   ```

### Updating Existing Documentation

1. **Modify docstrings** in the source code
2. **Run the build script**:
   ```bash
   python build_docs.py
   ```
3. **Documentation is automatically updated**

## Future Enhancements

Possible improvements to the system:

1. **Cross-references**: Automatic linking between related methods
2. **Type hints**: Enhanced parameter documentation from type annotations
3. **Examples extraction**: Automatic extraction of examples from test files
4. **Performance metrics**: Documentation of method performance characteristics
5. **Version tracking**: Documentation of API changes between versions

## Troubleshooting

### Common Issues

**Import errors during generation:**
- Check that `sys.path` includes the yambopy directory
- Ensure all dependencies are installed
- Verify module names are correct

**Missing methods in documentation:**
- Check that methods have proper docstrings
- Verify method names in the generation script
- Ensure methods are not filtered out by naming conventions

**Sphinx build errors:**
- Check that `conf.py` has correct extensions enabled
- Verify that autodoc can import the modules
- Check for syntax errors in generated RST directives

### Getting Help

For issues with the documentation system:

1. Check the build output for specific error messages
2. Verify that the source code docstrings are properly formatted
3. Test the generation script independently: `python generate_api_docs.py`
4. Check the Jupyter Book build logs for Sphinx errors

---
# Automatic API Documentation System

This directory contains an automatic API documentation system that generates documentation directly from Python docstrings in the source code.

## Overview

The API documentation is now **automatically generated** from the docstrings in the source code, ensuring that:

- ✅ Documentation stays synchronized with code changes
- ✅ No manual maintenance of API files required
- ✅ Comprehensive coverage of all methods and classes
- ✅ Consistent formatting and structure

## Files

### Auto-Generation Scripts

- **`generate_api_docs.py`**: Main script to generate API documentation from docstrings
- **`build_docs.py`**: Complete build script that generates API docs and builds the book
- **`README_API_DOCS.md`**: This documentation file

### Generated Documentation

- **`content/software/exciton_group_theory_api_auto.md`**: Auto-generated ExcitonGroupTheory API
- **`content/software/point_group_operations_api_auto.md`**: Auto-generated point group operations API

### Manual Documentation (Still Used)

- **`content/theory/exciton_group_theory.md`**: Theoretical background (manual)
- **`content/tutorials/exciton_group_theory_tutorial.md`**: Tutorial (manual)
- **`notebooks/exciton_group_theory_*.ipynb`**: Example notebooks (manual)

## Usage

### Quick Build

```bash
# Generate API docs and build the complete book
python build_docs.py

# Only generate API documentation
python build_docs.py --api-only

# Clean build and regenerate everything
python build_docs.py --clean
```

### Manual API Generation

```bash
# Generate only the API documentation
python generate_api_docs.py
```

## How It Works

### 1. Docstring Extraction

The `generate_api_docs.py` script:

1. **Imports the Python modules** using `importlib`
2. **Extracts class and method information** using `inspect`
3. **Generates Markdown files** with Sphinx `autodoc` directives
4. **Organizes methods** into public and internal categories

### 2. Sphinx Integration

The generated Markdown files use Sphinx `eval-rst` blocks with:

- `.. autoclass::` for complete class documentation
- `.. automethod::` for individual method documentation
- `.. automodule::` for module-level documentation

### 3. Jupyter Book Integration

The Jupyter Book configuration (`conf.py`) includes:

- `sphinx.ext.autodoc`: Automatic documentation generation
- `sphinx.ext.autosummary`: Summary tables
- `sys_path = ['../../yambopy']`: Python path configuration

## Benefits

### For Developers

- **No manual API maintenance**: Just write good docstrings
- **Automatic synchronization**: Documentation updates with code changes
- **Consistent formatting**: Standardized API documentation structure
- **Easy updates**: Run one script to regenerate everything

### For Users

- **Always up-to-date**: API docs reflect current code state
- **Comprehensive coverage**: All methods and classes documented
- **Rich formatting**: Full Sphinx features (math, cross-references, etc.)
- **Integrated experience**: API docs seamlessly integrated with tutorials

## Docstring Guidelines

To ensure high-quality auto-generated documentation, follow these guidelines:

### Class Docstrings

```python
class ExcitonGroupTheory(BaseOpticalProperties):
    """
    Brief one-line description.
    
    Longer description with multiple paragraphs explaining the purpose,
    theoretical background, and key features.
    
    **Theoretical Background**
    
    Mathematical formulations and equations using LaTeX:
    
        D^(n)_R = ⟨ψ_n(Rk)| U(R) |ψ_n(k)⟩
    
    Parameters
    ----------
    param1 : type
        Description of parameter 1.
    param2 : type, optional
        Description of optional parameter 2.
    
    Attributes
    ----------
    attr1 : type
        Description of attribute 1.
    
    Examples
    --------
    >>> egt = ExcitonGroupTheory(path='.')
    >>> results = egt.analyze_exciton_symmetry(iQ=1, nstates=10)
    """
```

### Method Docstrings

```python
def analyze_exciton_symmetry(self, iQ, nstates, degen_thres=0.001):
    """
    Brief description of what the method does.
    
    Longer description with implementation details and mathematical
    background if relevant.
    
    Parameters
    ----------
    iQ : int
        Q-point index (1-based indexing as in Yambo).
    nstates : int
        Number of exciton states to analyze.
    degen_thres : float, optional
        Degeneracy threshold in eV. Default is 0.001.
    
    Returns
    -------
    results : dict
        Dictionary containing analysis results with keys:
        - 'point_group_label': Point group symbol
        - 'irrep_decomposition': List of irrep labels
        - 'optical_activity': Optical selection rules
    
    Raises
    ------
    IOError
        If required database files cannot be read.
    ValueError
        If parameters are invalid.
    
    Examples
    --------
    >>> results = egt.analyze_exciton_symmetry(iQ=1, nstates=5)
    >>> print(results['point_group_label'])
    'D3h'
    
    Notes
    -----
    This method implements the complete group theory analysis workflow
    following the algorithm described in [1]_.
    
    References
    ----------
    .. [1] Author, "Title", Journal, Year.
    """
```

## Maintenance

### Adding New Classes/Modules

To add documentation for new classes or modules:

1. **Add to `generate_api_docs.py`**:
   ```python
   # In main() function
   generate_class_api_doc(
       'yambopy.new_module.new_class',
       'NewClass',
       output_dir / 'new_class_api_auto.md'
   )
   ```

2. **Add to `_toc.yml`**:
   ```yaml
   - file: content/software/new_class_api_auto
   ```

3. **Regenerate documentation**:
   ```bash
   python build_docs.py
   ```

### Updating Existing Documentation

1. **Modify docstrings** in the source code
2. **Run the build script**:
   ```bash
   python build_docs.py
   ```
3. **Documentation is automatically updated**

## Future Enhancements

Possible improvements to the system:

1. **Cross-references**: Automatic linking between related methods
2. **Type hints**: Enhanced parameter documentation from type annotations
3. **Examples extraction**: Automatic extraction of examples from test files
4. **Performance metrics**: Documentation of method performance characteristics
5. **Version tracking**: Documentation of API changes between versions