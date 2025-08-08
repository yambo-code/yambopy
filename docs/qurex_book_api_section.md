# QuREX-Book API Section: Exciton Group Theory

## ExcitonGroupTheory Class

Automated symmetry analysis of exciton states using crystallographic group theory. Integrates Yambo databases with spglib/spgrep for comprehensive point group analysis.

### Key Features

- **Automatic crystal structure reading** from Yambo SAVE database
- **spglib space group identification** from atomic structure  
- **spgrep point group classification** using crystallographic databases
- **Systematic irrep decomposition** with universal Γ_n labeling
- **Magnetic system support** through reduced symmetry analysis

### Constructor

```python
ExcitonGroupTheory(path=None, save='SAVE', BSE_dir='bse', 
                   LELPH_dir='lelph', bands_range=None)
```

### Main Method

```python
analyze_exciton_symmetry(iQ, nstates, degen_thres=0.001)
```

**Returns:** Dictionary with point group label, irrep decomposition, and symmetry analysis.

### Example

```python
from yambopy.optical_properties.exciton_group_theory import ExcitonGroupTheory

# Initialize with automatic database reading
egt = ExcitonGroupTheory(path='/calc/path', BSE_dir='bse')

# Analyze Γ-point excitons  
results = egt.analyze_exciton_symmetry(iQ=1, nstates=5)

print(f"Point group: {results['point_group_label']}")
# Output: Point group: 6/mmm (D6h)

for E, deg, irrep in zip(results['unique_energies'], 
                        results['degeneracies'],
                        results['irrep_decomposition']):
    print(f"{E:.3f} eV (deg={deg}) → {irrep}")
# Output: 7.364 eV (deg=2) → Γ₅
```

### Integration Architecture

- **Yambo symmetries**: Wavefunction rotations and matrix elements
- **spglib**: Space group identification from crystal structure
- **spgrep**: Point group classification and irrep generation

This dual approach ensures compatibility with Yambo while leveraging comprehensive crystallographic databases for accurate symmetry analysis.

### Supported Systems

- **Non-magnetic crystals**: Full point group symmetry analysis
- **Magnetic systems**: Reduced symmetry with magnetic point groups  
- **2D materials**: Automatic handling of layered structures
- **All crystal systems**: Triclinic to cubic, including non-symmorphic groups

### References

- spglib: Togo & Tanaka, arXiv:1808.01590 (2018)
- spgrep: Watanabe et al., J. Open Source Softw. 8, 5269 (2023)