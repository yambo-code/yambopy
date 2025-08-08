# Exciton Group Theory Analysis Tutorial

## Overview

This tutorial demonstrates how to use the **significantly improved** `ExcitonGroupTheory` class to analyze the symmetry properties of exciton states in crystalline materials. The implementation has been **completely rewritten** in 2024 to follow the original algorithm exactly, ensuring **maximum accuracy** and **enhanced performance**. We'll walk through a complete example showing how to set up the analysis, interpret the results, and understand the physical implications.

## Prerequisites

Before starting this tutorial, ensure you have:

1. **Yambo calculation completed**: A converged BSE calculation with exciton eigenvalues and eigenvectors
2. **LetzElPhC calculation**: Electron-phonon matrix elements and D-matrices for wavefunction rotation
3. **Required databases**: All necessary Yambo database files in the correct directories
4. **Yambopy installed**: The yambopy package with the ExcitonGroupTheory module
5. **spgrep library (recommended)**: For enhanced point group analysis
   ```bash
   pip install spgrep
   ```
   If not available, the system automatically falls back to the original implementation.

## Step-by-Step Tutorial

### Step 1: Setting Up the Calculation Directory

Organize your calculation directory structure:

```
calculation_directory/
├── SAVE/
│   ├── ns.db1          # Lattice database
│   └── ns.wf           # Wavefunction database
├── bse/
│   └── ndb.BS_diago_Q1 # BSE database for Q-point 1
└── lelph/
    ├── ndb.elph        # Electron-phonon database
    └── ndb.Dmats       # D-matrices
```

### Step 2: Basic Initialization

```python
import numpy as np
from yambopy.optical_properties.exciton_group_theory import ExcitonGroupTheory

# Set up paths and parameters
calculation_path = '.'  # Current directory
SAVE_dir = 'SAVE'
BSE_dir = 'bse'
LELPH_dir = 'lelph'

# Define the band range for analysis
bands_range = [1, 20]  # Analyze bands 1-20

# Initialize the ExcitonGroupTheory class
egt = ExcitonGroupTheory(
    path=calculation_path,
    save=SAVE_dir,
    BSE_dir=BSE_dir,
    LELPH_dir=LELPH_dir,
    bands_range=bands_range,
    read_symm_from_ns_db_file=False  # Use symmetries from elph file
)

print(f"Successfully initialized ExcitonGroupTheory!")
print(f"Number of IBZ k-points: {egt.nibz}")
print(f"Number of symmetry operations: {len(egt.symm_mats)}")
print(f"Time reversal symmetry: {egt.ele_time_rev}")
```

### Step 3: Performing the Group Theory Analysis

```python
# Analysis parameters
iQ = 1              # Q-point index (1-based, as in Yambo)
nstates = 10        # Number of exciton states to analyze
degen_thres = 0.001 # Degeneracy threshold in eV

print(f"Analyzing Q-point {iQ} with {nstates} states...")
print(f"Degeneracy threshold: {degen_thres} eV")

# Perform the analysis
results = egt.analyze_exciton_symmetry(
    iQ=iQ,
    nstates=nstates,
    degen_thres=degen_thres
)

print("Analysis completed successfully!")
```

### Step 4: Interpreting the Results

```python
# Extract key information
q_point = results['q_point']
point_group = results['point_group_label']
little_group = results['little_group']
unique_energies = results['unique_energies']
degeneracies = results['degeneracies']
irrep_decomposition = results['irrep_decomposition']

print("\n" + "="*50)
print("ANALYSIS RESULTS")
print("="*50)
print(f"Q-point coordinates: [{q_point[0]:.4f}, {q_point[1]:.4f}, {q_point[2]:.4f}]")
print(f"Point group: {point_group}")
print(f"Little group operations: {little_group}")
print(f"Number of unique energy levels: {len(unique_energies)}")

print("\nEnergy levels and their symmetries:")
print("-"*60)
print("Level | Energy (eV) | Degeneracy | Irreducible Representation")
print("-"*60)

for i, (energy, degen, irrep) in enumerate(zip(
    unique_energies, degeneracies, irrep_decomposition)):
    print(f"{i+1:5d} | {energy:10.4f} | {degen:10d} | {irrep}")
```

### Step 5: Advanced Analysis

#### Analyzing Optical Selection Rules

```python
def analyze_optical_selection_rules(results):
    """
    Analyze optical selection rules based on irreducible representations.
    """
    print("\n" + "="*50)
    print("OPTICAL SELECTION RULES ANALYSIS")
    print("="*50)
    
    irreps = results['irrep_decomposition']
    energies = results['unique_energies']
    
    # Common optically active irreps (this depends on the point group)
    optically_active = ['A1', 'E', 'T1u', 'T2u']  # Example for common groups
    
    print("Optical activity analysis:")
    print("-"*40)
    
    for i, (energy, irrep) in enumerate(zip(energies, irreps)):
        is_active = any(active in irrep for active in optically_active)
        activity = "BRIGHT" if is_active else "DARK"
        print(f"Level {i+1}: {energy:.4f} eV ({irrep}) -> {activity}")

# Perform optical analysis
analyze_optical_selection_rules(results)
```

#### Symmetry Breaking Analysis

```python
def analyze_symmetry_breaking(results, reference_point_group="D6h"):
    """
    Analyze potential symmetry breaking by comparing with a reference.
    """
    current_pg = results['point_group_label']
    
    print("\n" + "="*50)
    print("SYMMETRY BREAKING ANALYSIS")
    print("="*50)
    print(f"Current point group: {current_pg}")
    print(f"Reference point group: {reference_point_group}")
    
    if current_pg != reference_point_group:
        print("⚠️  Symmetry breaking detected!")
        print("Possible causes:")
        print("  - Strain effects")
        print("  - External fields")
        print("  - Surface/interface effects")
        print("  - Structural distortions")
    else:
        print("✅ No symmetry breaking detected")

# Analyze symmetry breaking
analyze_symmetry_breaking(results, "C2v")  # Compare with C2v symmetry
```

### Step 6: Saving and Visualizing Results

```python
# Save detailed results to file
output_filename = f"exciton_symmetry_Q{iQ}_analysis.txt"
egt.save_analysis_results(results, output_filename)
print(f"\nDetailed results saved to: {output_filename}")

# Create a summary plot (requires matplotlib)
try:
    import matplotlib.pyplot as plt
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot energy levels with degeneracies
    ax1.scatter(range(len(unique_energies)), unique_energies, 
                s=degeneracies*50, alpha=0.7)
    ax1.set_xlabel('Level Index')
    ax1.set_ylabel('Energy (eV)')
    ax1.set_title('Exciton Energy Levels\n(size ∝ degeneracy)')
    ax1.grid(True, alpha=0.3)
    
    # Plot degeneracy distribution
    ax2.bar(range(len(degeneracies)), degeneracies, alpha=0.7)
    ax2.set_xlabel('Level Index')
    ax2.set_ylabel('Degeneracy')
    ax2.set_title('Degeneracy Distribution')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'exciton_symmetry_Q{iQ}_plot.png', dpi=300, bbox_inches='tight')
    print(f"Plot saved as: exciton_symmetry_Q{iQ}_plot.png")
    
except ImportError:
    print("Matplotlib not available - skipping plot generation")
```

### Step 7: Multiple Q-point Analysis

```python
def analyze_multiple_q_points(egt, q_points, nstates=10, degen_thres=0.001):
    """
    Analyze multiple Q-points and compare their symmetries.
    """
    all_results = {}
    
    print("\n" + "="*50)
    print("MULTIPLE Q-POINT ANALYSIS")
    print("="*50)
    
    for iQ in q_points:
        print(f"\nAnalyzing Q-point {iQ}...")
        try:
            results = egt.analyze_exciton_symmetry(iQ, nstates, degen_thres)
            all_results[iQ] = results
            
            print(f"  Point group: {results['point_group_label']}")
            print(f"  Number of levels: {len(results['unique_energies'])}")
            
        except Exception as e:
            print(f"  Error analyzing Q-point {iQ}: {e}")
    
    return all_results

# Analyze multiple Q-points (if available)
q_points_to_analyze = [1, 2, 3]  # Adjust based on your calculation
try:
    multi_results = analyze_multiple_q_points(egt, q_points_to_analyze)
except:
    print("Multiple Q-point analysis not available")
```

## Understanding the Output

### Point Group Identification

The analysis identifies the little group of the Q-point, which determines the allowed symmetries. Common point groups include:

- **C1**: No symmetry (all states non-degenerate)
- **C2v**: Mirror symmetries (A1, A2, B1, B2 irreps)
- **D3h**: Hexagonal symmetry (A1', A2', E', A1", A2", E" irreps)
- **Oh**: Cubic symmetry (A1g, A2g, Eg, T1g, T2g, A1u, A2u, Eu, T1u, T2u irreps)

### Irreducible Representations

Each energy level is classified by its irreducible representation:

- **A-type**: Non-degenerate, symmetric under rotations
- **B-type**: Non-degenerate, antisymmetric under some rotations
- **E-type**: Doubly degenerate
- **T-type**: Triply degenerate
- **g/u subscripts**: Even/odd under inversion
- **'/'' subscripts**: Even/odd under horizontal reflection

### Physical Interpretation

1. **Bright vs. Dark Excitons**: Certain irreps are optically active (bright), others are forbidden (dark)
2. **Polarization Dependence**: Different irreps couple to different light polarizations
3. **Fine Structure**: Degeneracy lifting reveals fine structure splitting
4. **Selection Rules**: Determines allowed optical transitions

## Troubleshooting

### Common Issues

1. **File Not Found Errors**:
   ```python
   # Check if all required files exist
   import os
   required_files = [
       'SAVE/ns.db1', 'SAVE/ns.wf',
       f'{BSE_dir}/ndb.BS_diago_Q{iQ}',
       f'{LELPH_dir}/ndb.elph', f'{LELPH_dir}/ndb.Dmats'
   ]
   
   for file in required_files:
       if not os.path.exists(file):
           print(f"Missing file: {file}")
   ```

2. **Memory Issues**:
   ```python
   # Reduce the number of states or bands
   nstates = 5  # Instead of 10
   bands_range = [1, 10]  # Instead of [1, 20]
   ```

3. **Convergence Issues**:
   ```python
   # Adjust degeneracy threshold
   degen_thres = 0.01  # Less strict threshold
   ```

### Performance Optimization

```python
# For large systems, consider:
# 1. Limiting the number of states
# 2. Using a coarser k-point grid
# 3. Reducing the band range
# 4. Using symmetry to reduce calculations

# Example optimized setup
egt_optimized = ExcitonGroupTheory(
    path=calculation_path,
    save=SAVE_dir,
    BSE_dir=BSE_dir,
    LELPH_dir=LELPH_dir,
    bands_range=[1, 10],  # Reduced band range
    read_symm_from_ns_db_file=True  # Faster symmetry reading
)
```

## Advanced Applications

### Custom Point Group Analysis

```python
def custom_point_group_analysis(results):
    """
    Implement custom point group analysis for specific materials.
    """
    # Example for 2D materials with specific symmetries
    pg = results['point_group_label']
    
    if 'D3h' in pg or 'D6h' in pg:
        print("Hexagonal 2D material detected")
        print("Expected irreps: A1', A2', E1', E2', A1'', A2'', E1'', E2''")
    elif 'C2v' in pg:
        print("Rectangular 2D material detected")
        print("Expected irreps: A1, A2, B1, B2")
    
    # Add material-specific analysis here
```

### Integration with Other Yambopy Modules

```python
# Example: Combine with optical absorption analysis
from yambopy.bse import YamboExcitonDB

# Read additional BSE data
bse_db = YamboExcitonDB.from_db_file(egt.ydb, folder=BSE_dir, 
                                     filename=f'ndb.BS_diago_Q{iQ}')

# Combine symmetry analysis with oscillator strengths
oscillator_strengths = bse_db.get_oscillator_strengths()

print("Combined analysis:")
for i, (energy, irrep, osc_str) in enumerate(zip(
    results['unique_energies'], 
    results['irrep_decomposition'],
    oscillator_strengths[:len(results['unique_energies'])])):
    print(f"Level {i+1}: {energy:.4f} eV, {irrep}, f = {osc_str:.4f}")
```

## Conclusion

The ExcitonGroupTheory class provides a powerful tool for understanding the symmetry properties of exciton states. By following this tutorial, you can:

1. Set up and perform group theory analysis
2. Interpret the symmetry classification of exciton states
3. Understand optical selection rules and fine structure
4. Identify symmetry breaking effects
5. Integrate the analysis with other computational tools

This analysis is essential for understanding the optical properties of materials and designing new materials with desired symmetry properties.