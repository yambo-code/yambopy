#!/usr/bin/env python3
"""
Exciton Group Theory Analysis Example

This script demonstrates how to use the ExcitonGroupTheory class to analyze 
the symmetry properties of exciton states in crystalline materials.

Features:
- Automatic point group identification using spglib
- Irreducible representation decomposition using spgrep
- Optical activity analysis (Raman, IR, electric dipole)
- LaTeX formatting for publication-quality plots
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Add development version to path (if needed)
# sys.path.insert(0, '/path/to/yambopy-project')

from yambopy.optical_properties import ExcitonGroupTheory

def main():
    """Main analysis function."""
    
    # Set up plotting parameters
    plt.rcParams['figure.figsize'] = (12, 8)
    plt.rcParams['font.size'] = 12
    
    print("=" * 60)
    print("EXCITON GROUP THEORY ANALYSIS")
    print("=" * 60)
    
    # Configuration
    # Adjust these paths according to your calculation setup
    path = './'  # Current directory or path to your calculation
    SAVE_dir = './SAVE'
    BSE_dir = './GW_BSE/bse'  # or 'GW_BSE' depending on your setup
    LELPH_dir = './lelph'  # Directory containing electron-phonon data
    
    # Exciton analysis parameters
    iQ = 1  # Exciton Q-point index (1-based, as in Yambo)
    nstates = 10  # Number of exciton states to analyze
    degen_thres = 0.001  # Degeneracy threshold in eV
    
    # Band range (adjust according to your calculation)
    bands_range = [6, 10]  # Example: bands 6 to 10
    
    print("Configuration:")
    print(f"  Path: {path}")
    print(f"  BSE directory: {BSE_dir}")
    print(f"  LELPH directory: {LELPH_dir}")
    print(f"  Analyzing Q-point: {iQ}")
    print(f"  Number of states: {nstates}")
    print("-" * 50)
    
    # Initialize the ExcitonGroupTheory class
    print("\nInitializing ExcitonGroupTheory class...")
    print("This will read the database files and set up symmetry operations.")
    
    try:
        egt = ExcitonGroupTheory(
            path=path,
            save='SAVE',
            BSE_dir=BSE_dir,
            LELPH_dir=LELPH_dir,
            bands_range=bands_range,
            read_symm_from_ns_db_file=True  # Read symmetries from ns.db1
        )
        
        print("\nInitialization completed successfully!")
        print(f"  Point group: {egt.point_group_label}")
        print(f"  Space group: {egt.spacegroup_label}")
        print(f"  Number of symmetry operations: {len(egt.symm_mats)}")
        print(f"  Number of IBZ k-points: {egt.nibz}")
        print(f"  Number of bands: {egt.nbands}")
        
    except Exception as e:
        print(f"Error during initialization: {e}")
        return
    
    # Perform the symmetry analysis
    print("\n" + "=" * 60)
    print("EXCITON SYMMETRY ANALYSIS")
    print("=" * 60)
    
    try:
        results = egt.analyze_exciton_symmetry(
            iQ=iQ, 
            nstates=nstates, 
            degen_thres=degen_thres
        )
        
        print(f"\nCrystal Structure:")
        print(f"   Space Group: {results['space_group']}")
        print(f"   Point Group: {results['point_group']}")
        
        print(f"\nExciton States at Gamma Point:")
        print("-" * 50)
        
        for i, result in enumerate(results['results']):
            print(f"\n   State {i+1}:")
            print(f"     Energy: {result['energy']:.4f} eV")
            print(f"     Degeneracy: {result['degeneracy']}")
            print(f"     Irrep: {result['irrep']}")
        
        # Create visualization
        create_plot(egt, results)
        
        # Demonstrate LaTeX conversion
        demonstrate_latex_conversion(egt)
        
        print("\n" + "=" * 60)
        print("Analysis Complete!")
        print("=" * 60)
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        import traceback
        traceback.print_exc()

def create_plot(egt, results):
    """Create a visualization of the results."""
    
    print("\nCreating visualization...")
    
    # Extract data for plotting
    energies = []
    irrep_labels = []
    activities = []
    
    for result in results['results']:
        energies.append(result['energy'].real)
        # Extract just the irrep part (before parentheses)
        irrep_text = result['irrep'].split('(')[0].strip()
        irrep_labels.append(irrep_text)
        # Extract activity (in parentheses)
        if '(' in result['irrep']:
            activity = result['irrep'].split('(')[1].replace(')', '')
            activities.append(activity)
        else:
            activities.append('unknown')
    
    # Convert to LaTeX for plotting
    latex_labels = []
    for label in irrep_labels:
        individual_irreps = [irrep.strip() for irrep in label.split('+')]
        latex_irreps = egt.get_latex_labels(individual_irreps)
        latex_labels.append(' + '.join(latex_irreps))
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Plot energy levels
    for i, (energy, latex_label, activity) in enumerate(zip(energies, latex_labels, activities)):
        color = 'red' if 'inactive' in activity else 'blue' if 'Raman' in activity else 'green'
        ax.barh(i, energy, color=color, alpha=0.7, height=0.6)
        
        # Add LaTeX-formatted labels
        ax.text(energy + 0.01, i, latex_label, 
                va='center', ha='left', fontsize=12)
    
    ax.set_xlabel('Energy (eV)')
    ax.set_ylabel('Exciton State')
    ax.set_title(f'Exciton States - {results["space_group"]} ({results["point_group"]})')
    ax.set_yticks(range(len(energies)))
    ax.set_yticklabels([f'State {i+1}' for i in range(len(energies))])
    ax.grid(True, alpha=0.3)
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='blue', alpha=0.7, label='Raman active'),
        Patch(facecolor='red', alpha=0.7, label='Optically inactive'),
        Patch(facecolor='green', alpha=0.7, label='IR active')
    ]
    ax.legend(handles=legend_elements, loc='upper right')
    
    plt.tight_layout()
    
    # Save the plot
    plt.savefig('exciton_symmetry_analysis.png', dpi=300, bbox_inches='tight')
    print("Plot saved as 'exciton_symmetry_analysis.png'")
    
    # Show the plot
    plt.show()

def demonstrate_latex_conversion(egt):
    """Demonstrate the LaTeX conversion functionality."""
    
    print("\nLaTeX Conversion Examples:")
    print("=" * 30)
    
    test_labels = ['A1g', 'A2u', 'E1u', 'E2g', 'B1u', 'B2g']
    latex_labels = egt.get_latex_labels(test_labels)
    
    for text, latex in zip(test_labels, latex_labels):
        print(f"  {text:4s} -> {latex}")
    
    print("\nThese LaTeX labels can be used directly in matplotlib plots")
    print("with proper LaTeX rendering enabled.")

if __name__ == "__main__":
    main()