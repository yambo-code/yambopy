#!/usr/bin/env python3
"""
Example script demonstrating the use of ExcitonGroupTheory class
for analyzing the symmetry properties of exciton states.

This script shows how to:
1. Initialize the ExcitonGroupTheory class
2. Perform group theory analysis for exciton states
3. Save the results to a file

Usage:
    python exciton_group_theory_example.py

Make sure you have the required Yambo database files in the appropriate directories.
"""

import os
import sys
import numpy as np

# Add yambopy to path if needed
sys.path.insert(0, '../')

from yambopy.optical_properties import ExcitonGroupTheory

def main():
    """Main function to demonstrate ExcitonGroupTheory usage."""
    
    # Basic input parameters
    # Adjust these paths according to your calculation setup
    path = '.'  # Current directory or path to your calculation
    SAVE_dir = 'SAVE'
    BSE_dir = 'bse'  # or 'GW_BSE' depending on your setup
    LELPH_dir = 'lelph'  # Directory containing electron-phonon data
    
    # Exciton analysis parameters
    iQ = 1  # Exciton Q-point index (1-based, as in Yambo)
    nstates = 10  # Number of exciton states to analyze
    degen_thres = 0.001  # Degeneracy threshold in eV
    
    # Band range (adjust according to your calculation)
    bands_range = [1, 20]  # Example: bands 1 to 20
    
    print("Initializing ExcitonGroupTheory class...")
    print(f"Path: {path}")
    print(f"BSE directory: {BSE_dir}")
    print(f"LELPH directory: {LELPH_dir}")
    print(f"Analyzing Q-point: {iQ}")
    print(f"Number of states: {nstates}")
    print("-" * 50)
    
    try:
        # Initialize the ExcitonGroupTheory class
        egt = ExcitonGroupTheory(
            path=path,
            save=SAVE_dir,
            BSE_dir=BSE_dir,
            LELPH_dir=LELPH_dir,
            bands_range=bands_range,
            read_symm_from_ns_db_file=False  # Read symmetries from elph file
        )
        
        print("Successfully initialized ExcitonGroupTheory class!")
        print(f"Number of IBZ k-points: {egt.nibz}")
        print(f"Number of symmetry operations: {len(egt.symm_mats)}")
        print(f"Time reversal symmetry: {egt.time_rev}")
        print("-" * 50)
        
        # Perform group theory analysis
        print("Performing group theory analysis...")
        results = egt.analyze_exciton_symmetry(
            iQ=iQ,
            nstates=nstates,
            degen_thres=degen_thres
        )
        
        print("\nAnalysis completed successfully!")
        print("-" * 50)
        
        # Print summary of results
        print("SUMMARY OF RESULTS:")
        print(f"Q-point: {results['q_point']}")
        print(f"Point group: {results['point_group_label']}")
        print(f"Little group operations: {results['little_group']}")
        print(f"Number of unique energy levels: {len(results['unique_energies'])}")
        
        print("\nEnergy levels and their symmetries:")
        for i, (energy, degen, irrep) in enumerate(zip(
            results['unique_energies'],
            results['degeneracies'], 
            results['irrep_decomposition'])):
            print(f"  Level {i+1}: {energy:.4f} eV (deg={degen}) -> {irrep}")
        
        # Save results to file
        output_file = f"exciton_symmetry_Q{iQ}_analysis.txt"
        egt.save_analysis_results(results, output_file)
        print(f"\nDetailed results saved to: {output_file}")
        
    except FileNotFoundError as e:
        print(f"Error: Required database file not found: {e}")
        print("Please make sure you have the following files in the correct directories:")
        print("  - SAVE/ns.db1 (lattice database)")
        print("  - SAVE/ns.wf (wavefunction database)")
        print(f"  - {BSE_dir}/ndb.BS_diago_Q{iQ} (BSE database)")
        print(f"  - {LELPH_dir}/ndb.elph (electron-phonon database)")
        print(f"  - {LELPH_dir}/ndb.Dmats (D-matrices)")
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        print("Please check your input parameters and database files.")
        
    print("\nExample completed.")


if __name__ == "__main__":
    main()