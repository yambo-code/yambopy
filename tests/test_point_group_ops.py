#!/usr/bin/env python3
"""
Simple test script to verify the improved point_group_ops.py implementation.
"""

import numpy as np
import sys
import os

# Add the yambopy path
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

def test_basic_functionality():
    """Test basic functionality of point group operations."""
    print("Testing basic point group operations...")
    
    try:
        from yambopy.optical_properties.point_group_ops import (
            get_pg_info, decompose_rep2irrep, normalize, 
            find_symm_axis, get_point_grp
        )
        
        # Test with C2v point group matrices
        symm_mats = np.array([
            [[1, 0, 0], [0, 1, 0], [0, 0, 1]],   # E
            [[1, 0, 0], [0, -1, 0], [0, 0, -1]],  # C2
            [[-1, 0, 0], [0, 1, 0], [0, 0, -1]],  # σv
            [[-1, 0, 0], [0, -1, 0], [0, 0, 1]]   # σv'
        ])
        
        print("  Testing find_symm_axis...")
        dets, nfold, axes = find_symm_axis(symm_mats)
        print(f"    Determinants: {dets}")
        print(f"    N-fold values: {nfold}")
        print(f"    Axes shape: {axes.shape}")
        
        print("  Testing get_point_grp...")
        pg_label = get_point_grp(symm_mats)
        print(f"    Point group: {pg_label}")
        
        print("  Testing get_pg_info...")
        pg_label, classes, class_dict, char_tab, irreps = get_pg_info(symm_mats)
        print(f"    Point group: {pg_label}")
        print(f"    Classes: {classes}")
        print(f"    Irreps: {irreps}")
        print(f"    Character table shape: {char_tab.shape if char_tab is not None else 'None'}")
        
        # Test decomposition
        if char_tab is not None and len(classes) > 0:
            print("  Testing decompose_rep2irrep...")
            test_rep = np.array([2, 0, 0, 2])  # Example reducible representation
            class_orders = np.ones(len(classes), dtype=int)
            decomp = decompose_rep2irrep(test_rep, char_tab, 4, class_orders, irreps)
            print(f"    Test decomposition: {decomp}")
        
        print("  Basic functionality test: PASSED")
        return True
        
    except Exception as e:
        print(f"  Basic functionality test: FAILED - {e}")
        import traceback
        traceback.print_exc()
        return False


def test_utility_functions():
    """Test utility functions."""
    print("Testing utility functions...")
    
    try:
        from yambopy.optical_properties.point_group_ops import (
            normalize, rotation_matrix, reflection_matrix, inversion_matrix
        )
        
        # Test normalize
        vec = np.array([3, 4, 0])
        norm_vec = normalize(vec)
        expected_norm = np.array([0.6, 0.8, 0])
        if not np.allclose(norm_vec, expected_norm):
            raise ValueError("Normalize function failed")
        print("    normalize: OK")
        
        # Test rotation matrix
        axis = np.array([0, 0, 1])
        angle = np.pi / 2
        rot_mat = rotation_matrix(axis, angle)
        test_vec = np.array([1, 0, 0])
        rotated = rot_mat @ test_vec
        expected = np.array([0, 1, 0])
        if not np.allclose(rotated, expected, atol=1e-10):
            raise ValueError("Rotation matrix function failed")
        print("    rotation_matrix: OK")
        
        # Test inversion matrix
        inv_mat = inversion_matrix()
        expected_inv = -np.eye(3)
        if not np.allclose(inv_mat, expected_inv):
            raise ValueError("Inversion matrix function failed")
        print("    inversion_matrix: OK")
        
        print("  Utility functions test: PASSED")
        return True
        
    except Exception as e:
        print(f"  Utility functions test: FAILED - {e}")
        import traceback
        traceback.print_exc()
        return False

def test_from_database():
    from yambopy.optical_properties import ExcitonGroupTheory
    # Initialize the class
    egt = ExcitonGroupTheory(
        path='/mnt/lscratch/users/rreho/exc-ph-workflow/hBN-3D/',
        save='./SAVE',
        BSE_dir='./GW_BSE/bse',
        LELPH_dir='./lelph',
        bands_range=[6, 10]
    )

    # Perform group theory analysis
    results = egt.analyze_exciton_symmetry(
        iQ=1,           # Q-point index
        nstates=2,     # Number of states
        degen_thres=0.001  # Degeneracy threshold in eV
    )

    # # Save results
    # egt.save_analysis_results(results, 'exciton_symmetry.txt')
    
    # Test passed if we got here without errors
    return True
    

def main():
    """Run all tests."""
    print("Running point group operations tests...")
    print("=" * 50)
    
    tests = [
        # test_basic_functionality,
        # test_utility_functions,
        test_from_database
    ]
    
    passed = 0
    total = len(tests)
    
    for test in tests:
        if test():
            passed += 1
        print()
    
    print("=" * 50)
    print(f"Test results: {passed}/{total} tests passed")
    
    if passed == total:
        print("All tests PASSED! ✓")
        return 0
    else:
        print("Some tests FAILED! ✗")
        return 1


if __name__ == "__main__":
    sys.exit(main())
