#!/usr/bin/env python3
"""
Comprehensive test suite for ExcitonGroupTheory class.

This module contains unit tests and integration tests for the ExcitonGroupTheory
class and related point group operations.
"""

import unittest
import numpy as np
import tempfile
import os
import sys
from unittest.mock import Mock, patch, MagicMock

# Add yambopy to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

class TestPointGroupOperations(unittest.TestCase):
    """Test the point group operations module."""
    
    def setUp(self):
        """Set up test fixtures."""
        from yambopy.optical_properties.point_group_ops import (
            get_pg_info, decompose_rep2irrep, normalize, rotation_matrix,
            reflection_matrix, inversion_matrix, find_axis_angle
        )
        self.get_pg_info = get_pg_info
        self.decompose_rep2irrep = decompose_rep2irrep
        self.normalize = normalize
        self.rotation_matrix = rotation_matrix
        self.reflection_matrix = reflection_matrix
        self.inversion_matrix = inversion_matrix
        self.find_axis_angle = find_axis_angle
    
    def test_normalize(self):
        """Test vector normalization."""
        # Test 3D vector
        vec = np.array([3.0, 4.0, 0.0])
        norm_vec = self.normalize(vec)
        expected = np.array([0.6, 0.8, 0.0])
        np.testing.assert_allclose(norm_vec, expected, rtol=1e-10)
        
        # Test that normalized vector has unit length
        self.assertAlmostEqual(np.linalg.norm(norm_vec), 1.0, places=10)
        
        # Test zero vector handling
        zero_vec = np.array([0.0, 0.0, 0.0])
        with self.assertRaises(ZeroDivisionError):
            self.normalize(zero_vec)
    
    def test_rotation_matrix(self):
        """Test rotation matrix generation."""
        # Test 90-degree rotation around z-axis
        axis = np.array([0.0, 0.0, 1.0])
        angle = np.pi / 2
        rot_mat = self.rotation_matrix(axis, angle)
        
        # Test rotation of x-axis to y-axis
        x_vec = np.array([1.0, 0.0, 0.0])
        rotated = rot_mat @ x_vec
        expected = np.array([0.0, 1.0, 0.0])
        np.testing.assert_allclose(rotated, expected, atol=1e-10)
        
        # Test that rotation matrix is orthogonal
        np.testing.assert_allclose(rot_mat @ rot_mat.T, np.eye(3), atol=1e-10)
        
        # Test determinant is 1 (proper rotation)
        self.assertAlmostEqual(np.linalg.det(rot_mat), 1.0, places=10)
    
    def test_reflection_matrix(self):
        """Test reflection matrix generation."""
        # Test reflection across xy-plane (normal = z-axis)
        normal = np.array([0.0, 0.0, 1.0])
        refl_mat = self.reflection_matrix(normal)
        
        # Test reflection of z-vector
        z_vec = np.array([0.0, 0.0, 1.0])
        reflected = refl_mat @ z_vec
        expected = np.array([0.0, 0.0, -1.0])
        np.testing.assert_allclose(reflected, expected, atol=1e-10)
        
        # Test that x and y components are unchanged
        x_vec = np.array([1.0, 0.0, 0.0])
        reflected_x = refl_mat @ x_vec
        np.testing.assert_allclose(reflected_x, x_vec, atol=1e-10)
    
    def test_inversion_matrix(self):
        """Test inversion matrix generation."""
        inv_mat = self.inversion_matrix()
        expected = -np.eye(3)
        np.testing.assert_allclose(inv_mat, expected, atol=1e-10)
        
        # Test inversion of a vector
        vec = np.array([1.0, 2.0, 3.0])
        inverted = inv_mat @ vec
        expected_inv = np.array([-1.0, -2.0, -3.0])
        np.testing.assert_allclose(inverted, expected_inv, atol=1e-10)
    
    def test_find_axis_angle(self):
        """Test finding rotation axis and angle from matrix."""
        # Test 90-degree rotation around z-axis
        original_axis = np.array([0.0, 0.0, 1.0])
        original_angle = np.pi / 2
        rot_mat = self.rotation_matrix(original_axis, original_angle)
        
        found_axis, found_angle = self.find_axis_angle(rot_mat)
        
        # Check angle
        self.assertAlmostEqual(found_angle, original_angle, places=6)
        
        # Check axis (may be opposite direction)
        axis_dot = np.abs(np.dot(found_axis, original_axis))
        self.assertAlmostEqual(axis_dot, 1.0, places=6)
    
    def test_get_pg_info_c1(self):
        """Test point group identification for C1."""
        # Identity matrix only
        symm_mats = np.array([[[1, 0, 0], [0, 1, 0], [0, 0, 1]]])
        
        pg_label, classes, class_dict, char_tab, irreps = self.get_pg_info(symm_mats)
        
        self.assertEqual(pg_label, "C1")
        self.assertEqual(classes, ["E"])
        self.assertEqual(class_dict, {0: [0]})
        np.testing.assert_array_equal(char_tab, np.array([[1]]))
        self.assertEqual(irreps, ["A"])
    
    def test_get_pg_info_c2v(self):
        """Test point group identification for C2v."""
        # C2v symmetry matrices
        symm_mats = np.array([
            [[1, 0, 0], [0, 1, 0], [0, 0, 1]],   # E
            [[1, 0, 0], [0, -1, 0], [0, 0, -1]],  # C2
            [[-1, 0, 0], [0, 1, 0], [0, 0, -1]],  # σv
            [[-1, 0, 0], [0, -1, 0], [0, 0, 1]]   # σv'
        ])
        
        pg_label, classes, class_dict, char_tab, irreps = self.get_pg_info(symm_mats)
        
        self.assertEqual(pg_label, "C2v")
        self.assertEqual(len(classes), 4)
        self.assertEqual(len(irreps), 4)
        self.assertEqual(irreps, ["A1", "A2", "B1", "B2"])
        
        # Check character table dimensions
        self.assertEqual(char_tab.shape, (4, 4))
    
    def test_decompose_rep2irrep(self):
        """Test irreducible representation decomposition."""
        # Simple test with C2v character table
        char_table = np.array([
            [1, 1, 1, 1],    # A1
            [1, 1, -1, -1],  # A2
            [1, -1, 1, -1],  # B1
            [1, -1, -1, 1]   # B2
        ])
        irreps = ["A1", "A2", "B1", "B2"]
        class_orders = np.array([1, 1, 1, 1])
        pg_order = 4
        
        # Test reducible representation that decomposes to A1 + B1
        red_rep = np.array([2, 0, 2, 0])
        decomp = self.decompose_rep2irrep(red_rep, char_table, pg_order, 
                                         class_orders, irreps)
        self.assertEqual(decomp, "A1 + B1")
        
        # Test single irrep
        red_rep = np.array([1, 1, 1, 1])  # A1 irrep
        decomp = self.decompose_rep2irrep(red_rep, char_table, pg_order,
                                         class_orders, irreps)
        self.assertEqual(decomp, "A1")


class TestExcitonGroupTheoryMocked(unittest.TestCase):
    """Test ExcitonGroupTheory class with mocked dependencies."""
    
    def setUp(self):
        """Set up test fixtures with mocked dependencies."""
        # Mock all the database classes
        self.mock_lattice_db = Mock()
        self.mock_lattice_db.lat = np.eye(3)
        self.mock_lattice_db.rlat = np.eye(3)
        self.mock_lattice_db.ibz_nkpoints = 10
        self.mock_lattice_db.sym_car = np.array([np.eye(3)])
        self.mock_lattice_db.time_rev = False
        
        self.mock_wf_db = Mock()
        self.mock_wf_db.nkpoints = 100
        self.mock_wf_db.nspin = 1
        self.mock_wf_db.nspinor = 1
        self.mock_wf_db.nbands = 20
        
        self.mock_lelph_db = Mock()
        self.mock_lelph_db.kpoints = np.random.random((100, 3))
        self.mock_lelph_db.qpoints = np.random.random((10, 3))
        self.mock_lelph_db.bands = [1, 20]
        
        # Mock netCDF4 Dataset
        self.mock_dataset = Mock()
        self.mock_dataset.__enter__ = Mock(return_value=self.mock_dataset)
        self.mock_dataset.__exit__ = Mock(return_value=None)
        
        # Mock D-matrices data
        mock_dmats = np.random.random((1, 100, 1, 20, 20, 2))
        self.mock_dataset.__getitem__ = Mock(return_value=Mock(data=mock_dmats))
        self.mock_dataset.close = Mock()
    
    @patch('yambopy.optical_properties.exciton_group_theory.Dataset')
    @patch('yambopy.optical_properties.exciton_group_theory.LetzElphElectronPhononDB')
    @patch('yambopy.optical_properties.exciton_group_theory.YamboWFDB')
    @patch('yambopy.optical_properties.exciton_group_theory.YamboLatticeDB')
    @patch('yambopy.optical_properties.exciton_group_theory.build_ktree')
    def test_initialization(self, mock_build_ktree, mock_lattice_db_class, 
                           mock_wf_db_class, mock_lelph_db_class, mock_dataset_class):
        """Test ExcitonGroupTheory initialization."""
        # Set up mocks
        mock_lattice_db_class.from_db_file.return_value = self.mock_lattice_db
        mock_wf_db_class.return_value = self.mock_wf_db
        mock_lelph_db_class.return_value = self.mock_lelph_db
        mock_dataset_class.return_value = self.mock_dataset
        mock_build_ktree.return_value = Mock()
        
        # Create temporary directory structure
        with tempfile.TemporaryDirectory() as temp_dir:
            save_dir = os.path.join(temp_dir, 'SAVE')
            bse_dir = os.path.join(temp_dir, 'bse')
            lelph_dir = os.path.join(temp_dir, 'lelph')
            
            os.makedirs(save_dir)
            os.makedirs(bse_dir)
            os.makedirs(lelph_dir)
            
            # Create dummy files
            open(os.path.join(save_dir, 'ns.db1'), 'w').close()
            open(os.path.join(save_dir, 'ns.wf'), 'w').close()
            open(os.path.join(lelph_dir, 'ndb.elph'), 'w').close()
            open(os.path.join(lelph_dir, 'ndb.Dmats'), 'w').close()
            
            # Test initialization
            from yambopy.optical_properties.exciton_group_theory import ExcitonGroupTheory
            
            egt = ExcitonGroupTheory(
                path=temp_dir,
                save='SAVE',
                BSE_dir='bse',
                LELPH_dir='lelph',
                bands_range=[1, 10]
            )
            
            # Check that attributes are set correctly
            self.assertEqual(egt.path, temp_dir)
            self.assertEqual(egt.SAVE_dir, save_dir)
            self.assertEqual(egt.BSE_dir, bse_dir)
            self.assertEqual(egt.LELPH_dir, lelph_dir)
    
    def test_file_not_found_error(self):
        """Test that appropriate errors are raised for missing files."""
        from yambopy.optical_properties.exciton_group_theory import ExcitonGroupTheory
        
        with tempfile.TemporaryDirectory() as temp_dir:
            # Try to initialize without required files
            with self.assertRaises(IOError):
                ExcitonGroupTheory(path=temp_dir)


class TestExcitonGroupTheoryIntegration(unittest.TestCase):
    """Integration tests for ExcitonGroupTheory (requires actual data files)."""
    
    def setUp(self):
        """Set up integration test fixtures."""
        # These tests require actual Yambo database files
        # Skip if not available
        self.test_data_available = False
        
        # Check for test data directory
        test_data_dir = os.path.join(os.path.dirname(__file__), 'test_data')
        if os.path.exists(test_data_dir):
            required_files = [
                'SAVE/ns.db1',
                'SAVE/ns.wf',
                'bse/ndb.BS_diago_Q1',
                'lelph/ndb.elph',
                'lelph/ndb.Dmats'
            ]
            
            if all(os.path.exists(os.path.join(test_data_dir, f)) for f in required_files):
                self.test_data_available = True
                self.test_data_dir = test_data_dir
    
    @unittest.skipUnless(False, "Integration tests require actual Yambo data files")
    def test_full_analysis_workflow(self):
        """Test complete analysis workflow with real data."""
        if not self.test_data_available:
            self.skipTest("Test data not available")
        
        from yambopy.optical_properties.exciton_group_theory import ExcitonGroupTheory
        
        # Initialize with test data
        egt = ExcitonGroupTheory(
            path=self.test_data_dir,
            save='SAVE',
            BSE_dir='bse',
            LELPH_dir='lelph',
            bands_range=[1, 10]
        )
        
        # Perform analysis
        results = egt.analyze_exciton_symmetry(
            iQ=1,
            nstates=5,
            degen_thres=0.001
        )
        
        # Check results structure
        required_keys = [
            'q_point', 'little_group', 'point_group_label',
            'unique_energies', 'degeneracies', 'irrep_decomposition'
        ]
        
        for key in required_keys:
            self.assertIn(key, results)
        
        # Check data types and shapes
        self.assertIsInstance(results['q_point'], np.ndarray)
        self.assertEqual(len(results['q_point']), 3)
        self.assertIsInstance(results['point_group_label'], str)
        self.assertIsInstance(results['unique_energies'], np.ndarray)
        self.assertIsInstance(results['degeneracies'], np.ndarray)
        self.assertIsInstance(results['irrep_decomposition'], list)
        
        # Check that energies are positive
        self.assertTrue(np.all(results['unique_energies'] > 0))
        
        # Check that degeneracies are positive integers
        self.assertTrue(np.all(results['degeneracies'] > 0))
        self.assertTrue(np.all(results['degeneracies'] == results['degeneracies'].astype(int)))


class TestUtilityFunctions(unittest.TestCase):
    """Test utility functions and helper methods."""
    
    def test_matrix_properties(self):
        """Test mathematical properties of generated matrices."""
        from yambopy.optical_properties.point_group_ops import (
            rotation_matrix, reflection_matrix, inversion_matrix
        )
        
        # Test rotation matrix properties
        axis = np.array([1.0, 1.0, 1.0]) / np.sqrt(3)  # Normalized
        angle = np.pi / 3
        rot_mat = rotation_matrix(axis, angle)
        
        # Should be orthogonal
        np.testing.assert_allclose(rot_mat @ rot_mat.T, np.eye(3), atol=1e-10)
        
        # Should have determinant 1
        self.assertAlmostEqual(np.linalg.det(rot_mat), 1.0, places=10)
        
        # Test reflection matrix properties
        normal = np.array([1.0, 0.0, 0.0])
        refl_mat = reflection_matrix(normal)
        
        # Should be symmetric
        np.testing.assert_allclose(refl_mat, refl_mat.T, atol=1e-10)
        
        # Should have determinant -1
        self.assertAlmostEqual(np.linalg.det(refl_mat), -1.0, places=10)
        
        # Should be its own inverse
        np.testing.assert_allclose(refl_mat @ refl_mat, np.eye(3), atol=1e-10)
        
        # Test inversion matrix
        inv_mat = inversion_matrix()
        self.assertAlmostEqual(np.linalg.det(inv_mat), -1.0, places=10)
        np.testing.assert_allclose(inv_mat @ inv_mat, np.eye(3), atol=1e-10)


class TestErrorHandling(unittest.TestCase):
    """Test error handling and edge cases."""
    
    def test_invalid_parameters(self):
        """Test handling of invalid parameters."""
        from yambopy.optical_properties.point_group_ops import decompose_rep2irrep
        
        # Test with None character table
        result = decompose_rep2irrep(
            np.array([1, 1, 1]), None, 4, np.array([1, 1, 1]), []
        )
        self.assertEqual(result, "Analysis not available")
        
        # Test with empty irreps list
        char_table = np.array([[1, 1], [1, -1]])
        result = decompose_rep2irrep(
            np.array([2, 0]), char_table, 2, np.array([1, 1]), []
        )
        self.assertEqual(result, "Analysis not available")
    
    def test_edge_cases(self):
        """Test edge cases and boundary conditions."""
        from yambopy.optical_properties.point_group_ops import get_pg_info
        
        # Test with large number of symmetries (unknown group)
        large_symm_mats = np.random.random((20, 3, 3))
        
        with self.assertWarns(UserWarning):
            pg_label, classes, class_dict, char_tab, irreps = get_pg_info(large_symm_mats)
        
        self.assertTrue(pg_label.startswith("Unknown_"))
        self.assertEqual(len(classes), 20)
        self.assertEqual(len(irreps), 20)


def create_test_suite():
    """Create a comprehensive test suite."""
    suite = unittest.TestSuite()
    
    # Add all test classes
    test_classes = [
        TestPointGroupOperations,
        TestExcitonGroupTheoryMocked,
        TestExcitonGroupTheoryIntegration,
        TestUtilityFunctions,
        TestErrorHandling
    ]
    
    for test_class in test_classes:
        tests = unittest.TestLoader().loadTestsFromTestCase(test_class)
        suite.addTests(tests)
    
    return suite


def run_tests(verbosity=2):
    """Run all tests with specified verbosity."""
    suite = create_test_suite()
    runner = unittest.TextTestRunner(verbosity=verbosity)
    result = runner.run(suite)
    
    return result.wasSuccessful()


if __name__ == '__main__':
    # Run tests when script is executed directly
    print("Running ExcitonGroupTheory test suite...")
    print("=" * 60)
    
    success = run_tests(verbosity=2)
    
    print("\n" + "=" * 60)
    if success:
        print("All tests PASSED! ✓")
        sys.exit(0)
    else:
        print("Some tests FAILED! ✗")
        sys.exit(1)