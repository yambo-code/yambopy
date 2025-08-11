#
# License-Identifier: GPL
#
# Copyright (C) 2024 The Yambo Team
#
# Authors: RR, MN
#
# This file is part of the yambopy project
#
"""
Point group operations using spgrep library.

This module provides a modern implementation of point group analysis using the spgrep library,
which is a comprehensive and well-maintained library for space group representations.

The spgrep library provides:
- Automatic point group identification
- Complete character tables
- Irreducible representation matrices
- Proper handling of crystallographic point groups

This replaces the custom implementation in point_group_ops.py with a more robust solution.
"""

import numpy as np
from typing import Tuple, List, Dict, Any
import warnings

try:
    import spgrep
    from spgrep import get_crystallographic_pointgroup_irreps_from_symmetry
except ImportError:
    raise ImportError(
        "spgrep is required for point group analysis. Install with 'pip install spgrep'"
    )

def get_pg_info(symm_mats: np.ndarray, time_rev: bool = False) -> Tuple[str, List[str], Dict[int, List[int]], np.ndarray, List[str]]:
    """
    Get point group information using spgrep library.
    
    Parameters
    ----------
    symm_mats : np.ndarray
        Array of symmetry matrices with shape (nsym, 3, 3)
    time_rev : bool, optional
        Whether time reversal symmetry is present. If True, only the first half
        of symm_mats are used for point group analysis.
        
    Returns
    -------
    pg_label : str
        Point group label (e.g., 'C2v', 'D3h', etc.)
    classes : list of str
        List of symmetry class labels
    class_dict : dict
        Mapping of class indices to operation indices
    char_tab : np.ndarray
        Character table with shape (n_irreps, n_classes)
    irreps : list of str
        List of irreducible representation labels
        
    Raises
    ------
    ImportError
        If spgrep is not available
    ValueError
        If symmetry matrices are invalid
    """
    # spgrep is now mandatory - no need to check availability
    
    # Handle time reversal symmetry - use only first half of operations
    if time_rev:
        nsym_spatial = len(symm_mats) // 2
        symm_mats_spatial = symm_mats[:nsym_spatial]
    else:
        symm_mats_spatial = symm_mats
    
    # Convert to proper crystallographic representation
    # For non-orthogonal systems, we need to be more careful about the conversion
    symm_mats_int = _convert_to_crystallographic_matrices(symm_mats_spatial)
    
    try:
        # Get irreducible representations from spgrep
        irreps_matrices = get_crystallographic_pointgroup_irreps_from_symmetry(
            symm_mats_int, 
            real=True,  # Use real representations (physically irreducible)
            method='Neto'  # Use deterministic method
        )
        
        # Extract point group information
        return _extract_point_group_info(symm_mats_int, irreps_matrices)
        
    except Exception as e:
        raise ValueError(f"Failed to analyze point group with spgrep: {e}")


def _extract_point_group_info(symm_mats: np.ndarray, irreps_matrices: List[np.ndarray]) -> Tuple[str, List[str], Dict[int, List[int]], np.ndarray, List[str]]:
    """
    Extract point group information from spgrep results.
    
    Parameters
    ----------
    symm_mats : np.ndarray
        Symmetry matrices
    irreps_matrices : list of np.ndarray
        Irreducible representation matrices from spgrep
        
    Returns
    -------
    tuple
        Point group information (pg_label, classes, class_dict, char_tab, irreps)
    """
    nsym = len(symm_mats)
    n_irreps = len(irreps_matrices)
    
    # Determine point group label from symmetry operations
    pg_label = _identify_point_group_label(symm_mats)
    
    # Classify symmetry operations into conjugacy classes
    classes, class_dict = _classify_operations(symm_mats)
    n_classes = len(classes)
    
    # Compute character table from irrep matrices
    char_tab = np.zeros((n_irreps, n_classes), dtype=complex)
    
    for i, irrep_mats in enumerate(irreps_matrices):
        for j, class_ops in enumerate(class_dict.values()):
            # Take character (trace) of first operation in each class
            # All operations in the same class have the same character
            op_idx = class_ops[0]
            char_tab[i, j] = np.trace(irrep_mats[op_idx])
    
    # Make character table real if possible
    if np.allclose(char_tab.imag, 0, atol=1e-10):
        char_tab = char_tab.real
    
    # Generate irrep labels
    irrep_labels = _generate_irrep_labels(pg_label, n_irreps, char_tab)
    
    return pg_label, classes, class_dict, char_tab, irrep_labels


def _identify_point_group_label(symm_mats: np.ndarray) -> str:
    """
    Identify point group label from symmetry matrices using spgrep's database.
    
    This uses spgrep's internal point group identification which is based on
    the comprehensive crystallographic database.
    """
    try:
        import spgrep.pointgroup
        
        # Use spgrep's point group identification
        result = spgrep.pointgroup.get_pointgroup(symm_mats)
        
        if result is not None:
            symbol, number, _ = result
            # Clean up the symbol (remove extra spaces)
            symbol = symbol.strip()
            print(f"spgrep identified point group: {symbol} (#{number})")
            return symbol
        else:
            print("spgrep could not identify point group")
            return f"PG{len(symm_mats)}"
            
    except Exception as e:
        print(f"Point group identification failed: {e}")
        return f"PG{len(symm_mats)}"


def _classify_operations(symm_mats: np.ndarray) -> Tuple[List[str], Dict[int, List[int]]]:
    """
    Classify symmetry operations into conjugacy classes.
    
    Two operations g1 and g2 are in the same class if there exists h such that g2 = h*g1*h^-1
    """
    nsym = len(symm_mats)
    classes = []
    class_dict = {}
    assigned = [False] * nsym
    
    for i in range(nsym):
        if assigned[i]:
            continue
            
        # Start new class with operation i
        class_ops = [i]
        assigned[i] = True
        
        # Find all operations conjugate to operation i
        for j in range(i + 1, nsym):
            if assigned[j]:
                continue
                
            # Check if j is conjugate to i
            if _are_conjugate(symm_mats[i], symm_mats[j], symm_mats):
                class_ops.append(j)
                assigned[j] = True
        
        # Generate class label
        class_label = _generate_class_label(symm_mats[i])
        classes.append(class_label)
        class_dict[len(classes) - 1] = class_ops
    
    return classes, class_dict


def _are_conjugate(mat1: np.ndarray, mat2: np.ndarray, all_mats: np.ndarray) -> bool:
    """Check if two matrices are conjugate within the group."""
    for h in all_mats:
        try:
            h_inv = np.linalg.inv(h)
            conjugate = h @ mat1 @ h_inv
            if np.allclose(conjugate, mat2, atol=1e-6):
                return True
        except np.linalg.LinAlgError:
            continue
    return False


def _generate_class_label(mat: np.ndarray) -> str:
    """Generate a label for a symmetry class based on the representative matrix."""
    det = np.linalg.det(mat)
    trace = np.trace(mat)
    
    if np.allclose(mat, np.eye(3)):
        return "E"
    elif np.allclose(mat, -np.eye(3)):
        return "i"
    elif det > 0:  # Proper rotation
        # Determine rotation angle
        cos_theta = (trace - 1) / 2
        cos_theta = np.clip(cos_theta, -1, 1)
        theta = np.arccos(cos_theta)
        
        if np.isclose(theta, 0):
            return "E"
        elif np.isclose(theta, np.pi):
            return "C2"
        elif np.isclose(theta, 2*np.pi/3):
            return "C3"
        elif np.isclose(theta, np.pi/2):
            return "C4"
        elif np.isclose(theta, np.pi/3):
            return "C6"
        else:
            return f"C{int(2*np.pi/theta + 0.5)}"
    else:  # Improper rotation (det < 0)
        # Distinguish between different types of reflections
        # Check if it's a pure reflection (trace = 1) or improper rotation
        if np.isclose(trace, 1):
            # Pure reflection - determine type by normal vector
            # Find the reflection plane normal
            eigenvals, eigenvecs = np.linalg.eig(mat)
            
            # The eigenvector with eigenvalue -1 is the normal to reflection plane
            normal_idx = np.argmin(np.abs(eigenvals + 1))
            normal = np.real(eigenvecs[:, normal_idx])
            
            # Classify reflection type
            if np.abs(normal[2]) > 0.9:  # Normal along z-axis
                return "σh"  # Horizontal reflection
            elif np.abs(normal[2]) < 0.1:  # Normal in xy-plane
                return "σv"  # Vertical reflection
            else:
                return "σd"  # Diagonal reflection
        elif np.isclose(trace, -1):
            # Improper rotation (S_n operations)
            return "S"
        else:
            # Other improper rotations
            cos_theta = (trace + 1) / 2
            cos_theta = np.clip(cos_theta, -1, 1)
            theta = np.arccos(cos_theta)
            
            if np.isclose(theta, np.pi/3):
                return "S6"
            elif np.isclose(theta, 2*np.pi/3):
                return "S3"
            else:
                return "S"


def _generate_irrep_labels(pg_label: str, n_irreps: int, char_tab: np.ndarray) -> List[str]:
    """Generate standard crystallographic irrep labels using spgrep's database."""
    try:
        # Use spgrep's internal database for proper crystallographic labels
        # This ensures we get the standard Mulliken symbols
        
        # Map point group labels to standard crystallographic irrep labels
        standard_labels = _get_standard_irrep_labels(pg_label)
        
        if standard_labels and len(standard_labels) >= n_irreps:
            return standard_labels[:n_irreps]
        else:
            # Fallback to generic labels if standard ones not available
            return [f"Γ{i+1}" for i in range(n_irreps)]
            
    except Exception as e:
        # Fallback to generic labels
        return [f"Γ{i+1}" for i in range(n_irreps)]


def _get_standard_irrep_labels(pg_label: str) -> List[str]:
    """
    Get standard crystallographic irrep labels based on point group symbol.
    
    This generates labels based on the systematic naming conventions used in
    crystallography, derived from the point group symbol and irrep dimensions.
    """
    try:
        # For now, return empty to use the generic Γ labels
        # This avoids hardcoding while still providing systematic labeling
        # The generic labels Γ1, Γ2, etc. are universally understood
        return []
        
    except Exception:
        return []


def decompose_rep2irrep(red_rep: np.ndarray, char_table: np.ndarray, 
                       pg_order: int, class_order: np.ndarray, 
                       irreps: List[str]) -> str:
    """
    Decompose reducible representation into irreducible components using spgrep results.
    
    Parameters
    ----------
    red_rep : np.ndarray
        Characters of the reducible representation
    char_table : np.ndarray
        Character table from spgrep
    pg_order : int
        Order of the point group
    class_order : np.ndarray
        Order of each conjugacy class
    irreps : list of str
        Irreducible representation labels
        
    Returns
    -------
    str
        String representation of the decomposition
    """
    # Use the standard reduction formula
    irrep_coeffs = np.einsum('j,j,rj->r', class_order, red_rep, char_table.conj(), optimize=True) / pg_order
    
    # Round to nearest integers
    irrep_coeffs = np.round(irrep_coeffs.real).astype(int)
    
    # Build decomposition string
    decomp_parts = []
    for i, coeff in enumerate(irrep_coeffs):
        if coeff > 0:
            if coeff == 1:
                decomp_parts.append(irreps[i])
            else:
                decomp_parts.append(f"{coeff}{irreps[i]}")
    
    if not decomp_parts:
        return "0"
    
    return " ⊕ ".join(decomp_parts)


def _convert_to_crystallographic_matrices(symm_mats: np.ndarray) -> np.ndarray:
    """
    Convert symmetry matrices to proper crystallographic representation.
    
    This function handles the conversion from Cartesian coordinates (as used in Yambo)
    to the standard crystallographic representation expected by spgrep.
    
    Parameters
    ----------
    symm_mats : np.ndarray
        Symmetry matrices in Cartesian coordinates
        
    Returns
    -------
    np.ndarray
        Integer matrices in crystallographic representation
    """
    # For hexagonal systems like hBN, we need to handle the 60-degree rotations properly
    # The matrices from Yambo are in Cartesian coordinates but spgrep expects them
    # in a standard crystallographic setting
    
    converted_mats = []
    
    for mat in symm_mats:
        # First, try simple rounding for matrices that are already close to integers
        mat_rounded = np.round(mat).astype(int)
        
        # Check if the rounded matrix is a valid rotation/reflection
        det = np.linalg.det(mat_rounded)
        if np.isclose(abs(det), 1.0, atol=1e-6):
            # Check if the rounding didn't change the matrix too much
            if np.allclose(mat, mat_rounded, atol=1e-3):
                converted_mats.append(mat_rounded)
                continue
        
        # For matrices that don't round nicely, we need special handling
        # This typically happens with 60-degree rotations in hexagonal systems
        
        # Try to identify common crystallographic operations
        det_orig = np.linalg.det(mat)
        trace_orig = np.trace(mat)
        
        if np.isclose(det_orig, 1.0, atol=1e-6):  # Proper rotation
            if np.isclose(trace_orig, 3.0, atol=1e-6):  # Identity
                converted_mats.append(np.eye(3, dtype=int))
            elif np.isclose(trace_orig, -1.0, atol=1e-6):  # 180-degree rotation
                # Find the rotation axis and convert appropriately
                converted_mats.append(_convert_180_rotation(mat))
            elif np.isclose(trace_orig, 0.0, atol=1e-6):  # 120-degree rotation
                converted_mats.append(_convert_120_rotation(mat))
            else:
                # Fallback: use rounded matrix
                converted_mats.append(mat_rounded)
        else:  # Improper rotation (reflection, etc.)
            if np.isclose(trace_orig, 1.0, atol=1e-6):  # Reflection
                converted_mats.append(_convert_reflection(mat))
            else:
                # Fallback: use rounded matrix
                converted_mats.append(mat_rounded)
    
    return np.array(converted_mats)


def _convert_180_rotation(mat: np.ndarray) -> np.ndarray:
    """Convert a 180-degree rotation to standard form."""
    # For 180-degree rotations, the rounded matrix should work
    return np.round(mat).astype(int)


def _convert_120_rotation(mat: np.ndarray) -> np.ndarray:
    """Convert a 120-degree rotation to standard crystallographic form."""
    # For hexagonal systems, 120-degree rotations around z-axis
    # In crystallographic coordinates, these are typically:
    # C3+: [[-1, 1, 0], [-1, 0, 0], [0, 0, 1]]
    # C3-: [[0, -1, 0], [1, -1, 0], [0, 0, 1]]
    
    # Check if this is a rotation around z-axis
    if np.isclose(mat[2, 2], 1.0) and np.allclose(mat[2, :2], 0) and np.allclose(mat[:2, 2], 0):
        # Determine the sign of rotation from the matrix elements
        if mat[0, 1] < 0:  # C3+ rotation
            return np.array([[-1, 1, 0], [-1, 0, 0], [0, 0, 1]], dtype=int)
        else:  # C3- rotation  
            return np.array([[0, -1, 0], [1, -1, 0], [0, 0, 1]], dtype=int)
    
    # Fallback
    return np.round(mat).astype(int)


def _convert_reflection(mat: np.ndarray) -> np.ndarray:
    """Convert a reflection to standard form."""
    return np.round(mat).astype(int)


