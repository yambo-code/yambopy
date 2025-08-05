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
Point group operations and character table analysis for exciton symmetry analysis.
This is a simplified version of the original point_group_ops.py module.
"""
"""
MIT License

Copyright (c) 2025 Muralidhar Nalabothula
Copyright (c) 2023 Stephen M. Goodlett, Nathaniel L. Kitzmiller

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""
# Adapated and modified from MolSym package
# MolSym: A Python package for handling symmetry in molecular quantum chemistry
## Cite the following paper :https://doi.org/10.1063/5.0216738


import numpy as np
import warnings

def get_pg_info(symm_mats):
    """
    Get point group information from symmetry matrices.
    
    This is a simplified implementation that provides basic point group
    identification and character table information.
    
    Parameters
    ----------
    symm_mats : numpy.ndarray
        Array of symmetry matrices.
        
    Returns
    -------
    pg_label : str
        Point group label.
    classes : list
        List of symmetry classes.
    class_dict : dict
        Dictionary mapping class indices to symmetry operations.
    char_tab : numpy.ndarray
        Character table.
    irreps : list
        List of irreducible representations.
    """
    n_symm = len(symm_mats)
    
    # Simple point group identification based on number of symmetries
    # This is a very basic implementation
    if n_symm == 1:
        pg_label = "C1"
        classes = ["E"]
        class_dict = {0: [0]}
        char_tab = np.array([[1]])
        irreps = ["A"]
    elif n_symm == 2:
        pg_label = "Ci"
        classes = ["E", "i"]
        class_dict = {0: [0], 1: [1]}
        char_tab = np.array([[1, 1], [1, -1]])
        irreps = ["Ag", "Au"]
    elif n_symm == 4:
        pg_label = "C2v"
        classes = ["E", "C2", "σv", "σv'"]
        class_dict = {0: [0], 1: [1], 2: [2], 3: [3]}
        char_tab = np.array([[1, 1, 1, 1], 
                            [1, 1, -1, -1], 
                            [1, -1, 1, -1], 
                            [1, -1, -1, 1]])
        irreps = ["A1", "A2", "B1", "B2"]
    elif n_symm == 8:
        pg_label = "D2h"
        classes = ["E", "C2z", "C2y", "C2x", "i", "σxy", "σxz", "σyz"]
        class_dict = {i: [i] for i in range(8)}
        char_tab = np.array([[1, 1, 1, 1, 1, 1, 1, 1],
                            [1, 1, -1, -1, 1, 1, -1, -1],
                            [1, -1, 1, -1, 1, -1, 1, -1],
                            [1, -1, -1, 1, 1, -1, -1, 1],
                            [1, 1, 1, 1, -1, -1, -1, -1],
                            [1, 1, -1, -1, -1, -1, 1, 1],
                            [1, -1, 1, -1, -1, 1, -1, 1],
                            [1, -1, -1, 1, -1, 1, 1, -1]])
        irreps = ["Ag", "B1g", "B2g", "B3g", "Au", "B1u", "B2u", "B3u"]
    else:
        # Default case for unknown point groups
        pg_label = f"Unknown_{n_symm}"
        classes = [f"Op{i}" for i in range(n_symm)]
        class_dict = {i: [i] for i in range(n_symm)}
        char_tab = np.eye(n_symm)
        irreps = [f"Irrep{i}" for i in range(n_symm)]
        
        warnings.warn(f"Point group with {n_symm} operations not recognized. "
                     "Using default classification.")
    
    return pg_label, classes, class_dict, char_tab, irreps


def decompose_rep2irrep(red_rep, char_table, pg_order, class_order, irreps):
    """
    Decompose a reducible representation into irreducible representations.
    
    Parameters
    ----------
    red_rep : numpy.ndarray
        Characters of the reducible representation.
    char_table : numpy.ndarray
        Character table of the point group.
    pg_order : int
        Order of the point group.
    class_order : numpy.ndarray
        Order of each class.
    irreps : list
        List of irreducible representation labels.
        
    Returns
    -------
    decomposition : str
        String representation of the decomposition.
    """
    if char_table is None or len(irreps) == 0:
        return "Analysis not available"
    
    n_irreps = len(irreps)
    coefficients = []
    
    for i in range(n_irreps):
        # Calculate coefficient using the reduction formula
        coeff = 0
        for j, char_red in enumerate(red_rep):
            if j < len(char_table[i]) and j < len(class_order):
                coeff += class_order[j] * char_red.real * char_table[i][j]
        coeff = int(round(coeff / pg_order))
        coefficients.append(coeff)
    
    # Build decomposition string
    decomposition_parts = []
    for i, coeff in enumerate(coefficients):
        if coeff > 0:
            if coeff == 1:
                decomposition_parts.append(irreps[i])
            else:
                decomposition_parts.append(f"{coeff}{irreps[i]}")
    
    if decomposition_parts:
        return " + ".join(decomposition_parts)
    else:
        return "0"



def normalize(vec):
    norm = np.linalg.norm(vec)
    if norm == 0:
        raise ZeroDivisionError("Cannot normalize a zero vector.")
    return vec / norm

def rotation_matrix(axis, theta):
    """
    Create a rotation matrix for rotation around an axis by angle theta.
    
    Parameters
    ----------
    axis : numpy.ndarray
        Rotation axis (3D vector).
    theta : float
        Rotation angle in radians.
        
    Returns
    -------
    numpy.ndarray
        3x3 rotation matrix.
    """
    axis = normalize(axis)
    a = np.cos(theta / 2.0)
    b, c, d = -axis * np.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])


def reflection_matrix(axis):
    """
    Create a reflection matrix for reflection across a plane with normal axis.
    
    Parameters
    ----------
    axis : numpy.ndarray
        Normal to the reflection plane (3D vector).
        
    Returns
    -------
    numpy.ndarray
        3x3 reflection matrix.
    """
    axis = normalize(axis)
    return np.eye(3) - 2 * np.outer(axis, axis)


def inversion_matrix():
    """Create an inversion matrix."""
    return -np.eye(3)


def find_axis_angle(Rmat):
    """
    Find the rotation axis and angle from a rotation matrix.
    
    Parameters
    ----------
    Rmat : numpy.ndarray
        3x3 rotation matrix.
        
    Returns
    -------
    axis : numpy.ndarray
        Rotation axis.
    angle : float
        Rotation angle in radians.
    """
    # Calculate the rotation angle
    trace = np.trace(Rmat)
    angle = np.arccos((trace - 1) / 2)
    
    if np.abs(angle) < 1e-6:  # Identity matrix
        return np.array([0, 0, 1]), 0.0
    elif np.abs(angle - np.pi) < 1e-6:  # 180-degree rotation
        # Find the eigenvector with eigenvalue 1
        eigenvals, eigenvecs = np.linalg.eig(Rmat)
        idx = np.argmin(np.abs(eigenvals - 1))
        axis = eigenvecs[:, idx].real
        return normalize(axis), angle
    else:
        # General case
        axis = np.array([Rmat[2, 1] - Rmat[1, 2],
                        Rmat[0, 2] - Rmat[2, 0],
                        Rmat[1, 0] - Rmat[0, 1]])
        axis = normalize(axis)
        return axis, angle


def point_group_classes(point_group, tol=1e-4):
    """
    Determine the conjugacy classes of a point group.
    
    Parameters
    ----------
    point_group : list
        List of symmetry matrices.
    tol : float
        Tolerance for matrix comparison.
        
    Returns
    -------
    classes : list
        List of conjugacy classes.
    """
    n = len(point_group)
    classes = []
    assigned = [False] * n
    
    for i in range(n):
        if assigned[i]:
            continue
            
        # Start a new class with element i
        current_class = [i]
        assigned[i] = True
        
        # Find all elements conjugate to element i
        for j in range(n):
            if assigned[j]:
                continue
                
            # Check if j is conjugate to i
            for k in range(n):
                conjugate = point_group[k] @ point_group[i] @ np.linalg.inv(point_group[k])
                if np.allclose(conjugate, point_group[j], atol=tol):
                    current_class.append(j)
                    assigned[j] = True
                    break
        
        classes.append(current_class)
    
    return classes