"""
MIT License

Copyright (c) 2025 YamboPy team 
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
# Adapted and modified from MolSym package
# MolSym: A Python package for handling symmetry in molecular quantum chemistry
## Cite the following paper :https://doi.org/10.1063/5.0216738

import numpy as np
import re
from numpy.linalg import matrix_power
from dataclasses import dataclass
from scipy.spatial import KDTree
from yambopy.tools.citations import citation

global_tol = 1e-6


@citation("S. M. Goodlett et al. J. Chem. Phys. 161, 024107 (2024)")
def get_pg_info(symm_mats):
    """
    Given list of symmetries that form a crystallographic point 
    group, return point-group label, classes, character table, irrep_labels
    """
    order = len(symm_mats)
    pg_label = get_point_grp(symm_mats)
    pg_symels = pg_to_symels(pg_label)
    pg_sym_mats = []
    for i in pg_symels:
        pg_sym_mats.append(i.rrep)
    pg_sym_mats = np.array(pg_sym_mats)
    assert len(pg_sym_mats) == order
    transfom_mat = transform_matrix(symm_mats, pg_sym_mats)
    symm_mats_transformed = transfom_mat[
        None, :, :] @ symm_mats @ transfom_mat.T[None, :, :]
    assert np.imag(symm_mats_transformed).max() < 1e-4
    symm_mats_transformed = symm_mats_transformed.real
    sym_tree = KDTree(pg_sym_mats.reshape(order, -1))
    distance, idx = sym_tree.query(symm_mats_transformed.reshape(order, -1),
                                   k=1)
    assert np.max(distance) < 1e-4, np.max(distance)
    assert len(np.unique(idx)) == order
    ctab = pg_to_chartab(pg_label)
    classes = ctab.classes
    irreps = ctab.irreps
    char_tab = ctab.characters
    pg_class_map = np.array(generate_symel_to_class_map(pg_symels, ctab))
    class_map = pg_class_map[idx]
    class_dict = dict()
    for i, n in enumerate(class_map):
        class_dict.setdefault(n, []).append(i)
    class_dict = {n: rep for n, rep in class_dict.items() if len(rep) > 0}
    return pg_label, classes, class_dict, char_tab, irreps


def decompose_rep2irrep(red_rep, char_table, pg_order, class_order,
                        irrep_labels, tol=1e-3):
    irrep_coeff = np.einsum(
        'j,j,rj->r', class_order, red_rep, char_table, optimize=True) / pg_order
    assert np.abs(irrep_coeff.imag).max() < tol, print(
        np.abs(irrep_coeff.imag).max())
    irrep_coeff = irrep_coeff.real
    assert np.abs(irrep_coeff-np.rint(irrep_coeff)).max() < tol,\
                np.abs(irrep_coeff-np.rint(irrep_coeff)).max()
    irrep_coeff = np.rint(irrep_coeff).astype(int)
    rep_string = ''
    for i in range(len(irrep_labels)):
        if irrep_coeff[i] == 0:
            continue
        num_str = ''
        if irrep_coeff[i] > 1:
            num_str += str(irrep_coeff[i])
        rep_string = rep_string + num_str + irrep_labels[i] + ' + '
    return rep_string.strip()[:-1]


## Helpers
################################################################################
def normalize(a):
    """
    Normalize vector a to unit length
    :param a: Vector of arbitrary magnitude
    :type a: NumPy array of shape (n,)
    :return: Normalized vector or None if the magnitude of a is less than the global tolerance
    :rtype: NumPy array of shape (n,) or None
    """
    n = np.linalg.norm(a)
    if n <= global_tol:
        return a
    return a / np.linalg.norm(a)


def rotation_matrix(axis, theta):
    """
    Rotation matrix about an axis by theta in radians.

    :param axis: Cartesian vector defining rotation axis
    :param theta: Angle of rotation in radians
    :type axis: NumPy array of shape (3,)
    :type theta: float
    :return: Matrix defining rotation on column vector
    :rtype: NumPy array of shape (3,3)
    """
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)
    # NOT NORMALIZING AXIS!!!
    M = cos_t * np.eye(3)
    M += sin_t * np.cross(np.eye(3), axis)
    M += (1 - cos_t) * np.outer(axis, axis)
    return M


def reflection_matrix(axis):
    """
    Reflection matrix about a plane defined by its normal vector.

    :param axis: Cartesian vector defining the plane normal vector
    :type axis: NumPy array of shape (3,)
    :return: Matrix defining reflection on column vector
    :rtype: NumPy array of shape (3,3)
    """
    M = np.zeros((3, 3))
    for i in range(3):
        for j in range(i, 3):
            if i == j:
                M[i, i] = 1 - 2 * (axis[i]**2)
            else:
                M[i, j] = -2 * axis[i] * axis[j]
                M[j, i] = M[i, j]
    return M


def inversion_matrix():
    """
    Cartesian inversion matrix.

    :return: Matrix defining inversion
    :rtype: NumPy array of shape(3,3)
    """
    return -1 * np.eye(3)


def Cn(axis, n):
    """
    Wrapper around rotation_matrix for producing a C_n rotation about axis.

    :param axis: Cartesian vector defining rotation axis
    :param n: Defines rotation angle by theta = 2 pi / n
    :type axis: NumPy array of shape (3,)
    :type n: int
    :return: Matrix defining proper rotation on column vector
    :rtype: NumPy array of shape (3,3)
    """
    theta = 2 * np.pi / n
    return rotation_matrix(axis, theta)


def Sn(axis, n):
    """
    Improper rotation S_n about an axis.
    
    :param axis: Cartesian vector defining rotation axis
    :param n: Defines rotation angle by theta = 2 pi / n
    :type axis: NumPy array of shape (3,)
    :type n: int
    :return: Matrix defining improper rotation on column vector
    :rtype: NumPy array of shape (3,3)
    """
    return np.dot(reflection_matrix(axis), Cn(axis, n))

def get_symmetry_axes(mats, tol=1e-3):
    """
    Extracts unique rotation axes and their orders.
    Robustly handles float32 noise and det scaling.
    """
    unique_axes = [] # List of [order, axis]
    
    for mat in mats:
        # 1. Use sign of determinant to handle Improper Rotations without scaling magnitude
        det = np.linalg.det(mat)
        sign = 1.0 if det > 0 else -1.0
        R_part = mat * sign
        
        # 2. Eigen decomposition
        w, v = np.linalg.eig(R_part)
        
        # 3. Find eigenvalue ~ 1.0
        idx = np.argmin(np.abs(w - 1.0))
        if np.abs(w[idx] - 1.0) > tol: continue
            
        # 4. Extract Real Axis Robustly
        # Handle case where eig returns complex vector for real eigenvalue
        vec = v[:, idx]
        if np.all(np.abs(vec.imag) < tol):
            axis = np.real(vec)
        else:
            # Try to rotate phase to make it real (unlikely for valid R_part, but safe)
            # or just take real part if imag is small noise
            axis = np.real(vec)
            
        norm = np.linalg.norm(axis)
        if norm < tol: continue
        axis = axis / norm
        
        # 5. Enforce unique hemisphere
        for k in range(3):
            if abs(axis[k]) > tol:
                if axis[k] < 0: axis = -axis
                break

        # 6. Calculate Order
        tr = np.trace(R_part)
        cos_t = np.clip((tr - 1.0) / 2.0, -1.0, 1.0)
        theta = np.arccos(cos_t)
        
        if abs(theta) < tol: continue # Identity
            
        n = int(np.round(2 * np.pi / theta))
        if n < 2: continue

        # 7. Check Uniqueness & Upgrade Order
        found = False
        for i in range(len(unique_axes)):
            old_n, old_ax = unique_axes[i]
            if 1.0 - abs(np.dot(axis, old_ax)) < tol:
                found = True
                if n > old_n:
                    unique_axes[i] = [n, axis]
                break
        
        if not found:
            unique_axes.append([n, axis])
            
    unique_axes.sort(key=lambda x: x[0], reverse=True)
    return [tuple(x) for x in unique_axes]

def align_basis(u1, v1, u2, v2):
    """Constructs rotation matrix R aligning (u1, v1) to (u2, v2)."""
    # Frame 1
    x1 = u1 / np.linalg.norm(u1)
    y1 = v1 - np.dot(v1, x1) * x1
    if np.linalg.norm(y1) < 1e-4: return None
    y1 = y1 / np.linalg.norm(y1)
    z1 = np.cross(x1, y1)
    F1 = np.column_stack((x1, y1, z1))
    
    # Frame 2
    x2 = u2 / np.linalg.norm(u2)
    y2 = v2 - np.dot(v2, x2) * x2
    if np.linalg.norm(y2) < 1e-4: return None
    y2 = y2 / np.linalg.norm(y2)
    z2 = np.cross(x2, y2)
    F2 = np.column_stack((x2, y2, z2))
    
    return F2 @ F1.T

def check_transform(R, mats1, mats2, tol=1e-3):
    """Verifies if R @ mats1 @ R.T is a subset of mats2."""
    for m1 in mats1:
        m1_trans = R @ m1 @ R.T
        found = False
        # Optimization: Filter by trace first (rotation angle invariant)
        tr1 = np.trace(m1_trans)
        for m2 in mats2:
            if abs(tr1 - np.trace(m2)) > tol: continue
            if np.allclose(m1_trans, m2, atol=tol):
                found = True
                break
        if not found: return False
    return True

def transform_matrix(sym_mats_old, sym_mats_new):
    """
    Finds transformation R mapping sym_mats_old to sym_mats_new.
    Uses Geometric Hashing and Relaxed Order Checks.
    """
    mats1 = np.array(sym_mats_old)[:, :3, :3]
    mats2 = np.array(sym_mats_new)[:, :3, :3]

    if len(mats1) <= 1: return np.eye(3)

    axes1 = get_symmetry_axes(mats1)
    axes2 = get_symmetry_axes(mats2)

    if not axes1 or not axes2: return np.eye(3)

    # Strategy 1: Two-Vector Alignment
    n1_a, u1 = axes1[0]

    v1 = None
    for n, ax in axes1[1:]:
        if abs(np.dot(u1, ax)) < 1.0 - 1e-3:
            v1 = ax
            break

    if v1 is not None:
        target_angle = np.arccos(np.clip(np.dot(u1, v1), -1.0, 1.0))

        # Search all compatible pairs in target group
        for n2_a, u2 in axes2:
            # RELAXED CHECK: Iterate all u2 to handle order mismatches (e.g. C3 vs C6)

            for n2_b, v2 in axes2:
                # Check angle match
                ang = np.arccos(np.clip(np.dot(u2, v2), -1.0, 1.0))
                if abs(ang - target_angle) > 1e-3: continue

                R = align_basis(u1, v1, u2, v2)
                if R is not None and check_transform(R, mats1, mats2):
                    return R

    # Strategy 2: Single Axis Alignment
    for n2_a, u2 in axes2:
        # RELAXED CHECK: Try aligning primary axes even if orders mismatch

        dot = np.dot(u1, u2)
        if dot > 1.0 - 1e-4:
            R = np.eye(3)
        elif dot < -1.0 + 1e-4:
            perp = np.array([1.0, 0.0, 0.0])
            if abs(np.dot(perp, u1)) > 0.9: perp = np.array([0.0, 1.0, 0.0])
            rot_ax = np.cross(u1, perp)
            rot_ax /= np.linalg.norm(rot_ax)
            R = rotation_matrix(rot_ax, np.pi)
        else:
            rot_ax = np.cross(u1, u2)
            rot_ax /= np.linalg.norm(rot_ax)
            angle = np.arccos(np.clip(dot, -1.0, 1.0))
            R = rotation_matrix(rot_ax, angle)

        if check_transform(R, mats1, mats2):
            return R

    return np.eye(3)


def reduce(n, i):
    """
    Divide n and i by their greatest common divisor g.

    :type n: int
    :type i: int
    :return: Tuple of n/g and i/g
    :rtype: (int, int)
    """
    g = gcd(n, i)
    return n // g, i // g
    # floor divide to get an int, there should never be a remainder since we are dividing by the gcd


def gcd(A, B):
    """
    A quick implementation of the Euclid algorithm for finding the greatest common divisor between A and B.
    
    :type A: int
    :type B: int
    :return: Greatest common divisor between A and B
    :rtype: int
    """
    a = max(A, B)
    b = min(A, B)
    if a == 0:
        return b
    elif b == 0:
        return a
    else:
        r = a % b
        return gcd(b, r)

class PointGroup():
    """
    Class for defining point group.
    """

    def __init__(self, s, family, n, subfamily):
        self.str = s
        self.family = family
        self.n = n
        if self.n == 0:
            self.is_linear = True
        else:
            self.is_linear = False
        self.subfamily = subfamily
        self.dumb_pg()

    @classmethod
    def from_string(cls, s):
        regex = r"([A-Z]+)(\d+)?([a-z]+)?"
        m = re.match(regex, s)
        family, n, subfamily = m.groups()
        if n is not None:
            n = int(n)
        if subfamily is not None:
            subfamily = str(subfamily)
        family = str(family)
        return cls(s, family, n, subfamily)

    def __str__(self):
        nstr = self.n
        sfstr = self.subfamily
        if self.n is None:
            nstr = ""
        elif self.n == 0:
            nstr = "inf"
        else:
            nstr = str(self.n)
        if self.subfamily is None:
            sfstr = ""
        return self.family + nstr + sfstr

    def __repr__(self) -> str:
        return self.__str__()

    def dumb_pg(self):
        # Check if a dumb point group has been made (e.g. D1h, D0v, C2i)
        argstr = (f"You have generated a dumb point group:" +
                  " {self.str}. Family {self.family}, n {self.n}," +
                  " subfamily {self.subfamily}. We aren't sure " +
                  "how you managed to do this but we aren't paid " +
                  "enough to proceed with any calculations. " +
                  "If you have any questions, " +
                  "feel free to email the CFOUR listserv.")
        if self.n is None:
            if self.family == "C":
                allowed = ["s", "i"]
                if self.subfamily in allowed:
                    return 0
            elif self.family == "T":
                allowed = [None, "h", "d"]
                if self.subfamily in allowed:
                    return 0
            elif self.family == "O" or self.family == "I":
                allowed = [None, "h"]
                if self.subfamily in allowed:
                    return 0
        elif self.n == 0:
            if self.family == "D" and self.subfamily == "h":
                return 0
            elif self.family == "C" and self.subfamily == "v":
                return 0
        elif self.n == 1:
            if self.family == "C" and self.subfamily is None:
                return 0
        elif self.n >= 2:
            if self.family == "C":
                allowed = [None, "v", "h"]
                if self.subfamily in allowed:
                    return 0
            elif self.family == "D":
                allowed = [None, "d", "h"]
                if self.subfamily in allowed:
                    return 0
            elif self.family == "S":
                if self.subfamily is None and self.n % 2 == 0:
                    return 0
        raise Exception(argstr)


class CharacterTable():

    def __init__(self, pg, irreps, classes, class_orders, chars,
                 irrep_dims) -> None:
        self.name = pg
        self.irreps = irreps
        self.classes = classes
        self.class_orders = class_orders
        self.characters = chars
        self.irrep_dims = irrep_dims

    def __repr__(self) -> str:
        return f"Character Table for {self.name}\nIrreps: \
                {self.irreps}\nClasses: {self.classes}\nCharacters:\n{self.characters}\n"

    def __eq__(self, other):
        if len(self.irreps) == len(other.irreps) and \
                len(self.classes) == len(other.classes) and \
                np.shape(self.characters)==np.shape(other.characters):
            return (self.irreps == other.irreps).all() and \
                    (self.classes == other.classes).all() and \
                    np.isclose(self.characters,other.characters,atol=1e-10).all()
        else:
            return False


def pg_to_chartab(PG):
    pg = PointGroup.from_string(PG)
    irreps = []
    if pg.family == "C":
        if pg.subfamily == "s":
            irreps = ["A'", "A''"]
            classes = ["E", "sigma_h"]
            chars = np.array([[1.0, 1.0], [1.0, -1.0]])
        elif pg.subfamily == "i":
            irreps = ["Ag", "Au"]
            classes = ["E", "i"]
            chars = np.array([[1.0, 1.0], [1.0, -1.0]])
        elif pg.subfamily == "v":
            irreps, classes, chars = Cnv_irr(pg.n)
        elif pg.subfamily == "h":
            irreps, classes, chars = Cnh_irr(pg.n)
        else:
            #irreps, classes, chars = Cn_irrmat(pg.n)
            irreps, classes, chars = Cn_irr_complex(pg.n)
    elif pg.family == "D":
        if pg.subfamily == "d":
            irreps, classes, chars = Dnd_irr(pg.n)
        elif pg.subfamily == "h":
            irreps, classes, chars = Dnh_irr(pg.n)
        else:
            irreps, classes, chars = Dn_irr(pg.n)
    elif pg.family == "S":
        irreps, classes, chars = Sn_irr_complex(pg.n)
    else:
        cp3 = np.cos(np.pi / 3)
        pr5 = 0.5 * (1.0 + np.sqrt(5.0))
        mr5 = 0.5 * (1.0 - np.sqrt(5.0))
        if pg.family == "T":
            if pg.subfamily == "h":
                irreps, classes, chars = (
                    ["Ag", "Au", "Eg", "Eu", "Tg", "Tu"], [
                        "E", "4C_3", "4C_3^2", "3C_2", "i", "S_6", "S_6^5",
                        "3sigma_h"
                    ],
                    np.array([[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
                              [1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0],
                              [2.0, cp3, cp3, 2.0, 2.0, cp3, cp3, 1.0],
                              [2.0, cp3, cp3, 2.0, -2.0, -cp3, -cp3, -1.0],
                              [3.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0, -1.0],
                              [3.0, 0.0, 0.0, -1.0, -1.0, 0.0, 0.0, 1.0]]))
            elif pg.subfamily == "d":
                irreps, classes, chars = (["A1", "A2", "E", "T1", "T2"], [
                    "E", "8C_3", "3C_2", "6S_4", "6sigma_d"
                ],
                                          np.array([[1.0, 1.0, 1.0, 1.0, 1.0],
                                                    [1.0, 1.0, 1.0, -1.0, -1.0],
                                                    [2.0, -1.0, 2.0, 0.0, 0.0],
                                                    [3.0, 0.0, -1.0, 1.0, -1.0],
                                                    [3.0, 0.0, -1.0, -1.0,
                                                     1.0]]))
            else:
                irreps, classes, chars = (["A", "E", "T"],
                                          ["E", "4C_3", "4C_3^2", "3C_2"],
                                          np.array([[1.0, 1.0, 1.0, 1.0],
                                                    [2.0, cp3, cp3, 2.0],
                                                    [3.0, 0.0, 0.0, -1.0]]))
        elif pg.family == "O":
            if pg.subfamily == "h":
                irreps, classes, chars = (
                    [
                        "A1g", "A2g", "Eg", "T1g", "T2g", "A1u", "A2u", "Eu",
                        "T1u", "T2u"
                    ], [
                        "E", "8C_3", "6C_2", "6C_4", "3C_2", "i", "6S_4",
                        "8S_6", "3sigma_h", "6sigma_d"
                    ],
                    np.array([
                        [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
                        [1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0, -1.0],
                        [2.0, -1.0, 0.0, 0.0, 2.0, 2.0, 0.0, -1.0, 2.0, 0.0],
                        [3.0, 0.0, -1.0, 1.0, -1.0, 3.0, 1.0, 0.0, -1.0, -1.0],
                        [3.0, 0.0, 1.0, -1.0, -1.0, 3.0, -1.0, 0.0, -1.0, 1.0],
                        [1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0],
                        [1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 1.0, -1.0, -1.0, 1.0],
                        [2.0, -1.0, 0.0, 0.0, 2.0, -2.0, 0.0, 1.0, -2.0, 0.0],
                        [3.0, 0.0, -1.0, 1.0, -1.0, -3.0, -1.0, 0.0, 1.0, 1.0],
                        [3.0, 0.0, 1.0, -1.0, -1.0, -3.0, 1.0, 0.0, 1.0, -1.0]
                    ]))
            else:
                irreps, classes, chars = (["A1", "A2", "E", "T1", "T2"],
                                          ["E", "6C_4", "3C_2", "8C_3", "6C_2"],
                                          np.array([[1.0, 1.0, 1.0, 1.0, 1.0],
                                                    [1.0, -1.0, 1.0, 1.0, -1.0],
                                                    [2.0, 0.0, 2.0, -1.0, 0.0],
                                                    [3.0, 1.0, -1.0, 0.0, -1.0],
                                                    [3.0, -1.0, -1.0, 0.0,
                                                     1.0]]))
        elif pg.family == "I":
            if pg.subfamily == "h":
                irreps, classes, chars = (
                    [
                        "Ag", "T1g", "T2g", "Gg", "Hg", "Au", "T1u", "T2u",
                        "Gu", "Hu"
                    ], [
                        "E", "12C_5", "12C_5^2", "20C_3", "15C_2", "i",
                        "12S_10", "12S_10^3", "20S_6", "15sigma"
                    ],
                    np.array([
                        [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
                        [3.0, pr5, mr5, 0.0, -1.0, 3.0, mr5, pr5, 0.0, -1.0],
                        [3.0, mr5, pr5, 0.0, -1.0, 3.0, pr5, mr5, 0.0, -1.0],
                        [4.0, -1.0, -1.0, 1.0, 0.0, 4.0, -1.0, -1.0, 1.0, 0.0],
                        [5.0, 0.0, 0.0, -1.0, 1.0, 5.0, 0.0, 0.0, -1.0, 1.0],
                        [1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0],
                        [3.0, pr5, mr5, 0.0, -1.0, -3.0, -mr5, -pr5, 0.0, 1.0],
                        [3.0, mr5, pr5, 0.0, -1.0, -3.0, -pr5, -mr5, 0.0, 1.0],
                        [4.0, -1.0, -1.0, 1.0, 0.0, -4.0, 1.0, 1.0, -1.0, 0.0],
                        [5.0, 0.0, 0.0, -1.0, 1.0, -5.0, 0.0, 0.0, 1.0, -1.0]
                    ]))
            else:
                irreps, classes, chars = (["A", "T1", "T2", "G", "H"], [
                    "E", "12C_5", "12C_5^2", "20C_3", "15C_2"
                ],
                                          np.array([[1.0, 1.0, 1.0, 1.0, 1.0],
                                                    [3.0, pr5, mr5, 0.0, -1.0],
                                                    [3.0, mr5, pr5, 0.0, -1.0],
                                                    [4.0, -1.0, -1.0, 1.0, 0.0],
                                                    [5.0, 0.0, 0.0, -1.0,
                                                     1.0]]))
        else:
            raise Exception(f"An invalid point group has been " +
                            "given or unexpected parsing of the " +
                            "point group string has occured: {pg.str}")
    class_orders = grab_class_orders(classes)
    irr_dims = {}
    for (irr_idx, irrep) in enumerate(irreps):
        if pg.n == 1:
            irr_dims[irrep] = int(chars[0])
        else:
            irr_dims[irrep] = int(np.real(chars[irr_idx, 0]))
    return CharacterTable(PG, irreps, classes, class_orders, chars, irr_dims)


def grab_class_orders(classes):
    ncls = len(classes)
    class_orders = np.zeros(ncls)
    for i in range(ncls):  # = 1:ncls
        class_orders[i] = grab_order(classes[i])
    return class_orders


def grab_order(class_str):
    regex = r"^(\d+)"
    m = re.match(regex, class_str)
    if m is not None:
        return int(m.groups()[0])
    else:
        return 1


def c_class_string(classes, pre, r, s):
    "Pushes class string to 'classes' for rotations, but ignores superscript if s is one"
    if s == 1:
        classes.append(pre + f"_{r}")
    else:
        classes.append(pre + f"_{r}^{s}")


def Cn_irrmat(n):
    names = ["A"]
    classes = ["E"]
    for c in range(1, n):
        r, s = reduce(n, c)
        c_class_string(classes, "C", r, s)
    chars = np.ones(n)
    if n % 2 == 0:
        names.append("B")
        bi = np.ones(n)
        #for i=1:n:
        for i in range(n):
            if (i + 1) % 2 == 0:
                bi[i] *= -1
        chars = np.vstack((chars, bi))
    if 2 < n < 5:
        # No label associated with E if n < 5
        names.append("E")
        theta = 2 * np.pi / n
        v = np.zeros(n)
        for j in range(n):
            v[j] += 2 * np.cos(j * theta)
        chars = np.vstack((chars, v))
    elif n >= 5:
        theta = 2 * np.pi / n
        l = round(((n - len(names)) / 2))
        for i in range(l):  #= 1:l:
            names.append(f"E{i+1}")
            v = np.zeros(n)
            for j in range(n):
                v[j] += 2 * np.cos((i + 1) * j * theta)
            chars = np.vstack((chars, v))
    return names, classes, chars


def Cn_irr_complex(n):
    names = ["A"]
    classes = ["E"]
    for c in range(1, n):
        r, s = reduce(n, c)
        c_class_string(classes, "C", r, s)
    chars = np.ones(n)
    if n % 2 == 0:
        names.append("B")
        bi = np.ones(n)
        #for i=1:n:
        for i in range(n):
            if (i + 1) % 2 == 0:
                bi[i] *= -1
        chars = np.vstack((chars, bi))
    if 2 < n < 5:
        # No label associated with E if n < 5
        names.append("E_1")
        names.append("E_2")
        theta = np.exp(2 * np.pi * 1j / n)
        v1 = np.zeros(n, dtype=np.complex128)
        v2 = np.zeros(n, dtype=np.complex128)
        for a in range(n):
            v1[a] += theta**a
            v2[a] += np.conj(theta**a)
        chars = np.vstack((chars, v1))
        chars = np.vstack((chars, v2))
    elif n >= 5:
        theta = 2 * np.pi / n
        l = round(((n - len(names)) / 2))
        for a in range(l):  #= 1:l:
            names.append(f"E{a+1}_1")
            names.append(f"E{a+1}_2")
            theta = np.exp(2 * np.pi * 1j * (a + 1) / n)
            v1 = np.zeros(n, dtype=np.complex128)
            v2 = np.zeros(n, dtype=np.complex128)
            for b in range(n):
                v1[b] += theta**b
                v2[b] += np.conj(theta**b)
            chars = np.vstack((chars, v1))
            chars = np.vstack((chars, v2))
    return names, classes, chars


def Cnv_irr(n):
    names, classes, chars = Cn_irrmat(n)
    classes = ["E"]
    if n % 2 == 0:
        for c in range(1, n >> 1):  # = 1:(n>>1)-1:
            r, s = reduce(n, c)
            c_class_string(classes, "2C", r, s)
        classes.append("C_2")
        if n >> 1 == 1:
            classes.append("sigma_v(xz)")
            classes.append("sigma_d(yz)")
        else:
            classes.append(f"{n>>1}sigma_v")
            classes.append(f"{n>>1}sigma_d")
    else:
        for c in range(1, (n >> 1) + 1):  # = 1:n>>1:
            r, s = reduce(n, c)
            c_class_string(classes, "2C", r, s)
        classes.append(f"{n}sigma_v")
    names[0] = "A1"
    names.insert(1, "A2")
    chars = np.vstack((chars[0, :], chars[0, :], chars[1:, :]))
    for i in range(1, n):  # = 2:n:
        if i >= n - i:
            break
        #deleteat!(classes, n-i+2)
        chars = chars[:, [j for j in range(np.shape(chars)[1]) if j != n - i
                         ]]  #chars[:,1:-1.!=n-i+2]
    if n % 2 == 0:
        nirr = round((n / 2) + 3)
        names[2] = "B1"
        names.insert(3, "B2")
        chars = np.vstack((chars[0:3, :], chars[2, :], chars[3:, :]))
        sigma_v = np.zeros(nirr)
        sigma_d = np.zeros(nirr)
        sigma_v[0:4] = np.array([1, -1, 1, -1])
        sigma_d[0:4] = np.array([1, -1, -1, 1])
        chars = np.hstack((chars, sigma_v[:, None], sigma_d[:, None]))
    else:
        nirr = round((n - 1) / 2 + 2)
        sigma_v = np.zeros(nirr)
        sigma_v[0:2] = np.array([1, -1])
        chars = np.hstack((chars, sigma_v[:, None]))
    return names, classes, chars


def Cnh_irr(n):
    #names, classes, cnchars = Cn_irrmat(n)
    names, classes, cnchars = Cn_irr_complex(n)
    if n % 2 == 0:
        classes.append("i")
        for i in range(1, n):  # = 1:n-1:
            if i == n >> 1:
                classes.append("sigma_h")
            else:
                r, s = reduce(n, (i + (n >> 1)) % n)
                if s % 2 == 0:
                    s += r
                c_class_string(classes, "S", r, s)
    else:
        classes.append("sigma_h")
        for i in range(1, n):  # = 1:n-1:
            r, s = reduce(n, i)
            if i % 2 == 0:
                c_class_string(classes, "S", r, s + n)
            else:
                c_class_string(classes, "S", r, s)
    if n % 2 == 0:
        newnames = []
        for i in range(len(names)):  # = 1:length(names):
            newnames.append(names[i] + "u")
            names[i] = names[i] + "g"
        names += newnames
        cncharsi = -1 * cnchars
        top = np.hstack((cnchars, cnchars))
        bot = np.hstack((cnchars, cncharsi))
        chars = np.vstack((top, bot))
    else:
        newnames = []
        for i in range(len(names)):  # = 1:length(names):
            newnames.append(names[i] + "''")
            names[i] = names[i] + "'"
        names += newnames
        cncharsi = -1 * cnchars
        top = np.hstack((cnchars, cnchars))
        bot = np.hstack((cnchars, cncharsi))
        chars = np.vstack((top, bot))
    return names, classes, chars


def Sn_irr(n):
    if n % 4 == 0:
        names, classes, chars = Cn_irrmat(n)
        for i in range(n):  # = 1:n:
            if (i + 1) % 2 == 0:
                classes[i] = "S" + classes[i][1:]
    elif n % 2 == 0:
        ni = round(n / 2)
        names, classes, cnchars = Cn_irrmat(ni)
        classes = ["E"]
        for i in range(1, n >> 1):  # = 1:n>>1-1:
            r, s = reduce(n >> 1, i)
            c_class_string(classes, "C", r, s)
        classes.append("i")
        for i in range(1, n >> 1):  # = 1:n>>1-1:
            r, s = reduce(n, ((i << 1) + (n >> 1)) % n)
            c_class_string(classes, "S", r, s)
        newnames = []
        for i in range(len(names)):  # = 1:length(names):
            newnames.append(names[i] + "u")
            names[i] = names[i] + "g"
        names += newnames
        cncharsi = -1 * cnchars
        top = np.hstack((cnchars, cnchars))
        bot = np.hstack((cnchars, cncharsi))
        chars = np.vstack((top, bot))
    else:
        raise Exception("Odd number n for S group")
    return names, classes, chars


def Sn_irr_complex(n):
    if n % 4 == 0:
        names, classes, chars = Cn_irr_complex(n)
        for i in range(n):  # = 1:n:
            if (i + 1) % 2 == 0:
                classes[i] = "S" + classes[i][1:]
    elif n % 2 == 0:
        ni = round(n / 2)
        names, classes, cnchars = Cn_irr_complex(ni)
        classes = ["E"]
        for i in range(1, n >> 1):  # = 1:n>>1-1:
            r, s = reduce(n >> 1, i)
            c_class_string(classes, "C", r, s)
        classes.append("i")
        for i in range(1, n >> 1):  # = 1:n>>1-1:
            r, s = reduce(n, ((i << 1) + (n >> 1)) % n)
            c_class_string(classes, "S", r, s)
        newnames = []
        for i in range(len(names)):  # = 1:length(names):
            newnames.append(names[i] + "u")
            names[i] = names[i] + "g"
        names += newnames
        cncharsi = -1 * cnchars
        top = np.hstack((cnchars, cnchars))
        bot = np.hstack((cnchars, cncharsi))
        chars = np.vstack((top, bot))
    else:
        raise Exception("Odd number n for S group")
    return names, classes, chars


def Dn_irr(n):
    if n == 2:
        names = ["A", "B1", "B2", "B3"]
        classes = ["E", "C_2(z)", "C_2(y)", "C_2(x)"]
        chars = np.array([[1, 1, 1, 1], [1, 1, -1, -1], [1, -1, 1, -1],
                          [1, -1, -1, 1]])
        return names, classes, chars
    names, garbage_classes, chars = Cn_irrmat(n)
    garbage_names, classes, garbage_chars = Cnv_irr(n)
    if n % 2 == 0:
        classes[-2] = classes[-2][0] + "C_2'"
        classes[-1] = classes[-1][0] + "C_2''"
    else:
        classes[-1] = classes[-1][0] + "C_2"
    names[0] = "A1"
    names.insert(1, "A2")
    chars = np.vstack((chars[0, :], chars[0, :], chars[1:, :]))
    for i in range(1, n):  # = 2:n:
        if i == n - i or i > n - i:
            break
        #deleteat!(classes, n-i+2)
        chars = chars[:, [j for j in range(np.shape(chars)[1]) if j != n - i
                         ]]  # chars[:,0:-1.!=n-i+2]
    if n % 2 == 0:
        nirr = round((n / 2) + 3)
        names[2] = "B1"
        names.insert(3, "B2")
        chars = np.vstack((chars[0:3, :], chars[2, :], chars[3:, :]))
        C2p = np.zeros(nirr)
        C2pp = np.zeros(nirr)
        C2p[0:4] = np.array([1, -1, 1, -1])
        C2pp[0:4] = np.array([1, -1, -1, 1])
        chars = np.hstack((chars, C2p[:, None], C2pp[:, None]))
    else:
        nirr = round((n - 1) / 2 + 2)
        C2p = np.zeros(nirr)
        C2p[0:2] = np.array([1, -1])
        chars = np.hstack((chars, C2p[:, None]))
    return names, classes, chars


def Dnh_irr(n):
    names, classes, dnchars = Dn_irr(n)
    if n % 2 == 0:
        classes.append("i")
        if n == 2:
            classes.append("sigma(xy)")
            classes.append("sigma(xz)")
            classes.append("sigma(yz)")
        else:
            for i in range(1, n >> 1):  # = 1:n>>1-1:
                a = i + (n >> 1)
                if a > (n >> 1):
                    a = n - a
                r, s = reduce(n, a)
                c_class_string(classes, "2S", r, s)
            classes.append("sigma_h")
            if n % 4 == 0:
                classes.append(f"{n>>1}sigma_v")
                classes.append(f"{n>>1}sigma_d")
            else:
                classes.append(f"{n>>1}sigma_d")
                classes.append(f"{n>>1}sigma_v")
        newnames = []
        for i in range(len(names)):  # = 1:length(names):
            newnames.append(names[i] + "u")
            names[i] = names[i] + "g"
        names += newnames
        dncharsi = -1 * dnchars
        top = np.hstack((dnchars, dnchars))
        bot = np.hstack((dnchars, dncharsi))
        chars = np.vstack((top, bot))
    else:
        classes.append("sigma_h")
        for i in range(1, (n >> 1) + 1):  # = 1:n>>1:
            if i % 2 == 0:
                r, s = reduce(n, n - i)
            else:
                r, s = reduce(n, i)
            c_class_string(classes, "2S", r, s)
        classes.append(f"{n}sigma_v")
        newnames = []
        for i in range(len(names)):  # = 1:length(names):
            newnames.append(names[i] + "''")
            names[i] = names[i] + "'"
        names += newnames
        dncharsi = -1 * dnchars
        top = np.hstack((dnchars, dnchars))
        bot = np.hstack((dnchars, dncharsi))
        chars = np.vstack((top, bot))
    return names, classes, chars


def Dnd_irr(n):
    if n % 2 == 0:
        n2 = 2 * n
        names, classes, chars = Sn_irr(n2)
        #classes = collect(1:2*n2)
        classes = classes[0:n + 1]
        for i in range(1, n):  # = 2:n:
            classes[i] = "2" + classes[i]
        classes.append(f"{n}C_2'")
        classes.append(f"{n}sigma_d")
        names[0] = "A1"
        names.insert(1, "A2")
        chars = np.vstack((chars[0, :], chars[0, :], chars[1:, :]))
        for i in range(1, n2):  # = 2:n2:
            if i >= n2 - i:
                break
            chars = chars[:,
                          [j for j in range(np.shape(chars)[1]) if j != n2 - i
                          ]]  # chars[:,0:-1.!=n2-i+2]
        nirr = n + 3
        names[2] = "B1"
        names.insert(3, "B2")
        chars = np.vstack((chars[0:3, :], chars[2, :], chars[3:, :]))
        C2p = np.zeros(nirr)
        sigma_d = np.zeros(nirr)
        C2p[0:4] = np.array([1, -1, 1, -1])
        sigma_d[0:4] = np.array([1, -1, -1, 1])
        chars = np.hstack((chars, C2p[:, None], sigma_d[:, None]))
    else:
        names, classes, dnchars = Dn_irr(n)
        classes.append("i")
        for i in range(1, (n >> 1) + 1):  # = 1:n>>1:
            r, s = reduce(2 * n, 2 * i + n)
            if s > n:
                s = 2 * n - s
            c_class_string(classes, "2S", r, s)
        classes.append(f"{n}sigma_d")
        newnames = []
        for i in range(len(names)):  # = 1:length(names):
            newnames.append(names[i] + "u")
            names[i] = names[i] + "g"
        names += newnames
        dncharsi = -1 * dnchars
        top = np.hstack((dnchars, dnchars))
        bot = np.hstack((dnchars, dncharsi))
        chars = np.vstack((top, bot))
    return names, classes, chars


@dataclass
class Symel():
    """
    Deprecated data structure for symmetry elements.
    Still hanging around because of cubic/icosahedral groups and tests.
    """
    symbol: str
    vector: np.array  # Not defined for E or i, axis vector for Cn and Sn, plane normal vector for sigma
    rrep: np.array

    def __str__(self) -> str:
        with np.printoptions(precision=5,
                             suppress=True,
                             formatter={"all": lambda x: f"{x:8.5f}"}):
            return f"\nSymbol: {self.symbol:>10s}: [{self.rrep[0,:]},{self.rrep[1,:]},{self.rrep[2,:]}]"

    def __repr__(self) -> str:
        return self.__str__()

    def __eq__(self, other):
        return self.symbol == other.symbol and np.isclose(
            self.rrep, other.rrep, atol=1e-10).all()


def pg_to_symels(PG):
    pg = PointGroup.from_string(PG)
    argerr = f"An invalid point group has been given or "+\
        "unexpected parsing of the point group string has occured: {pg.str}"
    symels = [Symel("E", None, np.asarray([[1, 0, 0], [0, 1, 0], [0, 0, 1]]))]
    z_axis = np.array([0, 0, 1])
    sigma_h = np.asarray([[1, 0, 0], [0, 1, 0], [0, 0, -1]])
    if pg.family == "C":
        if pg.subfamily == "h":
            symels.append(Symel("sigma_h", z_axis, sigma_h))
            if pg.n % 2 == 0:
                symels.append(Symel("i", None, inversion_matrix()))
            cns = generate_Cn(pg.n)
            sns = generate_Sn(pg.n)
            symels = symels + cns + sns
        elif pg.subfamily == "v":
            cns = generate_Cn(pg.n)
            if pg.n % 2 == 0:
                n = pg.n >> 1
                sigma_ds = generate_sigma_d(n)
            else:
                n = pg.n
                sigma_ds = []
            sigma_vs = generate_sigma_v(pg.n)
            symels = symels + cns + sigma_vs + sigma_ds
        elif pg.subfamily == "s":
            symels.append(Symel("sigma_h", z_axis, sigma_h))
        elif pg.subfamily == "i":
            symels.append(Symel("i", None, inversion_matrix()))
        elif pg.subfamily is None:
            cns = generate_Cn(pg.n)
            symels = symels + cns
        else:
            raise Exception(argerr)
    elif pg.family == "D":
        if pg.subfamily == "h":
            symels.append(Symel("sigma_h", z_axis, sigma_h))
            if pg.n % 2 == 0:
                symels.append(Symel("i", None, inversion_matrix()))
                n = pg.n >> 1
                sigma_ds = generate_sigma_d(n)
                c2ps = generate_C2p(pg.n)
                c2pps = generate_C2pp(pg.n)
                c2s = c2ps + c2pps
            else:
                n = pg.n
                sigma_ds = []
                c2s = generate_C2p(pg.n)
            cns = generate_Cn(pg.n)
            sns = generate_Sn(pg.n)
            sigma_vs = generate_sigma_v(pg.n)
            #c2s = generate_C2(pg.n)
            symels = symels + cns + c2s + sns + sigma_vs + sigma_ds
        elif pg.subfamily == "d":
            if pg.n % 2 == 0:
                c2ps = generate_C2p(pg.n)
                c2pps = generate_C2pp(pg.n)
                c2s = c2ps + c2pps
            else:
                c2s = generate_C2p(pg.n)
                symels.append(Symel("i", None, inversion_matrix()))
            cns = generate_Cn(pg.n)
            sns = generate_Sn(pg.n * 2, S2n=True)
            sigma_ds = generate_sigma_d(pg.n)
            symels = symels + cns + sns + c2s + sigma_ds
        elif pg.subfamily is None:
            cns = generate_Cn(pg.n)
            if pg.n % 2 == 0:
                c2ps = generate_C2p(pg.n)
                c2pps = generate_C2pp(pg.n)
                c2s = c2ps + c2pps
            else:
                c2s = generate_C2p(pg.n)
            symels = symels + cns + c2s
        else:
            raise Exception(argerr)
    elif pg.family == "S":
        if pg.subfamily is None and (pg.n % 2 == 0):
            n = pg.n >> 1
            if n % 2 != 0:
                symels.append(Symel("i", None, inversion_matrix()))
            cns = generate_Cn(n)
            sns = generate_Sn(pg.n, S2n=True)
            symels = symels + cns + sns
        else:
            raise Exception(argerr)
    else:
        if pg.family == "T":
            if pg.subfamily == "h":
                Ths = generate_Th()
                symels = symels + Ths
            elif pg.subfamily == "d":
                Tds = generate_Td()
                symels = symels + Tds
            else:
                Ts = generate_T()
                symels = symels + Ts
        elif pg.family == "O":
            if pg.subfamily == "h":
                Ohs = generate_Oh()
                symels = symels + Ohs
            else:
                Os = generate_O()
                symels = symels + Os
        elif pg.family == "I":
            if pg.subfamily == "h":
                Ihs = generate_Ih()
                symels = symels + Ihs
            else:
                Is = generate_I()
                symels = symels + Is
        else:
            raise Exception(argerr)
    return symels


def generate_Cn(n):
    symels = []
    axis = np.asarray([0, 0, 1])
    #axis = [0 0 1]'
    cn_r = Cn(axis, n)
    for i in range(1, n):
        a, b = reduce(n, i)
        symels.append(Symel(f"C_{a:d}^{b:d}", axis, matrix_power(cn_r,
                                                                 i)))  # Cns
    return symels


def generate_Sn(n, S2n=False):
    symels = []
    axis = np.asarray([0, 0, 1])
    sigma_h = reflection_matrix(axis)
    cn_r = Cn(axis, n)
    if S2n:  # Generating improper rotations for S2n PG
        for i in range(1, n):
            if i % 2 == 0:
                continue
            else:
                a, b = reduce(n, i)
                if a == 2:
                    continue
                else:
                    symels.append(
                        Symel(f"S_{a}^{b}", axis,
                              np.dot(matrix_power(cn_r, i), sigma_h)))
        return symels
    for i in range(1, n):
        a, b = reduce(n, i)
        if b % 2 == 0:
            b += a
        if a == 2:
            continue
        else:
            symels.append(
                Symel(f"S_{a}^{b}", axis, np.dot(matrix_power(cn_r, i),
                                                 sigma_h)))  # Sns
    return symels


def generate_sigma_v(n):
    if n % 2 == 0:
        nsigma_vs = n >> 1
    else:
        nsigma_vs = n
    symels = []
    x_axis = np.asarray([1, 0, 0])  # Orient C2 and sigma_v along x-axis
    z_axis = np.asarray([0, 0, 1])
    rot_mat = Cn(z_axis, n)
    for i in range(nsigma_vs):
        axis = np.cross(np.dot(matrix_power(rot_mat, i), x_axis), z_axis)
        symels.append(Symel(f"sigma_v({i+1})", axis, reflection_matrix(axis)))
    return symels


def generate_sigma_d(n):
    symels = []
    x_axis = np.asarray([1, 0, 0])  # Orient C2 and sigma_v along x-axis
    z_axis = np.asarray([0, 0, 1])
    rot_mat = Cn(z_axis, 2 * n)
    base_axis = np.dot(
        Cn(z_axis, 4 * n),
        x_axis)  # Rotate x-axis by Cn/2 to produce an axis for sigma_d's
    for i in range(n):
        axis = np.cross(np.dot(matrix_power(rot_mat, i), base_axis), z_axis)
        symels.append(Symel(f"sigma_d({i+1})", axis, reflection_matrix(axis)))
    return symels


def generate_C2p(n):
    if n % 2 == 0:
        nn = n >> 1
    else:
        nn = n
    symels = []
    x_axis = np.asarray([1, 0, 0])  # Orient C2 and sigma_v along x-axis
    rot_mat = Cn([0, 0, 1], n)
    for i in range(nn):
        axis = np.dot(matrix_power(rot_mat, i), x_axis)
        symels.append(Symel(f"C_2'({i+1})", axis, Cn(axis, 2)))
    return symels


def generate_C2pp(n):
    nn = n >> 1
    symels = []
    x_axis = np.asarray([1, 0, 0])
    rot_mat = Cn([0, 0, 1], n)
    base_axis = np.dot(Cn([0, 0, 1], 2 * n), x_axis)
    for i in range(nn):
        axis = np.dot(matrix_power(rot_mat, i), base_axis)
        symels.append(Symel(f"C_2''({i+1})", axis, Cn(axis, 2)))
    return symels


def generate_T():
    """
    Generate symmetry elements for the T point group.
    Assume a tetrahedron contained in a cube, then we can easily generate
    the vectors for the rotation elements.

    :rtype: List[molsym.Symel]
    """
    # Generate C3's
    symels = []
    C3_1v = normalize(np.array([1.0, 1.0, 1.0]))
    C3_2v = normalize(np.array([-1.0, 1.0, -1.0]))
    C3_3v = normalize(np.array([-1.0, -1.0, 1.0]))
    C3_4v = normalize(np.array([1.0, -1.0, -1.0]))
    C3list = [C3_1v, C3_2v, C3_3v, C3_4v]
    namelist = ("alpha", "beta", "gamma", "delta")
    for i in range(4):
        C3 = Cn(C3list[i], 3)
        C3s = matrix_power(C3, 2)
        symels.append(Symel(f"C_3({namelist[i]})", C3list[i], C3))
        symels.append(Symel(f"C_3^2({namelist[i]})", C3list[i], C3s))
    # Generate C2's
    C2_x = np.array([1.0, 0.0, 0.0])
    C2_y = np.array([0.0, 1.0, 0.0])
    C2_z = np.array([0.0, 0.0, 1.0])
    C2list = [C2_x, C2_y, C2_z]
    namelist = ["x", "y", "z"]
    for i in range(3):
        C2 = Cn(C2list[i], 2)
        symels.append(Symel(f"C_2({namelist[i]})", C2list[i], C2))
    return symels


def generate_Td():
    """
    Generate symmetry elements for the Td point group.
    Assume a tetrahedron contained in a cube, then we can easily generate
    the vectors for the rotation elements.

    :rtype: List[molsym.Symel]
    """
    symels = generate_T()
    # σd's
    sigma_d_1v = normalize(np.array([1.0, 1.0, 0.0]))
    sigma_d_2v = normalize(np.array([1.0, -1.0, 0.0]))
    sigma_d_3v = normalize(np.array([1.0, 0.0, 1.0]))
    sigma_d_4v = normalize(np.array([1.0, 0.0, -1.0]))
    sigma_d_5v = normalize(np.array([0.0, 1.0, 1.0]))
    sigma_d_6v = normalize(np.array([0.0, 1.0, -1.0]))
    sigmas = [
        sigma_d_1v, sigma_d_2v, sigma_d_3v, sigma_d_4v, sigma_d_5v, sigma_d_6v
    ]
    namelist = ["xyp", "xym", "xzp", "xzm", "yzp", "yzm"]
    for i in range(6):
        sigma_d = reflection_matrix(sigmas[i])
        symels.append(Symel(f"sigma_d({namelist[i]})", sigmas[i], sigma_d))
    # S4's
    S4_1v = np.array([1.0, 0.0, 0.0])
    S4_2v = np.array([0.0, 1.0, 0.0])
    S4_3v = np.array([0.0, 0.0, 1.0])
    S4vlist = [S4_1v, S4_2v, S4_3v]
    namelist = ["x", "y", "z"]
    for i in range(3):
        S4 = Sn(S4vlist[i], 4)
        S43 = matrix_power(S4, 3)
        symels.append(Symel(f"S_4({namelist[i]})", S4vlist[i], S4))
        symels.append(Symel(f"S_4^3({namelist[i]})", S4vlist[i], S43))
    return symels


def generate_Th():
    """
    Generate symmetry elements for the Th point group.
    Assume a tetrahedron contained in a cube, then we can easily generate
    the vectors for the rotation elements.

    :rtype: List[molsym.Symel]
    """
    symels = generate_T()
    # i
    symels.append(Symel("i", None, inversion_matrix()))
    # S6
    S6_1v = normalize(np.array([1.0, 1.0, 1.0]))
    S6_2v = normalize(np.array([-1.0, 1.0, -1.0]))
    S6_3v = normalize(np.array([-1.0, -1.0, 1.0]))
    S6_4v = normalize(np.array([1.0, -1.0, -1.0]))
    S6list = [S6_1v, S6_2v, S6_3v, S6_4v]
    namelist = ["alpha", "beta", "gamma", "delta"]
    for i in range(4):
        S6 = Sn(S6list[i], 6)
        S65 = matrix_power(S6, 5)
        symels.append(Symel(f"S_6({namelist[i]})", S6list[i], S6))
        symels.append(Symel(f"S_6^5({namelist[i]})", S6list[i], S65))
    # 3sigma_h
    sigma_h_xv = np.array([1.0, 0.0, 0.0])
    sigma_h_yv = np.array([0.0, 1.0, 0.0])
    sigma_h_zv = np.array([0.0, 0.0, 1.0])
    sigma_list = [sigma_h_xv, sigma_h_yv, sigma_h_zv]
    namelist = ["x", "y", "z"]
    for i in range(3):
        sigma_h = reflection_matrix(sigma_list[i])
        symels.append(Symel(f"sigma_h({namelist[i]})", sigma_list[i], sigma_h))
    return symels


def generate_O():
    """
    Generate symmetry elements for the O point group.
    Assume operations on a cube.

    :rtype: List[molsym.Symel]
    """
    symels = []
    # C4
    C4_xv = np.array([1.0, 0.0, 0.0])
    C4_yv = np.array([0.0, 1.0, 0.0])
    C4_zv = np.array([0.0, 0.0, 1.0])
    C4list = [C4_xv, C4_yv, C4_zv]
    namelist = ["x", "y", "z"]
    for i in range(3):
        C4 = Cn(C4list[i], 4)
        C42 = matrix_power(C4, 2)
        C43 = matrix_power(C4, 3)
        symels.append(Symel(f"C_4({namelist[i]})", C4list[i], C4))
        symels.append(Symel(f"C_2({namelist[i]})", C4list[i], C42))
        symels.append(Symel(f"C_4^3({namelist[i]})", C4list[i], C43))
    # C3
    C3_1v = normalize(np.array([1.0, 1.0, 1.0]))
    C3_2v = normalize(np.array([1.0, -1.0, 1.0]))
    C3_3v = normalize(np.array([1.0, 1.0, -1.0]))
    C3_4v = normalize(np.array([1.0, -1.0, -1.0]))
    C3list = [C3_1v, C3_2v, C3_3v, C3_4v]
    namelist = ["alpha", "beta", "gamma", "delta"]
    for i in range(4):
        C3 = Cn(C3list[i], 3)
        C32 = matrix_power(C3, 2)
        symels.append(Symel(f"C_3({namelist[i]})", C3list[i], C3))
        symels.append(Symel(f"C_3^2({namelist[i]})", C3list[i], C32))
    # C2
    C2_1v = normalize(np.array([1.0, 0.0, 1.0]))
    C2_2v = normalize(np.array([1.0, 0.0, -1.0]))
    C2_3v = normalize(np.array([1.0, 1.0, 0.0]))
    C2_4v = normalize(np.array([1.0, -1.0, 0.0]))
    C2_5v = normalize(np.array([0.0, 1.0, 1.0]))
    C2_6v = normalize(np.array([0.0, -1.0, 1.0]))

    C2list = [C2_1v, C2_2v, C2_3v, C2_4v, C2_5v, C2_6v]
    namelist = ["xzp", "xzm", "xyp", "xym", "yzp", "yzm"]
    for i in range(6):
        C2 = Cn(C2list[i], 2)
        symels.append(Symel(f"C_2({namelist[i]})", C2list[i], C2))
    return symels


def generate_Oh():
    """
    Generate symmetry elements for the Oh point group.
    Assume operations on a cube.

    :rtype: List[molsym.Symel]
    """
    symels = generate_O()
    symels.append(Symel("i", None, inversion_matrix()))
    # S4 and σh
    S4_xv = np.array([1.0, 0.0, 0.0])
    S4_yv = np.array([0.0, 1.0, 0.0])
    S4_zv = np.array([0.0, 0.0, 1.0])
    S4list = [S4_xv, S4_yv, S4_zv]
    namelist = ["x", "y", "z"]
    for i in range(3):
        S4 = Sn(S4list[i], 4)
        sigma_h = reflection_matrix(S4list[i])
        S43 = matrix_power(S4, 3)
        symels.append(Symel(f"S_4({namelist[i]})", S4list[i], S4))
        symels.append(Symel(f"sigma_h({namelist[i]})", S4list[i], sigma_h))
        symels.append(Symel(f"S_4^3({namelist[i]})", S4list[i], S43))
    # S6
    S6_1v = normalize(np.array([1.0, 1.0, 1.0]))
    S6_2v = normalize(np.array([1.0, -1.0, 1.0]))
    S6_3v = normalize(np.array([1.0, 1.0, -1.0]))
    S6_4v = normalize(np.array([1.0, -1.0, -1.0]))
    S6list = [S6_1v, S6_2v, S6_3v, S6_4v]
    namelist = ["alpha", "beta", "gamma", "delta"]
    for i in range(4):
        S6 = Sn(S6list[i], 6)
        S65 = matrix_power(S6, 5)
        symels.append(Symel(f"S_6({namelist[i]})", S6list[i], S6))
        symels.append(Symel(f"S_6^5({namelist[i]})", S6list[i], S65))
    # C2
    sigma_d_1v = normalize(np.array([1.0, 0.0, 1.0]))
    sigma_d_2v = normalize(np.array([1.0, 0.0, -1.0]))
    sigma_d_3v = normalize(np.array([1.0, 1.0, 0.0]))
    sigma_d_4v = normalize(np.array([1.0, -1.0, 0.0]))
    sigma_d_5v = normalize(np.array([0.0, 1.0, 1.0]))
    sigma_d_6v = normalize(np.array([0.0, -1.0, 1.0]))

    sigma_dlist = [
        sigma_d_1v, sigma_d_2v, sigma_d_3v, sigma_d_4v, sigma_d_5v, sigma_d_6v
    ]
    namelist = ["xzp", "xzm", "xyp", "xym", "yzp", "yzm"]
    for i in range(6):
        sigma_d = reflection_matrix(sigma_dlist[i])
        symels.append(Symel(f"sigma_d({namelist[i]})", sigma_dlist[i], sigma_d))
    return symels


def generate_I():
    """
    Generate symmetry elements for the I point group.

    :rtype: List[molsym.Symel]
    """
    symels = []
    faces, vertices, edgecenters = generate_I_vectors()
    # C5 (face vectors)
    for i in range(6):
        C5 = Cn(faces[i], 5)
        C52 = matrix_power(C5, 2)
        C53 = matrix_power(C5, 3)
        C54 = matrix_power(C5, 4)
        symels.append(Symel(f"C_5({i})", faces[i], C5))
        symels.append(Symel(f"C_5^2({i})", faces[i], C52))
        symels.append(Symel(f"C_5^3({i})", faces[i], C53))
        symels.append(Symel(f"C_5^4({i})", faces[i], C54))

    # C3 (vertex vectors)
    for i in range(10):
        C3 = Cn(vertices[i], 3)
        C32 = matrix_power(C3, 2)
        symels.append(Symel(f"C_3({i})", vertices[i], C3))
        symels.append(Symel(f"C_3^2({i})", vertices[i], C32))

    # C2 (edge vectors)
    for i in range(15):
        C2 = Cn(edgecenters[i], 2)
        symels.append(Symel(f"C_2({i})", edgecenters[i], C2))

    return symels


def generate_Ih():
    """
    Generate symmetry elements for the Ih point group.

    :rtype: List[molsym.Symel]
    """
    symels = generate_I()
    faces, vertices, edgecenters = generate_I_vectors()
    symels.append(Symel("i", None, inversion_matrix()))
    # S10 (face vectors)
    for i in range(6):
        S10 = Sn(faces[i], 10)
        S103 = matrix_power(S10, 3)
        S107 = matrix_power(S10, 7)
        S109 = matrix_power(S10, 9)
        symels.append(Symel(f"S_10({i})", faces[i], S10))
        symels.append(Symel(f"S_10^3({i})", faces[i], S103))
        symels.append(Symel(f"S_10^7({i})", faces[i], S107))
        symels.append(Symel(f"S_10^9({i})", faces[i], S109))

    # S6 (vertex vectors)
    for i in range(10):
        S6 = Sn(vertices[i], 6)
        S65 = matrix_power(S6, 5)
        symels.append(Symel(f"S_6({i})", vertices[i], S6))
        symels.append(Symel(f"S_6^5({i})", vertices[i], S65))

    # σ (edge vectors)
    for i in range(15):
        sigma_i = reflection_matrix(edgecenters[i])
        symels.append(Symel(f"sigma({i})", edgecenters[i], sigma_i))

    return symels


def generate_I_vectors():
    """
    Vectors defining the faces, vertices, and edges of a regular dodecahedron.

    :rtype: (List[NumPy array of shape (3,)], List[NumPy array of shape (3,)], List[NumPy array of shape (3,)])
    """
    #phi = (1+(np.sqrt(5.0)))/2.0
    #phi_i = 1.0/phi
    #faces_i = np.array([[1.0, phi, 0.0],[1.0, -phi, 0.0],[-1.0, phi, 0.0],[-1.0, -phi, 0.0],
    #         [0.0, 1.0, phi],[0.0, 1.0, -phi],[0.0, -1.0, phi],[0.0, -1.0, -phi],
    #         [phi, 0.0, 1.0],[-phi, 0.0, 1.0],[phi, 0.0, -1.0],[-phi, 0.0, -1.0]])
    #vertices_i = np.array([[1.0, 1.0, 1.0],[1.0, 1.0, -1.0],[1.0, -1.0, 1.0],[-1.0, 1.0, 1.0],
    #            [1.0, -1.0, -1.0],[-1.0, 1.0, -1.0],[-1.0, -1.0, 1.0],[-1.0, -1.0, -1.0],
    #            [0.0, phi, phi_i],[0.0, phi, -phi_i],[0.0, -phi, phi_i],[0.0, -phi, -phi_i],
    #            [phi_i, 0.0, phi],[-phi_i, 0.0, phi],[phi_i, 0.0, -phi],[-phi_i, 0.0, -phi],
    #            [phi, phi_i, 0.0],[phi, -phi_i, 0.0],[-phi, phi_i, 0.0],[-phi, -phi_i, 0.0]])
    ## Reorienting vectors such that one face is on the z-axis with "pentagon" pointing at the POSITIVE?  y-axis
    #theta = -np.arccos(phi/np.sqrt(1+(phi**2)))
    #rmat = rotation_matrix(np.array([1, 0, 0]), theta)
    #lf = np.shape(faces_i)[0]
    #l = np.shape(vertices_i)[0]
    #faces = [np.dot(rmat,faces_i[i,:]) for i in range(lf)]
    #vertices = [np.dot(rmat,vertices_i[i,:]) for i in range(l)]
    #faces = np.dot(rmat, faces_i.T).T
    #vertices = np.dot(rmat, vertices_i.T).T

    #l = np.shape(vertices)[0]
    #edglen = 2*phi_i
    #edgecenters = []
    #for i in range(l):
    #    for j in range(i+1,l):#= i+1:l
    #        if np.isclose(distance(vertices[i], vertices[j]), edglen):
    #            v = normalize(vertices[i]+vertices[j])
    #            same = False
    #            for k in edgecenters:
    #                if np.isclose(abs(np.dot(k,v)), 1.0):
    #                    same = True
    #                    break
    #            if not same:
    #                edgecenters.append(v)
    #for (idx,face) in enumerate(faces):
    #    faces[idx] = normalize(face)
    #face_vectors = []
    #same = False
    #for i in faces:
    #    v = normalize(i)
    #    for j in face_vectors:
    #        same = False
    #        if np.isclose(abs(np.dot(v,j)), 1.0):
    #            same = True
    #            break
    #    if not same:
    #        face_vectors.append(v)
    #
    #vertex_vectors = []
    #same = False
    #for i in vertices:
    #    v = normalize(i)
    #    for j in vertex_vectors:
    #        same = False
    #        if np.isclose(abs(np.dot(v,j)), 1.0):
    #            same = True
    #            break
    #    if not same:
    #        vertex_vectors.append(v)
    face_vectors = [
        np.array([0.0, 0.0, 1.0]),
        np.array([0.0, 2 / np.sqrt(5), 1 / np.sqrt(5)]),
        np.array([
            -np.sqrt((5 + np.sqrt(5)) / 10), (5 - np.sqrt(5)) / 10,
            1 / np.sqrt(5)
        ]),
        np.array([
            np.sqrt((5 + np.sqrt(5)) / 10), (5 - np.sqrt(5)) / 10,
            1 / np.sqrt(5)
        ]),
        np.array([
            -np.sqrt((5 - np.sqrt(5)) / 10), -(5 + np.sqrt(5)) / 10,
            1 / np.sqrt(5)
        ]),
        np.array([
            np.sqrt((5 - np.sqrt(5)) / 10), -(5 + np.sqrt(5)) / 10,
            1 / np.sqrt(5)
        ])
    ]
    vertex_vectors = [
        np.array([
            0.0, -np.sqrt((2 * (5 + np.sqrt(5))) / 15),
            np.sqrt((5 - 2 * np.sqrt(5)) / 15)
        ]),
        np.array([
            0.0, -np.sqrt((2 * (5 - np.sqrt(5))) / 15),
            np.sqrt((5 + 2 * np.sqrt(5)) / 15)
        ]),
        np.array([
            -1 / np.sqrt(3),
            np.sqrt((5 + 2 * np.sqrt(5)) / 15),
            np.sqrt((5 - 2 * np.sqrt(5)) / 15)
        ]),
        np.array([
            -1 / np.sqrt(3), -np.sqrt((5 - 2 * np.sqrt(5)) / 15),
            np.sqrt((5 + 2 * np.sqrt(5)) / 15)
        ]),
        np.array([
            1 / np.sqrt(3),
            np.sqrt((5 + 2 * np.sqrt(5)) / 15),
            np.sqrt((5 - 2 * np.sqrt(5)) / 15)
        ]),
        np.array([
            1 / np.sqrt(3), -np.sqrt((5 - 2 * np.sqrt(5)) / 15),
            np.sqrt((5 + 2 * np.sqrt(5)) / 15)
        ]),
        np.array([
            -np.sqrt((3 * np.sqrt(5) + 5) / (6 * np.sqrt(5))), -np.sqrt(
                (5 - np.sqrt(5)) / 30),
            np.sqrt((5 - 2 * np.sqrt(5)) / 15)
        ]),
        np.array([
            np.sqrt((3 * np.sqrt(5) - 5) / (6 * np.sqrt(5))),
            np.sqrt((5 + np.sqrt(5)) / 30),
            np.sqrt((5 + 2 * np.sqrt(5)) / 15)
        ]),
        np.array([
            -np.sqrt((3 * np.sqrt(5) - 5) / (6 * np.sqrt(5))),
            np.sqrt((5 + np.sqrt(5)) / 30),
            np.sqrt((5 + 2 * np.sqrt(5)) / 15)
        ]),
        np.array([
            np.sqrt((3 * np.sqrt(5) + 5) / (6 * np.sqrt(5))), -np.sqrt(
                (5 - np.sqrt(5)) / 30),
            np.sqrt((5 - 2 * np.sqrt(5)) / 15)
        ]),
    ]
    edgecenters = [
        np.array([
            0.0,
            np.sqrt((5 + np.sqrt(5)) / 10), -np.sqrt((5 - np.sqrt(5)) / 10)
        ]),
        np.array([
            0.0,
            np.sqrt((5 - np.sqrt(5)) / 10),
            np.sqrt((5 + np.sqrt(5)) / 10)
        ]),
        np.array([
            0.5, -np.sqrt(1 + (2 / np.sqrt(5))) / 2, -np.sqrt(
                (5 - np.sqrt(5)) / 10)
        ]),
        np.array([
            0.5,
            np.sqrt(1 - (2 / np.sqrt(5))) / 2,
            np.sqrt((5 + np.sqrt(5)) / 10)
        ]),
        np.array([
            0.5, -np.sqrt(1 - (2 / np.sqrt(5))) / 2, -np.sqrt(
                (5 + np.sqrt(5)) / 10)
        ]),
        np.array([
            0.5,
            np.sqrt(1 + (2 / np.sqrt(5))) / 2,
            np.sqrt((5 - np.sqrt(5)) / 10)
        ]),
        np.array([1.0, 0.0, 0.0]),
        np.array([(np.sqrt(5) - 1) / 4, -np.sqrt((5 + np.sqrt(5)) / 10) / 2,
                  np.sqrt((5 + np.sqrt(5)) / 10)]),
        np.array([(np.sqrt(5) - 1) / 4,
                  np.sqrt((5 + np.sqrt(5)) / 10) / 2, -np.sqrt(
                      (5 + np.sqrt(5)) / 10)]),
        np.array([(np.sqrt(5) - 1) / 4, -np.sqrt((5 + np.sqrt(5)) / 2) / 2,
                  0.0]),
        np.array([(np.sqrt(5) - 1) / 4,
                  np.sqrt((5 + np.sqrt(5)) / 2) / 2, 0.0]),
        np.array([(np.sqrt(5) + 1) / 4, -np.sqrt((5 - np.sqrt(5)) / 10) / 2,
                  np.sqrt((5 - np.sqrt(5)) / 10)]),
        np.array([(np.sqrt(5) + 1) / 4,
                  np.sqrt((5 - np.sqrt(5)) / 10) / 2, -np.sqrt(
                      (5 - np.sqrt(5)) / 10)]),
        np.array([(np.sqrt(5) + 1) / 4, -np.sqrt((5 - np.sqrt(5)) / 2) / 2,
                  0.0]),
        np.array([(np.sqrt(5) + 1) / 4,
                  np.sqrt((5 - np.sqrt(5)) / 2) / 2, 0.0]),
    ]

    return (face_vectors, vertex_vectors, edgecenters)


def point_group_classes(point_group, tol=1e-4):
    """
    Given a list of symmetries that form a point
    group. Identify the classes
    point_group : list of symmetries
    """
    # point_group = (nsym,3,3)
    pg_classes = []

    nsym = len(point_group)
    pg_tree = KDTree(point_group.reshape(nsym, -1))

    for isym in range(nsym):
        sym_conjugates = point_group.transpose(
            0, 2, 1) @ point_group[isym][None, :, :] @ point_group
        dd, ii = pg_tree.query(sym_conjugates.reshape(nsym, -1), k=1)
        if len(dd[dd > tol]) != 0:
            exit("Not a point group")
        s1_tmp = set(ii)
        if s1_tmp not in pg_classes:
            pg_classes.append(s1_tmp)

    return [list(item) for item in set(tuple(row) for row in pg_classes)]

def find_axis_angle(Rmat):
    """
    Given an orthogonal matrix 
    find its axes and angle.
    incase of +-np.eye(3), axes is [0,0,0]
    and nfold = 0
    """
    det = np.linalg.det(Rmat)
    A = Rmat * det
    if np.abs(A - np.eye(3)).sum() < 1e-4:
        return 0, np.array([0, 0, 0])
    ###
    w, v = np.linalg.eig(A)
    one_idx = np.argmin(np.abs(w - 1))
    assert np.abs(w[one_idx] - 1) < 1e-4
    axis = v[:, one_idx]
    axis = axis / np.linalg.norm(axis)
    # find angle
    cos_angle = np.trace(A) * 0.5 - 0.5
    Kmat = np.zeros(9, dtype=complex)
    Kmat[1] = -axis[2]
    Kmat[2] = axis[1]
    Kmat[3] = axis[2]
    Kmat[5] = -axis[0]
    Kmat[6] = -axis[1]
    Kmat[7] = axis[0]
    sin_angle = -0.5 * np.trace(Kmat.reshape(3, 3) @ A).real
    if abs(sin_angle) < 1e-3:
        sin_angle = abs(sin_angle)
    angle = np.arctan2(sin_angle, cos_angle)
    return angle, axis


def fix_axis_angle_gauge(axis, nfold):
    """
    Fix the gauge of the provided axis. and changes sign
    of nfold according to that
    The gauge is fixed such that the real part of [z/y/x]
    is real.
    """
    axis_tmp = np.real(axis)
    for i in [2, 1, 0]:
        if abs(axis_tmp[i]) > 1e-4:
            if axis_tmp[i] < 0:
                return -np.array(axis), -nfold
            else:
                return axis, nfold
    return axis, nfold


def find_symm_axis(sym_mats):
    ## find symmmats and nfold degeneracy
    ### in case not found, nfold is returned as None
    """
    Given a list of symmetry matrices that form a point,
    return their determinents, axis and nfold
    for mirrors : nfold == 2
    """
    group_order = len(sym_mats)
    dets = np.linalg.det(sym_mats)

    axes = []
    nfold = []

    for isym in sym_mats:
        ## check if this is S:
        w = np.linalg.eigvals(isym)
        impro_rot = False
        if np.abs(w - 1).min() > 1e-3:
            impro_rot = True
        t, ax = find_axis_angle(isym)
        if impro_rot and abs(t) > 1e-2:
            if t > 0:
                t -= np.pi
            else:
                t += np.pi
        if (abs(t) > 1e-2):
            t = int(np.rint(2 * np.pi / t))
            ax, t = fix_axis_angle_gauge(ax, t)
        else:
            t = 0
        if abs(t) <= 2:
            t = abs(t)
        axes.append(ax)
        nfold.append(t)

    axes = np.array(axes)
    nfold = np.array(nfold, dtype=int)
    return [dets, nfold, axes]

def get_point_grp(symm_mats):
    """
    Identifies the point group label from symmetry matrices.
    Correctly handles D2d vs D2 ambiguity by checking improper axes.
    """
    dets, nfold, axes = find_symm_axis(symm_mats)
    order = len(symm_mats)
    if order == 1: return "C1"
    
    # Identify Proper Rotations
    proper_indices = np.where(np.logical_and(dets > 0, nfold > 0))[0]
    if len(proper_indices) == 0:
        if np.any(dets < 0):
             return "Cs" if np.abs(nfold).max() > 0 else "Ci"
        return "C1"

    nfold_proper = nfold[proper_indices]
    axes_proper = axes[proper_indices]
    n = np.max(nfold_proper)
    
    # --- Principal Axis Selection (Critical for D2d) ---
    candidate_mask = (nfold_proper == n)
    candidate_axes = axes_proper[candidate_mask]
    paxis = candidate_axes[0]
    
    # If multiple max-order axes exist (e.g., D2, D2h, D2d), favor one collinear with S_2n
    if n == 2 and len(candidate_axes) > 1:
        improper_indices = np.where(dets < 0)[0]
        if len(improper_indices) > 0:
            nfold_imp = np.abs(nfold[improper_indices])
            axes_imp = axes[improper_indices]
            # S4 has nfold=4 in rotation part
            if np.any(nfold_imp > n): 
                target_axes = axes_imp[nfold_imp > n]
                for cand in candidate_axes:
                    if np.any(np.abs(np.abs(target_axes @ cand) - 1.0) < 1e-3):
                        paxis = cand
                        break

    # --- Classification ---
    nfold_abs = np.abs(nfold)
    inversion_present = len(nfold_abs[nfold_abs == 0]) > 1
    
    # Cubic/Icosahedral checks
    c3_axes = axes[np.logical_and(dets > 0, nfold_abs == 3)]
    # Filter unique C3 axes
    c3_unique = []
    for ax in c3_axes:
        if not any(np.abs(np.dot(u, ax)) > 1-1e-3 for u in c3_unique): c3_unique.append(ax)
        
    if len(c3_unique) > 1:
        has_c5 = np.any(np.logical_and(dets > 0, nfold_abs == 5))
        has_c4 = np.any(np.logical_and(dets > 0, nfold_abs == 4))
        if has_c5: return "Ih" if inversion_present else "I"
        if has_c4: return "Oh" if inversion_present else "O"
        if inversion_present: return "Th"
        return "Td" if np.any(np.logical_and(dets < 0, nfold_abs == 4)) else "T"

    # Axial Checks
    perp_c2_axes = []
    c2_candidates = axes[np.logical_and(dets > 0, nfold_abs == 2)]
    for ax in c2_candidates:
        if np.abs(np.dot(ax, paxis)) < 1e-3:
            if not any(np.abs(np.dot(u, ax)) > 1-1e-3 for u in perp_c2_axes):
                perp_c2_axes.append(ax)
    nC_perp = len(perp_c2_axes)
    
    mirror_axes = axes[np.logical_and(dets < 0, nfold_abs == 2)]
    has_sigma_h = any(np.abs(np.abs(np.dot(ax, paxis)) - 1.0) < 1e-3 for ax in mirror_axes)
    nsig_v = sum(1 for ax in mirror_axes if np.abs(np.dot(ax, paxis)) < 1e-3)

    if n == nC_perp: # D groups
        if has_sigma_h: return "D%dh" % n
        if nsig_v > 0: return "D%dd" % n
        return "D%d" % n
    else: # C/S groups
        if has_sigma_h: return "C%dh" % n
        if nsig_v > 0: return "C%dv" % n
        if np.any(np.logical_and(dets < 0, nfold_abs == 2 * n)): return "S%d" % (2 * n)
        return "C%d" % n



def generate_symel_to_class_map(symels, ctab):
    """
    Deprecated. Generate a map of the symels to their corresponding classes.
    """
    pg = PointGroup.from_string(ctab.name)
    if pg.n is not None:
        ns = pg.n >> 1  # pg.n floor divided by 2
    ncls = len(ctab.classes)
    nsymel = len(symels)
    class_map = np.zeros(nsymel, dtype=np.int32)
    class_map[0] = 0  # E is always first
    if pg.family == "C":
        if pg.subfamily == "s" or pg.subfamily == "i":
            class_map[1] = 1
        elif pg.subfamily == "h":
            if pg.n % 2 == 0:
                class_map[3:pg.n + 2] = [i for i in range(1, pg.n)
                                        ]  # [2:pg.n] # C_n
                class_map[2] = pg.n  # i
                class_map[1] = pg.n + ns  # σh
                for i in range(pg.n + 2, 2 * pg.n):  # = pg.n+3:2*pg.n # S_n
                    if i > 3 * ns:
                        class_map[i] = i - ns
                    else:
                        class_map[i] = i + ns - 1
            else:
                for i in range(1, pg.n):  # = 2:pg.n+1 # C_n
                    class_map[i + 1] = i
                class_map[1] = pg.n  # σh
                for i in range(pg.n + 1, 2 * pg.n):  # = pg.n+2:2*pg.n # S_n
                    class_map[i] = i
        elif pg.subfamily == "v":
            # The last class is σv (and then σd if n is even), and the last symels are also these!
            cn_class_map(class_map, pg.n, 0, 0)
            if pg.n % 2 == 0:
                class_map[-pg.n:-ns] = ncls - 2
                class_map[-ns:] = ncls - 1
            else:
                class_map[-pg.n:] = ncls - 1
        else:
            class_map[1:] = [i for i in range(1, nsymel)]  # 2:nsymel
    elif pg.family == "S":
        if pg.n % 4 == 0:
            for i in range(1, pg.n):  # = 2:pg.n
                if i <= ns - 1:
                    class_map[i] = 2 * i
                else:
                    class_map[i] = 2 * (i - ns) + 1
        else:
            class_map[1] = ns  # i
            class_map[2:ns + 1] = [i for i in range(1, ns)]  # 2:ns # C_n
            for i in range(ns + 1, pg.n):  # = ns+2:pg.n # S_n
                if i > ns + (pg.n >> 2):
                    class_map[i] = i - (pg.n >> 2)
                else:
                    class_map[i] = i + (pg.n >> 2)
    elif pg.family == "D":
        if pg.subfamily == "h":
            if pg.n % 2 == 0:
                class_map[1] = ncls - 3  # σh
                class_map[2] = (ncls >> 1)  # i
                cn_class_map(class_map, pg.n, 2, 0)  # Cn
                class_map[pg.n + 2:3 * ns + 2] = ns + 1  # C2'
                class_map[3 * ns + 2:2 * pg.n + 2] = ns + 2  # C2''
                for i in range(2 * pg.n + 2,
                               3 * pg.n):  # = 2*pg.n+3:3*pg.n+1 # Sn
                    if i > 3 * pg.n - ns:
                        class_map[i] = i - 2 * pg.n + 3
                    else:
                        class_map[i] = 3 * pg.n + 4 - i
                # The result of C2'×i changes depending if pg.n ≡ 0 (mod 4)
                # but also D2h doesn't need to be flipped because I treated it special
                if pg.n % 4 == 0 or pg.n == 2:
                    class_map[-pg.n:-ns] = ncls - 2  # σv
                    class_map[-ns:] = ncls - 1  # σd
                else:
                    class_map[-pg.n:-ns] = ncls - 1  # σv
                    class_map[-ns:] = ncls - 2  # σd
            else:
                class_map[1] = (ncls >> 1)
                cn_class_map(class_map, pg.n, 1, 0)
                class_map[pg.n + 1:2 * pg.n + 1] = ns + 1
                cn_class_map(class_map, pg.n, 2 * pg.n, ns + 2)
                class_map[-1 - pg.n + 1:] = ncls - 1
        elif pg.subfamily == "d":
            if pg.n % 2 == 0:
                cn_class_map(class_map, pg.n, 0, 0)  # Cn
                class_map[1:pg.n] = 2 * class_map[
                    1:pg.n]  # 2*class_map[2:pg.n].-1 # Reposition Cn
                cn_class_map(class_map, pg.n + 1, pg.n - 1, 0)  # Sn
                class_map[pg.n:2 * pg.n] = 2 * class_map[pg.n:2 * pg.n] - 1
                # 2*(class_map[pg.n+1:2*pg.n].-1) # Reposition Sn
                class_map[-2 * pg.n:-pg.n] = ncls - 2  # C2'
                class_map[-pg.n:] = ncls - 1  # σd
            else:
                class_map[1] = ncls >> 1  # i
                cn_class_map(class_map, pg.n, 1, 0)  # Cn
                for i in range(pg.n + 1, 2 * pg.n):  # = pg.n+2:2*pg.n # Sn
                    if i > pg.n + ns:
                        class_map[i] = i + 2 - pg.n
                    else:
                        class_map[i] = 2 * pg.n + 2 - i
                class_map[-2 * pg.n:-pg.n] = ns + 1
                class_map[-pg.n:] = ncls - 1  # σd
        else:
            cn_class_map(class_map, pg.n, 0, 0)  # Cn
            if pg.n % 2 == 0:
                class_map[-pg.n:-ns] = ncls - 2  # Cn'
                class_map[-ns:] = ncls - 1  # Cn''
            else:
                class_map[-pg.n:] = ncls - 1  # Cn
    else:
        if pg.family == "T":
            if pg.subfamily == "h":
                class_map = np.array([
                    0, 1, 2, 1, 2, 1, 2, 1, 2, 3, 3, 3, 4, 5, 6, 5, 6, 5, 6, 5,
                    6, 7, 7, 7
                ])
            elif pg.subfamily == "d":
                class_map = np.array([
                    0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 4, 4, 4, 4, 4, 4, 3, 3,
                    3, 3, 3, 3
                ])
            else:
                class_map = np.array([0, 1, 2, 1, 2, 1, 2, 1, 2, 3, 3, 3])
        elif pg.family == "O":
            if pg.subfamily == "h":
                class_map = np.array([
                    0, 3, 4, 3, 3, 4, 3, 3, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2,
                    2, 2, 2, 2, 5, 6, 8, 6, 6, 8, 6, 6, 8, 6, 7, 7, 7, 7, 7, 7,
                    7, 7, 9, 9, 9, 9, 9, 9
                ])
            else:
                class_map = np.array([
                    0, 1, 2, 1, 1, 2, 1, 1, 2, 1, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4,
                    4, 4, 4, 4
                ])
        elif pg.family == "I":
            if pg.subfamily == "h":
                class_map = np.array([
                    0, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2,
                    1, 1, 2, 2, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                    3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    5, 6, 7, 7, 6, 6, 7, 7, 6, 6, 7, 7, 6, 6, 7, 7, 6, 6, 7, 7,
                    6, 6, 7, 7, 6, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                    8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9
                ])
            else:
                class_map = np.array([
                    0, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2,
                    1, 1, 2, 2, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                    3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
                ])
        else:
            raise Exception(f"An invalid point group has been " +
                            "given or unexpected parsing of the " +
                            "point group string has occured: {pg.str}")
    return class_map


def cn_class_map(class_map, n, idx_offset, cls_offset):
    """
    Deprecated
    """
    for i in range(1, n):  # = 2:n
        if i > (n >> 1):
            class_map[i + idx_offset] = n - i + cls_offset
        else:
            class_map[i + idx_offset] = i + cls_offset
    return class_map






### ============ Unit test ======================
### Generated using chatgpt ####################
### =============================================
if __name__ == '__main__':
    import unittest
    import sys

    ntrials = 50
    class TestPointGroupRobustness(unittest.TestCase):

        def setUp(self):
            # List of groups to test (Crystallographic + Icosahedral)
            self.test_groups = [
                # Triclinic
                "C1", "Ci",
                # Monoclinic
                "C2", "Cs", "C2h",
                # Orthorhombic
                "D2", "C2v", "D2h",
                # Tetragonal
                "C4", "S4", "C4h", "D4", "C4v", "D2d", "D4h",
                # Trigonal
                "C3", "S6", "D3", "C3v", "D3d",
                # Hexagonal
                "C6", "C3h", "C6h", "D6", "C6v", "D3h", "D6h",
                # Cubic
                "T", "Th", "Td", "O", "Oh"
            ]

            self.num_trials = ntrials


        def test_get_pg_info_robustness(self):
            """
            Tests get_pg_info by rotating the molecule randomly in 3D space.
            """
            print("\n")
            print("="*50)
            print("="*20 + "Test 1" + "="*20)
            print("="*50)
            print(f"\nRunning Robustness Test on {len(self.test_groups)} Point Groups")
            print(f"Configurations per group: {self.num_trials}")
            print("=" * 70)

            for label in self.test_groups:
                with self.subTest(group=label):
                    # ----------------------------------------------------------
                    # 1. Generate Ground Truth (Standard Orientation)
                    # ----------------------------------------------------------
                    try:
                        symels = pg_to_symels(label)
                        # Use float64 for precision
                        std_mats = np.array([s.rrep for s in symels], dtype=np.float64)
                    except Exception as e:
                        self.fail(f"[{label}] Setup Error: Could not generate standard symels. {e}")

                    # ----------------------------------------------------------
                    # 2. Run Random Trials
                    # ----------------------------------------------------------
                    for i in range(self.num_trials):
                        R, _ = np.linalg.qr(np.random.randn(3, 3))
                        transformed_mats = np.einsum('ij,njk,kl->nil', R, std_mats, R.T,optimize=True)
                        np.random.shuffle(transformed_mats)
                        mats_clean = []
                        for m in transformed_mats:
                            U, S, Vt = np.linalg.svd(m)
                            m_clean = U @ Vt
                            if np.linalg.det(m_clean) * np.linalg.det(m) < 0:
                                m_clean = -m_clean
                            mats_clean.append(m_clean)
                        transformed_mats = np.array(mats_clean)

                        try:
                            res_label, classes, class_dict, char_tab, irreps = get_pg_info(transformed_mats)
                        except AssertionError as e:
                            # This catches "assert len(pg_sym_mats) == order" and "assert distance < 1e-4"
                            self.fail(f"[{label}] Trial {i+1} FAILED: Library Assertion Error. "
                                  f"The algorithm failed to identify or align the rotated group. "
                                  f"\nMatrix Det: {np.linalg.det(R):.2f}")
                        except Exception as e:
                            self.fail(f"[{label}] Trial {i+1} CRASHED: {e}")

                        self.assertEqual(res_label, label,
                                     f"[{label}] Trial {i+1}: Misclassified as '{res_label}'")

                        if char_tab.ndim == 1:
                            ct_rows, ct_cols = 1, char_tab.shape[0]
                        else:
                            ct_rows, ct_cols = char_tab.shape

                        self.assertEqual(len(irreps), ct_rows,
                                     f"[{label}] CharTable rows {ct_rows} != Irreps {len(irreps)}")
                        self.assertEqual(len(classes), ct_cols,
                                     f"[{label}] CharTable cols {ct_cols} != Classes {len(classes)}")

                        assigned_indices = []
                        for idx_list in class_dict.values():
                            assigned_indices.extend(idx_list)

                        self.assertEqual(len(assigned_indices), len(std_mats),
                                     f"[{label}] Class dictionary missing matrices")
                        self.assertEqual(len(set(assigned_indices)), len(std_mats),
                                     f"[{label}] Class dictionary contains duplicates")

            print("\nSUCCESS: All groups passed Test 1")

    class SymmetryBaseXX:
        """
        Port of the symmetry definitions from Quantum Espresso's symm_base.f90.
        """
        def __init__(self):
            # Constants
            sin3 = 0.866025403784438597
            cos3 = 0.5
            msin3 = -0.866025403784438597
            mcos3 = -0.5

            # Raw data from symm_base.f90
            # This is a flat list of the 32 matrices (3x3).
            # Fortran DATA statements fill arrays in column-major order.
            raw_s0 = [
                # 1-9 (Identity, 180s, 90s)
                1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,   # 1: Identity
                -1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0, # 2: 180 z
                -1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0, # 3: 180 y
                1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0, # 4: 180 x
                0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0,  # 5: 180 [1,1,0]
                0.0, -1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, -1.0,# 6: 180 [1,-1,0]
                0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,  # 7: 90 z
                0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0,  # 8: -90 z

                # 10-17 (More 180s, 90s)
                0.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, 0.0,  # 9: 180 [1,0,1]
                0.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, 0.0,# 10: 180 [-1,0,1]
                0.0, 0.0, 1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0,  # 11: 90 y
                0.0, 0.0, -1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0,  # 12: -90 y
                0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,   # 13: 180 [0,1,1]
                0.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 1.0, 0.0, # 14: 180 [0,1,-1]
                1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 0.0,  # 15: 90 x
                1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0,  # 16: -90 x

                # 18-25 (Cubic Diagonals / 120 deg)
                0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0,   # 17: 120 [-1,-1,-1]
                0.0, -1.0, 0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, # 18: 120 [-1,1,1]
                0.0, -1.0, 0.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0, # 19: 120 [1,1,-1]
                0.0, 1.0, 0.0, 0.0, 0.0, -1.0, -1.0, 0.0, 0.0, # 20: 120 [1,-1,1]
                0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,   # 21: 120 [1,1,1]
                0.0, 0.0, -1.0, 1.0, 0.0, 0.0, 0.0, -1.0, 0.0, # 22: 120 [-1,1,-1]
                0.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 1.0, 0.0, # 23: 120 [1,-1,-1]
                0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, -1.0, 0.0, # 24: 120 [-1,-1,1]

                # 26-29 (Hexagonal 60/120 z)
                cos3, msin3, 0.0, sin3, cos3, 0.0, 0.0, 0.0, 1.0,   # 25: 60 z
                cos3, sin3, 0.0, msin3, cos3, 0.0, 0.0, 0.0, 1.0,   # 26: -60 z
                mcos3, msin3, 0.0, sin3, mcos3, 0.0, 0.0, 0.0, 1.0, # 27: 120 z
                mcos3, sin3, 0.0, msin3, mcos3, 0.0, 0.0, 0.0, 1.0, # 28: -120 z

                # 30-32 (Hexagonal 180s)
                cos3, msin3, 0.0, msin3, mcos3, 0.0, 0.0, 0.0, -1.0,# 29: 180 [1,-1,0]
                cos3, sin3, 0.0, sin3, mcos3, 0.0, 0.0, 0.0, -1.0,  # 30: 180 [2,1,0]
                mcos3, msin3, 0.0, msin3, cos3, 0.0, 0.0, 0.0, -1.0,# 31: 180 [0,1,0]
                mcos3, sin3, 0.0, sin3, cos3, 0.0, 0.0, 0.0, -1.0   # 32: 180 [1,1,0]
            ]

            # Reshape: (32, 3, 3).
            # Note: Fortran fills columns first. A reshape(32,3,3) in numpy fills rows.
            # We must reshape to (3,3,32) with Fortran order, then transpose axes.
            s0_fortran = np.array(raw_s0).reshape((3, 3, 32), order='F')

            # Move the 32 axis to the front: (32, 3, 3)
            self.s0 = np.moveaxis(s0_fortran, -1, 0)

            # Create Inversions (Indices 32-63)
            self.s_full = np.concatenate((self.s0, -self.s0), axis=0)

        def get_group_by_generators(self, gen_indices):
            """
            Constructs a full group by closing the set of generators.
            gen_indices: List of indices into s_full (0-63).
            """
            generators = [self.s_full[i] for i in gen_indices]

            # Closure: Keep multiplying until set size is constant
            group = list(generators)

            # Add Identity explicitly if not present (Index 0)
            id_mat = self.s_full[0]
            if not any(np.allclose(g, id_mat, atol=1e-5) for g in group):
                group.append(id_mat)

            while True:
                current_len = len(group)
                new_mats = []
                for m1 in group:
                    for m2 in group:
                        prod = m1 @ m2
                        # Check if prod exists in group
                        is_new = True
                        for existing in group:
                            if np.allclose(prod, existing, atol=1e-5):
                                is_new = False
                                break
                        if is_new:
                            # Also check against new_mats to avoid duplicates in this batch
                            for nm in new_mats:
                                if np.allclose(prod, nm, atol=1e-5):
                                    is_new = False
                                    break
                        if is_new:
                            new_mats.append(prod)

                if not new_mats:
                    break
                group.extend(new_mats)

            return np.array(group)

    class TestAll32Groups(unittest.TestCase):
        @classmethod
        def setUpClass(cls):
            cls.qe = SymmetryBaseXX()

            # =========================================================================
            # Map 32 Point Groups 
            # 0-31: Proper Rotations | 32-63: Inversions (Inv * s0)
            # =========================================================================
            I = 32 # Offset for inversion

            # Indices derived from symm_base.f90 comments
            # 0: E
            # 1: C2z, 2: C2y, 3: C2x
            # 6: C4z (90 z), 7: C4z^-1
            # 24: C6z (60 z) -> s0 index 25 in Fortran (1-based) -> 24 Python
            # 26: C3z (120 z)
            # 16: C3 diagonal (120 [-1,-1,-1])

            cls.generators = {
                # Triclinic
                "C1": [0],
                "Ci": [0 + I],

                # Monoclinic (Standard: z-axis)
                "C2": [1],       # C2z
                "Cs": [1 + I],   # sigma_h (Inv * C2z)
                "C2h": [1, 0+I], # C2z, i

                # Orthorhombic
                "D2": [1, 3],      # C2z, C2x
                "C2v": [1, 3 + I], # C2z, sigma_v(x) (Inv * C2x)
                "D2h": [1, 3, 0+I],# C2z, C2x, i

                # Tetragonal
                "C4": [6],         # C4z
                "S4": [6 + I],     # S4z (Inv * C4z)
                "C4h": [6, 0+I],   # C4z, i
                "D4": [6, 3],      # C4z, C2x
                "C4v": [6, 3 + I], # C4z, sigma_v(x)
                "D2d": [6 + I, 3], # S4z, C2x (Classic D2d generators)
                "D4h": [6, 3, 0+I],

                # Trigonal (Using Hexagonal basis indices 24-31)
                "C3": [26],           # C3z (120 deg)
                "S6": [26, 0 + I],    # C3z, i (S6 = C3 x i)
                "D3": [26, 28],       # C3z, C2' (180 perp)
                "C3v": [26, 28 + I],  # C3z, sigma_v
                "D3d": [26, 28, 0+I], # C3z, C2', i

                # Hexagonal
                "C6": [24],           # C6z
                "C3h": [26, 1 + I],   # C3z, sigma_h (Inv * C2z is sigma_h perp to z)
                "C6h": [24, 0 + I],   # C6z, i
                "D6": [24, 28],       # C6z, C2'
                "C6v": [24, 28 + I],  # C6z, sigma_v
                "D3h": [26, 28, 1+I], # C3z, C2', sigma_h
                "D6h": [24, 28, 0+I],

                # Cubic
                "T": [3, 16],      # C2x, C3_diag
                "Th": [3, 16, 0+I],# T + i
                "O": [6, 16],      # C4z, C3_diag
                "Td": [16, 6 + I], # C3_diag, S4z
                "Oh": [6, 16, 0+I] # O + i
            }

        def test_randomized_32_groups(self):
            """
            Generates all 32 groups and applies random rotations,
            and verifies detection.
            """
            n_trials = ntrials # Number of random rotations per group

            print("\n")
            print("="*50)
            print("="*20 + "Test 2" + "="*20)
            print("="*50)
            print(f"Testing All 32 Point Groups with {n_trials} Random Rotations")
            print(f"{'='*70}")

            for label, gen_idx in self.generators.items():
                with self.subTest(group=label):
                    # 1. Generate Standard Group
                    mats = self.qe.get_group_by_generators(gen_idx)

                    # Verify Order (Sanity Check)
                    # print(f"Testing {label:<4} (Order {len(mats)})")

                    for i in range(n_trials):
                        # 2. Random Rotation
                        # qr decomposition of random matrix gives a uniform random rotation
                        R, _ = np.linalg.qr(np.random.randn(3, 3))

                        # Apply rotation: R * M * R.T
                        mats_rot = np.einsum('ij,njk,kl->nil', R, mats, R.T)

                        # Clean numerical noise: Force matrices to be perfectly orthogonal
                        mats_clean = []
                        for m in mats_rot:
                            U, S, Vt = np.linalg.svd(m)
                            # Reconstruct as a perfect rotation matrix (U * Vt)
                            m_clean = U @ Vt
                            # Fix determinant sign (ensure it doesn't flip improper/proper)
                            # (Ideally SVD preserves it, but checking det is safer)
                            if np.linalg.det(m_clean) * np.linalg.det(m) < 0:
                                m_clean = -m_clean
                            mats_clean.append(m_clean)
                        mats_rot = np.array(mats_clean)

                    # 3. Detect
                        try:
                            detected, _, _, _, _ = get_pg_info(mats_rot)
                        except Exception as e:
                            self.fail(f"[{label}] Crash on trial {i}: {e}")

                        # 4. Assert
                        # Note: S6 is often detected as C3i (equivalent)
                        #       C3h often C6... no, C3h is distinct.
                        #       Cs can be Cs or C1h (synonyms).

                        expected = label
                        # Handle common synonyms in output
                        if label == "S6" and detected == "C3i": detected = "S6"
                        if label == "C3i" and detected == "S6": detected = "C3i"
                        if label == "Cs" and detected == "C1h": detected = "Cs"

                        self.assertEqual(detected, expected,
                                         f"Trial {i}: Expected {expected}, got {detected}. "
                                         f"\n(Group Order: {len(mats)})")

            print("\nSUCCESS: All 32 groups passed Test 2")


    unittest.main(verbosity=2)

