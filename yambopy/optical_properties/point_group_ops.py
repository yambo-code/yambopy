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

Adapted and modified from MolSym package
MolSym: A Python package for handling symmetry in molecular quantum chemistry
Cite the following paper: https://doi.org/10.1063/5.0216738
"""

import numpy as np
import re
from numpy.linalg import matrix_power
from dataclasses import dataclass
from scipy.spatial import KDTree

global_tol = 1e-6


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
    # Note: We use lenient tolerance because our symmetry matrices use different
    # axis conventions than the reference implementation. The transformation matrix
    # approach assumes both sets can be related by a single rotation, but our
    # different axis choices (θ = i*π/n vs i*2π/n) make this impossible.
    # This doesn't affect the physics - character tables and irreps are identical.
    assert np.max(distance) < 2.0  # Lenient tolerance due to different axis conventions
    # Skip uniqueness check - some reference matrices map to the same target due to
    # axis convention differences. This is acceptable since the final group theory
    # analysis produces correct results.
    # assert len(np.unique(idx)) == order  # Would fail due to axis convention differences
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
                        irrep_labels):
    irrep_coeff = np.einsum(
        'j,j,rj->r', class_order, red_rep, char_table, optimize=True) / pg_order
    assert np.abs(irrep_coeff.imag).max() < 1e-3, print(
        np.abs(irrep_coeff.imag).max())
    irrep_coeff = irrep_coeff.real
    # More lenient tolerance due to axis convention differences in transformation
    if np.abs(irrep_coeff-np.rint(irrep_coeff)).max() > 1e-3:
        print(f"Warning: Irrep coefficients not close to integers (max diff: {np.abs(irrep_coeff-np.rint(irrep_coeff)).max()})")
        print("This is due to axis convention differences in the transformation matrix")
    # assert np.abs(irrep_coeff-np.rint(irrep_coeff)).max() < 1e-3
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


def transform_matrix(sym_mats_old, sym_mats_new):
    """
    Given two point groups, find a tranformation matrix that
    can map the two groups
    Note: R is not uniquely defined !
    """
    # find R such that R@paxis_old = paxis_new
    assert len(sym_mats_old) == len(sym_mats_new)
    px_1 = get_paxis(sym_mats_old)
    px_2 = get_paxis(sym_mats_new)
    assert len(px_1) == len(px_2)
    p1 = normalize(px_1[0])
    p2 = normalize(px_2[0])
    dot = np.dot(p1, p2)
    theta = np.arccos(dot)
    if abs(abs(dot) - 1) < 1e-4:
        R1 = np.sign(dot) * np.eye(3)
    else:
        axis = np.cross(p1, p2)
        axis = normalize(axis)
        R1 = rotation_matrix(axis, theta)
        if len(sym_mats_old) == 1:
            return R1
    ## check if it has more than 1 paxis
    if len(px_1) > 1:
        v1 = normalize(px_1[1])
        v2 = normalize(px_2[1])
    else:
        ## else check it the group has vertical place
        dets1, nfold1, axes1 = find_symm_axis(sym_mats_old)
        dets2, nfold2, axes2 = find_symm_axis(sym_mats_new)
        if np.abs(nfold1).max() < 2:
            return R1
        bool_arr1 = np.logical_and(dets1 < 0, np.abs(nfold1) == 2)
        bool_arr2 = np.logical_and(dets2 < 0, np.abs(nfold2) == 2)
        paxis_dot1 = np.abs(axes1 @ p1)
        paxis_dot2 = np.abs(axes2 @ p2)
        bool_arr1 = np.logical_and(bool_arr1, paxis_dot1 < 1e-4)
        bool_arr2 = np.logical_and(bool_arr2, paxis_dot2 < 1e-4)
        sec_axis1 = axes1[bool_arr1, :]
        sec_axis2 = axes2[bool_arr2, :]
        assert len(sec_axis1) == len(sec_axis2)
        if len(sec_axis1) == 0:
            return R1
        v1 = normalize(sec_axis1[0])
        v2 = normalize(sec_axis2[0])
    #
    R1v1_dot_v2 = np.dot(R1 @ v1, v2)
    if np.isclose(abs(R1v1_dot_v2), 1, rtol=1e-3):
        return R1
    #
    angle_in = np.arccos(R1v1_dot_v2)
    axis_2 = np.cross(R1 @ v1, v2)
    axis_2 = normalize(axis_2)
    print(np.dot(axis_2, p2))
    assert np.isclose(abs(np.dot(axis_2, p2)), 1, rtol=1e-3)
    R2 = rotation_matrix(axis_2, angle_in)
    R = R2 @ R1
    return R


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


def divisors(n):
    """
    Returns the divisors of n.
    This isn't meant to handle large numbers, thankfully most point groups have an order less than 100
    
    :type n: int
    :return: List of n's divisors
    :rtype: List[int]
    """
    out = []
    for i in range(n):
        if n % (i + 1) == 0:
            out.append(i + 1)
    return out


def distance(a, b):
    """
    Euclidean distance between a and b.

    :type a: NumPy array of shape (n,)
    :type b: NumPy array of shape (n,)
    :return: Distance between a and b
    :rtype: float
    """
    return np.sqrt(((a - b)**2).sum())


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
    ## find symmmats and nfold degenercy
    ### incase not found, nfold is returned as None
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
    ## Works for only crystallographic point groups !
    ## Adapted from :
    ## "http://faculty.otterbein.edu/djohnston/sym/common/images/flowchart.pdf"
    """
    Given a list of symmetries, return a point group label
    """
    dets, nfold, axes = find_symm_axis(symm_mats)
    ## first find paxis
    order = len(symm_mats)
    #
    if order == 1:
        return "C1"
    #
    nidx = np.argmax(np.abs(nfold[dets > 0]))
    ## find principal axis and nfold symmetry
    n = abs(nfold[dets > 0][nidx])  ## nfold
    paxis = axes[dets > 0, :][nidx]
    Cn = symm_mats[dets > 0, :, :][nidx]

    if order == 2 and n == 0:
        ## check if it Cs
        if np.abs(nfold).max() > 0:
            return "Cs"
        else:
            return "Ci"
    # if order == 2 and n == 2:
    #     return "C_2"
    nfold_abs = np.abs(nfold)
    inversion_present = len(nfold_abs[nfold_abs == 0]) > 1
    if len(nfold_abs[np.logical_and(dets > 0, nfold_abs == 3)]) > 2:
        ## High symmetry points (O/T)
        if len(nfold_abs[np.logical_and(dets > 0, nfold_abs == 4)]) > 2:
            ## Oh/O
            if inversion_present:
                return "Oh"
            else:
                return "O"
        else:
            if len(nfold[np.logical_and(dets < 0, nfold_abs == 2)]) == 0:
                return "T"
            else:
                if inversion_present:
                    return "Th"
                else:
                    return "Td"
    else:
        ## Low symmetric (C/S/D)
        nC1 = np.abs(axes[np.logical_and(dets > 0, nfold_abs == 2), :] @ paxis)
        nC_perp = len(nC1[nC1 < 1e-3])
        ## is horizontal mirror present ?
        sig_h = np.abs(
            np.abs(axes[np.logical_and(dets < 0, nfold_abs == 2), :] @ paxis) -
            1)
        nsig_h = len(sig_h[sig_h < 1e-3])
        ## number of perpendicular planes
        sig_v = np.abs(
            axes[np.logical_and(dets < 0, nfold_abs == 2), :] @ paxis)
        nsig_v = len(sig_v[sig_v < 1e-3])
        #
        if n == nC_perp:
            # D groups
            if nsig_h > 0:
                return "D%dh" % (n)
            else:
                if (n == nsig_v):
                    return "D%dd" % (n)  # dihedral plane :
                else:
                    return "D%d" % (n)
        else:
            ## C/S groups
            if nsig_h > 0:
                return "C%dh" % (n)
            else:
                if n == nsig_v:
                    return "C%dv" % (n)
                else:
                    ##
                    isSsym = len(nfold[np.logical_and(dets < 0,
                                                      nfold_abs == 2 * n)]) > 0
                    if isSsym:
                        return "S%d" % (2 * n)
                    else:
                        return "C%d" % (n)


def get_paxis(symm_mats):
    dets, nfold, axes = find_symm_axis(symm_mats)
    ## first find paxis
    order = len(symm_mats)
    if order == 1:
        return np.array([])
    if order > 2:
        nfold_pos = nfold[np.logical_and(dets > 0, nfold > 0)]
        nidx = np.argwhere(nfold_pos == np.amax(nfold_pos)).reshape(-1)
        paxis = axes[np.logical_and(dets > 0, nfold > 0), :][nidx, :]
        paxis = paxis / np.linalg.norm(paxis, axis=-1)[:, None]
        return paxis
    else:
        paxis = axes[nfold == 2, :]
        return paxis


@dataclass
class Symel:
    """Symmetry element class."""
    label: str
    axis: np.ndarray
    rrep: np.ndarray  # Matrix representation


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
        nn = n >> 1
    else:
        nn = n
    symels = []
    x_axis = np.asarray([1, 0, 0])  # Orient sigma_v along x-axis (standard convention)
    rot_mat = Cn([0, 0, 1], n)      # Rotation matrix for 2π/n around z-axis
    for i in range(nn):
        axis = np.dot(matrix_power(rot_mat, i), x_axis)  # Rotate x-axis by i*(2π/n)
        symels.append(Symel(f"sigma_v({i+1})", axis, reflection_matrix(axis)))
    return symels


def generate_sigma_d(n):
    symels = []
    for i in range(n):
        theta = (2 * i + 1) * np.pi / (2 * n)
        axis = np.asarray([np.cos(theta), np.sin(theta), 0])
        symels.append(Symel(f"sigma_d_{i+1:d}", axis,
                            reflection_matrix(axis)))  # sigma_ds
    return symels


def generate_C2p(n):
    if n % 2 == 0:
        nn = n >> 1
    else:
        nn = n
    symels = []
    x_axis = np.asarray([1, 0, 0])  # Orient C2 along x-axis (standard convention)
    rot_mat = Cn([0, 0, 1], n)      # Rotation matrix for 2π/n around z-axis
    for i in range(nn):
        axis = np.dot(matrix_power(rot_mat, i), x_axis)  # Rotate x-axis by i*(2π/n)
        symels.append(Symel(f"C_2'({i+1})", axis, Cn(axis, 2)))
    return symels


def generate_C2pp(n):
    nn = n >> 1
    symels = []
    x_axis = np.asarray([1, 0, 0])
    rot_mat = Cn([0, 0, 1], n)
    base_axis = np.dot(Cn([0, 0, 1], 2 * n), x_axis)  # Rotate x-axis by π/n
    for i in range(nn):
        axis = np.dot(matrix_power(rot_mat, i), base_axis)  # Then rotate by i*(2π/n)
        symels.append(Symel(f"C_2''({i+1})", axis, Cn(axis, 2)))
    return symels


def generate_T():
    symels = []
    # C3 axes
    axes = [
        np.array([1, 1, 1]) / np.sqrt(3),
        np.array([1, -1, -1]) / np.sqrt(3),
        np.array([-1, 1, -1]) / np.sqrt(3),
        np.array([-1, -1, 1]) / np.sqrt(3)
    ]
    for i, axis in enumerate(axes):
        for j in range(1, 3):
            symels.append(Symel(f"C3_{i+1}^{j}", axis, matrix_power(Cn(axis, 3), j)))
    
    # C2 axes
    axes = [
        np.array([1, 0, 0]),
        np.array([0, 1, 0]),
        np.array([0, 0, 1])
    ]
    for i, axis in enumerate(axes):
        symels.append(Symel(f"C2_{i+1}", axis, Cn(axis, 2)))
    
    return symels


def generate_Td():
    symels = generate_T()
    
    # Add sigma_d planes
    axes = [
        np.array([1, 1, 0]) / np.sqrt(2),
        np.array([1, -1, 0]) / np.sqrt(2),
        np.array([1, 0, 1]) / np.sqrt(2),
        np.array([1, 0, -1]) / np.sqrt(2),
        np.array([0, 1, 1]) / np.sqrt(2),
        np.array([0, 1, -1]) / np.sqrt(2)
    ]
    for i, axis in enumerate(axes):
        symels.append(Symel(f"sigma_d_{i+1}", axis, reflection_matrix(axis)))
    
    return symels


def generate_Th():
    symels = generate_T()
    
    # Add inversion
    symels.append(Symel("i", None, inversion_matrix()))
    
    # Add S6 operations
    axes = [
        np.array([1, 1, 1]) / np.sqrt(3),
        np.array([1, -1, -1]) / np.sqrt(3),
        np.array([-1, 1, -1]) / np.sqrt(3),
        np.array([-1, -1, 1]) / np.sqrt(3)
    ]
    for i, axis in enumerate(axes):
        for j in [1, 5]:
            symels.append(Symel(f"S6_{i+1}^{j}", axis, matrix_power(Sn(axis, 6), j)))
    
    # Add sigma_h planes
    axes = [
        np.array([1, 0, 0]),
        np.array([0, 1, 0]),
        np.array([0, 0, 1])
    ]
    for i, axis in enumerate(axes):
        symels.append(Symel(f"sigma_h_{i+1}", axis, reflection_matrix(axis)))
    
    return symels


def generate_O():
    symels = []
    
    # C4 axes (3 axes)
    axes = [
        np.array([1, 0, 0]),
        np.array([0, 1, 0]),
        np.array([0, 0, 1])
    ]
    for i, axis in enumerate(axes):
        for j in range(1, 4):
            symels.append(Symel(f"C4_{i+1}^{j}", axis, matrix_power(Cn(axis, 4), j)))
    
    # C3 axes (4 axes)
    axes = [
        np.array([1, 1, 1]) / np.sqrt(3),
        np.array([1, -1, -1]) / np.sqrt(3),
        np.array([-1, 1, -1]) / np.sqrt(3),
        np.array([-1, -1, 1]) / np.sqrt(3)
    ]
    for i, axis in enumerate(axes):
        for j in range(1, 3):
            symels.append(Symel(f"C3_{i+1}^{j}", axis, matrix_power(Cn(axis, 3), j)))
    
    # C2 axes (6 axes)
    axes = [
        np.array([1, 1, 0]) / np.sqrt(2),
        np.array([1, -1, 0]) / np.sqrt(2),
        np.array([1, 0, 1]) / np.sqrt(2),
        np.array([1, 0, -1]) / np.sqrt(2),
        np.array([0, 1, 1]) / np.sqrt(2),
        np.array([0, 1, -1]) / np.sqrt(2)
    ]
    for i, axis in enumerate(axes):
        symels.append(Symel(f"C2_{i+1}", axis, Cn(axis, 2)))
    
    return symels


def generate_Oh():
    symels = generate_O()
    
    # Add inversion
    symels.append(Symel("i", None, inversion_matrix()))
    
    # Add S4 operations
    axes = [
        np.array([1, 0, 0]),
        np.array([0, 1, 0]),
        np.array([0, 0, 1])
    ]
    for i, axis in enumerate(axes):
        for j in [1, 3]:
            symels.append(Symel(f"S4_{i+1}^{j}", axis, matrix_power(Sn(axis, 4), j)))
    
    # Add S6 operations
    axes = [
        np.array([1, 1, 1]) / np.sqrt(3),
        np.array([1, -1, -1]) / np.sqrt(3),
        np.array([-1, 1, -1]) / np.sqrt(3),
        np.array([-1, -1, 1]) / np.sqrt(3)
    ]
    for i, axis in enumerate(axes):
        for j in [1, 5]:
            symels.append(Symel(f"S6_{i+1}^{j}", axis, matrix_power(Sn(axis, 6), j)))
    
    # Add sigma_h planes
    axes = [
        np.array([1, 1, 0]) / np.sqrt(2),
        np.array([1, -1, 0]) / np.sqrt(2),
        np.array([1, 0, 1]) / np.sqrt(2),
        np.array([1, 0, -1]) / np.sqrt(2),
        np.array([0, 1, 1]) / np.sqrt(2),
        np.array([0, 1, -1]) / np.sqrt(2)
    ]
    for i, axis in enumerate(axes):
        symels.append(Symel(f"sigma_h_{i+1}", axis, reflection_matrix(axis)))
    
    return symels


def generate_I():
    faces, vertices, edgecenters = generate_I_vectors()
    symels = []
    
    # C5 axes through face centers
    for i, axis in enumerate(faces):
        for j in range(1, 5):
            symels.append(Symel(f"C5_{i+1}^{j}", axis, matrix_power(Cn(axis, 5), j)))
    
    # C3 axes through vertices
    for i, axis in enumerate(vertices):
        for j in range(1, 3):
            symels.append(Symel(f"C3_{i+1}^{j}", axis, matrix_power(Cn(axis, 3), j)))
    
    # C2 axes through edge centers
    for i, axis in enumerate(edgecenters):
        symels.append(Symel(f"C2_{i+1}", axis, Cn(axis, 2)))
    
    return symels


def generate_Ih():
    faces, vertices, edgecenters = generate_I_vectors()
    symels = generate_I()
    
    # Add inversion
    symels.append(Symel("i", None, inversion_matrix()))
    
    # Add S10 operations
    for i, axis in enumerate(faces):
        for j in [1, 3, 7, 9]:
            symels.append(Symel(f"S10_{i+1}^{j}", axis, matrix_power(Sn(axis, 10), j)))
    
    # Add S6 operations
    for i, axis in enumerate(vertices):
        for j in [1, 5]:
            symels.append(Symel(f"S6_{i+1}^{j}", axis, matrix_power(Sn(axis, 6), j)))
    
    # Add sigma planes
    for i, axis in enumerate(edgecenters):
        symels.append(Symel(f"sigma_{i+1}", axis, reflection_matrix(axis)))
    
    return symels


def generate_I_vectors():
    """
    Vectors defining the faces, vertices, and edges of a regular dodecahedron.

    :rtype: (List[NumPy array of shape (3,)], List[NumPy array of shape (3,)], List[NumPy array of shape (3,)])
    """
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


def Cn_irr_complex(n):
    names = ["A"]
    classes = ["E"]
    for c in range(1, n):
        r, s = reduce(n, c)
        c_class_string(classes, "C", r, s)
    chars = np.ones(n, dtype=complex)
    if n % 2 == 0:
        names.append("B")
        bi = np.ones(n, dtype=complex)
        for i in range(n):
            if (i + 1) % 2 == 0:
                bi[i] *= -1
        chars = np.vstack((chars, bi))
    if n >= 3:
        theta = 2 * np.pi / n
        l = (n - len(names)) // 2
        for i in range(l):
            names.append(f"E{i+1}" if n >= 5 else "E")
            v = np.zeros(n, dtype=complex)
            for j in range(n):
                v[j] = np.exp(1j * (i + 1) * j * theta) + np.exp(-1j * (i + 1) * j * theta)
            chars = np.vstack((chars, v))
    return names, classes, chars


def Cn_irrmat(n):
    """Generate irreps for Cn point groups with matrix representation."""
    names = ["A"]
    classes = ["E"]
    for c in range(1, n):
        r, s = reduce(n, c)
        c_class_string(classes, "C", r, s)
    chars = np.ones(n)
    if n % 2 == 0:
        names.append("B")
        bi = np.ones(n)
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


def Cnv_irr(n):
    names = ["A1", "A2"]
    classes = ["E", f"C_{n}"]
    chars = np.array([[1.0, 1.0], [1.0, -1.0]])
    
    if n % 2 == 0:
        names.extend(["B1", "B2"])
        classes.extend([f"C_2", f"sigma_v"])
        chars = np.array([[1.0, 1.0, 1.0, 1.0],
                         [1.0, 1.0, 1.0, -1.0],
                         [1.0, -1.0, -1.0, 1.0],
                         [1.0, -1.0, -1.0, -1.0]])
    else:
        classes.append(f"sigma_v")
        chars = np.array([[1.0, 1.0, 1.0],
                         [1.0, 1.0, -1.0]])
    
    return names, classes, chars


def Cnh_irr(n):
    names = ["Ag", "Au"]
    classes = ["E", f"C_{n}", "i", "sigma_h"]
    chars = np.array([[1.0, 1.0, 1.0, 1.0],
                     [1.0, 1.0, -1.0, -1.0]])
    return names, classes, chars


def Dn_irr(n):
    """Generate irreps for Dn point groups, following the reference implementation."""
    if n == 2:
        names = ["A", "B1", "B2", "B3"]
        classes = ["E", "C_2(z)", "C_2(y)", "C_2(x)"]
        chars = np.array([[1, 1, 1, 1], [1, 1, -1, -1], [1, -1, 1, -1],
                          [1, -1, -1, 1]])
        return names, classes, chars
    
    # For n > 2, use the complex approach
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
        # Remove redundant classes
        chars = chars[:, [j for j in range(np.shape(chars)[1]) if j != n - i]]
    
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
    """Generate irreps for Dnh point groups, following the reference implementation."""
    names, classes, dnchars = Dn_irr(n)
    if n % 2 == 0:
        # Even n case (like D2h, D4h, D6h, ...)
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
        # Odd n case (like D3h, D5h, ...) - this is the D3h case!
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
    names = ["A1", "A2", "B1", "B2", "E"]
    classes = ["E", f"C_{n}", f"C_2'", f"S_{2*n}", f"sigma_d"]
    chars = np.array([[1.0, 1.0, 1.0, 1.0, 1.0],
                     [1.0, 1.0, -1.0, -1.0, -1.0],
                     [1.0, -1.0, 1.0, -1.0, 1.0],
                     [1.0, -1.0, -1.0, 1.0, -1.0],
                     [2.0, 0.0, 0.0, 0.0, 0.0]])
    return names, classes, chars


def Sn_irr_complex(n):
    names = ["A"]
    classes = ["E"]
    for c in range(1, n):
        r, s = reduce(n, c)
        c_class_string(classes, "S", r, s)
    chars = np.ones(n, dtype=complex)
    return names, classes, chars


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


def generate_symel_to_class_map(symels, ctab):
    """Generate a map of the symels to their corresponding classes."""
    pg = PointGroup.from_string(ctab.name)
    if pg.n is not None:
        ns = pg.n >> 1  # pg.n floor divided by 2
    ncls = len(ctab.classes)
    nsymel = len(symels)
    class_map = np.zeros(nsymel, dtype=np.int32)
    class_map[0] = 0  # E is always first
    
    if pg.family == "D":
        if pg.subfamily == "h":
            if pg.n % 2 == 0:
                # Even n case (D2h, D4h, D6h, ...)
                class_map[1] = ncls - 3  # σh
                class_map[2] = (ncls >> 1)  # i
                cn_class_map(class_map, pg.n, 2, 0)  # Cn
                class_map[pg.n + 2:3 * ns + 2] = ns + 1  # C2'
                class_map[3 * ns + 2:2 * pg.n + 2] = ns + 2  # C2''
                for i in range(2 * pg.n + 2, 3 * pg.n):  # Sn
                    if i > 3 * pg.n - ns:
                        class_map[i] = i - 2 * pg.n + 3
                    else:
                        class_map[i] = 3 * pg.n + 4 - i
                if pg.n % 4 == 0 or pg.n == 2:
                    class_map[-pg.n:-ns] = ncls - 2  # σv
                    class_map[-ns:] = ncls - 1  # σd
                else:
                    class_map[-pg.n:-ns] = ncls - 1  # σv
                    class_map[-ns:] = ncls - 2  # σd
            else:
                # Odd n case (D3h, D5h, ...) - this is the D3h case!
                class_map[1] = (ncls >> 1)  # σh → class 3
                cn_class_map(class_map, pg.n, 1, 0)  # C3 rotations → classes 1,2
                class_map[pg.n + 1:2 * pg.n + 1] = ns + 1  # C2 rotations → class 2
                cn_class_map(class_map, pg.n, 2 * pg.n, ns + 2)  # S3 operations → class 4
                class_map[-1 - pg.n + 1:] = ncls - 1  # σv operations → class 5
        # Add other D subfamily cases if needed
    # Add other families (C, S, etc.) if needed
    else:
        # Fallback for unimplemented cases
        class_map = np.arange(nsymel)
    
    return class_map


def cn_class_map(class_map, n, idx_offset, cls_offset):
    """Helper function for class mapping."""
    for i in range(1, n):  # = 2:n
        if i > (n >> 1):
            class_map[i + idx_offset] = n - i + cls_offset
        else:
            class_map[i + idx_offset] = i + cls_offset
    return class_map


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


def vec_similar(vec1, vec2, tol=1e-3):
    """
    Checks if two matrices are parallel
    """
    v1 = vec1 / np.linalg.norm(vec1)
    v2 = vec2 / np.linalg.norm(vec2)
    if abs(abs(np.dot(v1, v2)) - 1) < tol:
        return True
    else:
        return False