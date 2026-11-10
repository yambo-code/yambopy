"""
Tests for the Brillouin-zone package (yambopy.bz).

They check geometric properties rather than stored reference numbers:
special points lie on the Wigner-Seitz BZ boundary, paths built from
equivalent cells have the same length, interpolated points are evenly
spread, and the band-structure x-axis, ticks and indices agree.
"""

import itertools

import numpy as np
import pytest

from yambopy.bz.brillouin import BrillouinZone
from yambopy.bz.lattice_data import (
    _standard_lattice_data,
    find_basis_change,
    get_lattice_data,
)
from yambopy.plot.bandstructure import YambopyBandStructure

# (ibrav, parameters, expected variant).  Unordered constants for ibrav
# 9, -9, 10 and 11 exercise the reordering in get_lattice_data.
LATTICES = [
    (1, {"a": 3}, "CUB"),
    (2, {"a": 3}, "FCC"),
    (3, {"a": 3}, "BCC"),
    (-3, {"a": 3}, "BCC"),
    (4, {"a": 3, "c": 5}, "HEX"),
    (5, {"a": 3, "gamma": 70}, "RHL1"),
    (5, {"a": 3, "gamma": 100}, "RHL2"),
    (-5, {"a": 3, "gamma": 70}, "RHL1"),
    (-5, {"a": 3, "gamma": 100}, "RHL2"),
    (6, {"a": 3, "c": 5}, "TET"),
    (7, {"a": 3, "c": 2}, "BCT1"),
    (7, {"a": 3, "c": 5}, "BCT2"),
    (8, {"a": 3, "b": 4, "c": 5}, "ORC"),
    (8, {"a": 5, "b": 3, "c": 4}, "ORC"),
    (9, {"a": 3, "b": 4, "c": 5}, "ORCC"),
    (9, {"a": 4, "b": 3, "c": 5}, "ORCC"),
    (-9, {"a": 3, "b": 4, "c": 5}, "ORCC"),
    (-9, {"a": 4, "b": 3, "c": 5}, "ORCC"),
    (10, {"a": 3, "b": 4, "c": 5}, "ORCF1"),
    (10, {"a": 5, "b": 4, "c": 3}, "ORCF1"),
    (10, {"a": 4, "b": 5, "c": 3}, "ORCF1"),
    (10, {"a": 3, "b": 3.5, "c": 4}, "ORCF2"),
    (10, {"a": 4, "b": 3, "c": 3.5}, "ORCF2"),
    (11, {"a": 3, "b": 4, "c": 5}, "ORCI"),
    (11, {"a": 5, "b": 4, "c": 3}, "ORCI"),
    (11, {"a": 4, "b": 5, "c": 3}, "ORCI"),
]
LATTICE_IDS = [f"ibrav{ib}-" + "-".join(f"{k}{v}" for k, v in p.items()) for ib, p, _ in LATTICES]

HEX = {"a": 3, "c": 5}
MULTI_SECTION = [
    (4, HEX),
    (2, {"a": 5}),
    (10, {"a": 3, "b": 4, "c": 5}),
    (1, {"a": 3}),
]

_G = np.array(list(itertools.product(range(-4, 5), repeat=3)))
_G = _G[np.any(_G != 0, axis=1)]


def bz_boundary_violations(cell, points):
    """Labels of the non-Gamma points not on the Wigner-Seitz BZ boundary."""
    rcell = np.linalg.inv(np.asarray(cell, dtype=float).T)
    gvecs = _G @ rcell
    bad = []
    for label, k in points.items():
        if label == "G":
            continue
        kc = np.asarray(k, dtype=float) @ rcell
        nearest = np.min(np.linalg.norm(kc - gvecs, axis=1))
        if abs(nearest - np.linalg.norm(kc)) > 1e-6:
            bad.append(label)
    return bad


def random_rotation(rng):
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    return q * np.sign(np.linalg.det(q))


def bz_points(bz):
    """Special points of a BrillouinZone (first star member)."""
    return {label: v[0] for label, v in bz.high_symmetry_points.items()}


# -----------------------------------------------------------------------------
# get_lattice_data
# -----------------------------------------------------------------------------

@pytest.mark.parametrize("ibrav, parameters, variant", LATTICES, ids=LATTICE_IDS)
def test_special_points_on_bz_boundary(ibrav, parameters, variant):
    cell, var, points, path = get_lattice_data(ibrav, parameters)
    assert var == variant
    assert bz_boundary_violations(cell, points) == []


@pytest.mark.parametrize("ibrav, parameters, variant", LATTICES, ids=LATTICE_IDS)
def test_cell_follows_qe_definition(ibrav, parameters, variant):
    """Reordering the constants must not change the returned cell."""
    cell = get_lattice_data(ibrav, parameters)[0]
    assert np.allclose(cell, _standard_lattice_data(ibrav, parameters)[0])


def test_find_basis_change():
    cell = get_lattice_data(2, {"a": 3})[0]
    U = np.array([[1, 1, 0], [0, 1, 0], [-1, 0, 1]])
    rot = random_rotation(np.random.default_rng(1))
    M = find_basis_change(cell, U @ cell @ rot.T)
    assert M is not None and round(abs(np.linalg.det(M))) == 1
    G_std = cell @ cell.T
    assert np.allclose(M @ G_std @ M.T, U @ G_std @ U.T)

    other = get_lattice_data(3, {"a": 3})[0]
    assert find_basis_change(cell, other) is None


# -----------------------------------------------------------------------------
# BrillouinZone from a cell
# -----------------------------------------------------------------------------

@pytest.mark.parametrize("transform", ["as_is", "rotated", "rebased"])
@pytest.mark.parametrize("ibrav, parameters, variant", LATTICES, ids=LATTICE_IDS)
def test_cell_input(ibrav, parameters, variant, transform):
    cell = get_lattice_data(ibrav, parameters)[0]
    if transform == "rotated":
        cell = cell @ random_rotation(np.random.default_rng(abs(ibrav))).T
    elif transform == "rebased":
        cell = np.array([[1, 1, 0], [0, 1, 0], [0, 0, 1]]) @ cell

    bz = BrillouinZone(cell=cell, npoints=50)
    ref = BrillouinZone(ibrav=ibrav, parameters=parameters, npoints=50)

    assert np.allclose(bz.cell, cell)
    assert bz.variant == variant
    assert bz_boundary_violations(bz.cell, bz_points(bz)) == []
    assert np.isclose(
        bz.special_kpoints_distances(True)[-1],
        ref.special_kpoints_distances(True)[-1],
    )


def test_hexagonal_60_degrees():
    a, c = HEX["a"], HEX["c"]
    cell60 = np.array([[a, 0, 0], [a / 2, a * np.sqrt(3) / 2, 0], [0, 0, c]])
    bz = BrillouinZone(cell=cell60, path_string="GMKG", npoints=30)
    ref = BrillouinZone(ibrav=4, parameters=HEX, path_string="GMKG", npoints=30)

    assert bz.variant == "HEX"
    assert np.allclose(bz.cell, cell60)
    assert bz_boundary_violations(cell60, bz_points(bz)) == []
    assert np.allclose(bz.special_kpoints_distances(True), ref.special_kpoints_distances(True))


def test_unsupported_cell_raises():
    monoclinic = np.array([[3.0, 0, 0], [0, 4.0, 0], [1.0, 0, 5.0]])
    with pytest.raises(ValueError, match="not supported"):
        BrillouinZone(cell=monoclinic)


# -----------------------------------------------------------------------------
# Interpolation
# -----------------------------------------------------------------------------

def segment_lengths(bz):
    return np.array([
        np.linalg.norm(s[i + 1] - s[i])
        for s in bz.special_kpoints("car")
        for i in range(len(s) - 1)
    ])


def test_density_is_uniform():
    density = 200
    bz = BrillouinZone(ibrav=4, parameters=HEX, density=density)
    lengths = segment_lengths(bz)
    per_length = np.array(bz.intervals) / lengths
    assert np.allclose(per_length, density, rtol=0.1)
    assert abs(len(bz.kpoints()) - lengths.sum() * density) <= len(bz.path_sections) + 1


@pytest.mark.parametrize("npoints", [50, 100, 301])
def test_npoints_is_exact(npoints):
    bz = BrillouinZone(ibrav=2, parameters={"a": 5}, npoints=npoints)
    assert len(bz.kpoints()) == npoints


def test_intervals():
    bz = BrillouinZone(ibrav=4, parameters=HEX, path_string="GMK,AL", intervals=[3, 4, 5])
    # 3 + 4 points + K closes section 1; 5 points + L closes section 2.
    assert len(bz.kpoints()) == 3 + 4 + 1 + 5 + 1
    assert bz.intervals == [3, 4, 5]
    with pytest.raises(ValueError):
        BrillouinZone(ibrav=4, parameters=HEX, path_string="GMK", intervals=[3])


# -----------------------------------------------------------------------------
# Band-structure x-axis, ticks and indices
# -----------------------------------------------------------------------------

@pytest.mark.parametrize("ibrav, parameters", MULTI_SECTION)
def test_axis_ticks_and_indices_agree(ibrav, parameters):
    bz = BrillouinZone(ibrav=ibrav, parameters=parameters, npoints=120)
    kpts = bz.kpoints()
    x = bz.kpoints_distances()
    ticks = bz.special_kpoints_distances(True)
    indices = bz.get_indices()

    assert len(bz.path_sections) > 1
    assert len(x) == len(kpts)
    assert [label for _, label in indices] == bz.path_labels_list(True)
    assert np.allclose([x[i] for i, _ in indices], ticks)
    for i, label in indices:
        first = label.split(" — ")[0]
        assert np.allclose(kpts[i], bz.high_symmetry_points[first][0])


@pytest.mark.parametrize("ibrav, parameters", MULTI_SECTION)
def test_bandstructure_distances_match_bz(ibrav, parameters):
    bz = BrillouinZone(ibrav=ibrav, parameters=parameters, npoints=120)
    nk = len(bz.kpoints())
    bands = YambopyBandStructure(np.zeros((nk, 2)), bz.kpoints("car"), bz=bz)
    assert np.allclose(bands.distances, bz.kpoints_distances())


# -----------------------------------------------------------------------------
# Path optimisation and collinear k-points
# -----------------------------------------------------------------------------

def hexagonal_sym_red(cell):
    """Point-group operations in yambopy's sym_red convention (x' = x @ S)."""
    spglib = pytest.importorskip("spglib")
    rotations = spglib.get_symmetry((cell, [[0, 0, 0]], [1]))["rotations"]
    return [w.T for w in rotations]


@pytest.mark.filterwarnings("ignore::DeprecationWarning:spglib")
def test_optimise_on_cell():
    a, c = HEX["a"], HEX["c"]
    cell60 = np.array([[a, 0, 0], [a / 2, a * np.sqrt(3) / 2, 0], [0, 0, c]])
    bz = BrillouinZone(cell=cell60, path_string="GMKGALHA,LM,KH", npoints=100,
                       optimise=True, sym_red=hexagonal_sym_red(cell60))
    ref = BrillouinZone(ibrav=4, parameters=HEX, npoints=100)

    stars = {label: len(v) for label, v in bz.high_symmetry_points.items()}
    assert stars == {"G": 1, "A": 1, "H": 2, "K": 2, "L": 3, "M": 3}
    assert len(bz.kpoints()) == 100
    assert np.isclose(bz.special_kpoints_distances(True)[-1],
                      ref.special_kpoints_distances(True)[-1])


def test_collinear_kpoints_independent_of_cell_choice():
    a, c = HEX["a"], HEX["c"]
    cell60 = np.array([[a, 0, 0], [a / 2, a * np.sqrt(3) / 2, 0], [0, 0, c]])
    bz = BrillouinZone(cell=cell60, npoints=50)
    ref = BrillouinZone(ibrav=4, parameters=HEX, npoints=50)
    grid = np.array([[i / 6, j / 6, k / 2] for i in range(6) for j in range(6) for k in range(2)])

    _, _, d = bz.find_collinear_kpoints(grid @ bz.rcell)
    _, _, d_ref = ref.find_collinear_kpoints(grid @ ref.rcell)
    assert len(d) > 0
    assert np.allclose(d, d_ref)
