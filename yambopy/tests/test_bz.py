"""
Tests for the Brillouin-zone package (yambopy.bz).

They check geometric properties rather than stored reference numbers:
special points lie on the Wigner-Seitz BZ boundary, paths built from
equivalent cells have the same length, interpolated points are evenly
spread, and the band-structure x-axis, ticks and indices agree.
"""

import itertools
import os

import numpy as np
import pytest

from yambopy.bz.brillouin import BrillouinZone
from yambopy.bz.lattice_data import (
    _standard_lattice_data,
    find_basis_change,
    get_lattice_data,
)
from yambopy.plot.bandstructure import YambopyBandStructure

# spglib >= 2.5 warns about its legacy error handling on every call.
pytestmark = pytest.mark.filterwarnings("ignore::DeprecationWarning:spglib")

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


# -----------------------------------------------------------------------------
# Grid k-points from a yambo database (units and cell choice)
# -----------------------------------------------------------------------------

BOHR_TO_ANGSTROM = 0.529177210903
NS_DB1 = os.path.join(os.path.dirname(__file__), "..", "data", "refs", "gw", "SAVE", "ns.db1")


@pytest.fixture(scope="module")
def latticedb():
    from yambopy.dbs.latticedb import YamboLatticeDB
    return YamboLatticeDB.from_db_file(NS_DB1, Expand=False)


@pytest.mark.parametrize("expand", [True, False])
def test_collinear_grid_independent_of_units_and_cell(latticedb, expand):
    lat = np.asarray(latticedb.lat)
    rotated = (np.array([[1, 1, 0], [0, 1, 0], [0, 0, 1]]) @ lat) @ random_rotation(
        np.random.default_rng(3)).T
    cells = [lat, lat * BOHR_TO_ANGSTROM]
    if expand:
        # Only the full grid is independent of which (equivalent)
        # orientation the path takes in a differently built cell; the
        # irreducible points on the path depend on it.
        cells.append(rotated * BOHR_TO_ANGSTROM)
    results = []
    for cell in cells:
        bz = BrillouinZone(cell=cell, npoints=60)
        _, idx, d = bz.find_collinear_grid(latticedb, expand_kpoints=expand)
        results.append((idx, d / bz.special_kpoints_distances(True)[-1]))

    idx_ref, x_ref = results[0]
    assert len(idx_ref) > 0
    assert np.all(idx_ref < latticedb.ibz_nkpoints)
    for idx, x in results[1:]:
        assert np.array_equal(idx, idx_ref)
        assert np.allclose(x, x_ref)


def test_collinear_grid_full_bz_finds_more_points(latticedb):
    bz = BrillouinZone(cell=latticedb.lat, npoints=60)
    _, idx_full, _ = bz.find_collinear_grid(latticedb, expand_kpoints=True)
    _, idx_ibz, _ = bz.find_collinear_grid(latticedb, expand_kpoints=False)
    assert len(idx_full) >= len(idx_ibz)
    # The database itself is not modified by the expansion.
    assert latticedb.ibz_nkpoints == latticedb.nkpoints


def test_collinear_grid_rejects_other_lattice(latticedb):
    bz = BrillouinZone(ibrav=3, parameters={"a": 5}, npoints=20)
    with pytest.raises(ValueError, match="does not match"):
        bz.find_collinear_grid(latticedb)


# -----------------------------------------------------------------------------
# Optimisation edge cases
# -----------------------------------------------------------------------------

def test_optimise_single_point_section():
    cell = get_lattice_data(4, HEX)[0]
    bz = BrillouinZone(cell=cell, path_string="GMK,A", npoints=40,
                       optimise=True, sym_red=hexagonal_sym_red(cell))
    assert bz.star_indices[1] == [0]
    assert len(bz.kpoints()) == 40


def test_optimise_requires_symmetries():
    with pytest.raises(ValueError, match="sym_red"):
        BrillouinZone(ibrav=4, parameters=HEX, optimise=True)
    with pytest.raises(ValueError, match="sym_red"):
        BrillouinZone(cell=get_lattice_data(4, HEX)[0], optimise=True)


# -----------------------------------------------------------------------------
# Serialisation
# -----------------------------------------------------------------------------

def test_bz_json_round_trip():
    import json

    cell = get_lattice_data(4, HEX)[0]
    bz = BrillouinZone(cell=cell, path_string="GMKGQ", extra_points={"Q": np.array([0.1, 0, 0])},
                       npoints=np.int64(50), optimise=True, sym_red=hexagonal_sym_red(cell))
    clone = BrillouinZone.from_dict(json.loads(json.dumps(bz.as_dict())))
    assert np.allclose(clone.kpoints("car"), bz.kpoints("car"))
    assert clone.path_labels_list(True) == bz.path_labels_list(True)


def test_bandstructure_json_round_trip(tmp_path):
    bz = BrillouinZone(ibrav=4, parameters=HEX, npoints=60)
    nk = len(bz.kpoints())
    bands = np.linspace(0, 1, nk * 2).reshape(nk, 2)
    original = YambopyBandStructure(bands, bz.kpoints("car"), bz=bz, fermie=0.5)
    path = tmp_path / "bands.json"
    original.write_json(str(path))

    loaded = YambopyBandStructure.from_json(str(path))
    assert np.allclose(loaded.bands, bands)
    assert loaded.fermie == 0.5
    assert np.allclose(loaded.distances, original.distances)
    assert loaded.bz.path_labels_list(True) == bz.path_labels_list(True)

    no_bz = YambopyBandStructure(bands, bz.kpoints("car"))
    no_bz.write_json(str(path))
    assert YambopyBandStructure.from_json(str(path)).bz is None


# -----------------------------------------------------------------------------
# Band structure from a yambo electrons database
# -----------------------------------------------------------------------------

@pytest.fixture(scope="module")
def electronsdb():
    from yambopy.dbs.electronsdb import YamboElectronsDB
    return YamboElectronsDB.from_db_file(os.path.dirname(NS_DB1), Expand=False)


def test_electrons_path_matches_lattice_grid(electronsdb, latticedb):
    lat = np.asarray(electronsdb.lat)
    bz = BrillouinZone(cell=lat, npoints=60)
    x, ibz = electronsdb.get_bandstructure_path(bz)
    _, ibz_ref, x_ref = bz.find_collinear_grid(latticedb, expand_kpoints=True)
    assert len(ibz) > 0
    assert np.array_equal(ibz, ibz_ref)
    assert np.allclose(x, x_ref)

    for cell in (lat * BOHR_TO_ANGSTROM,
                 (np.array([[1, 1, 0], [0, 1, 0], [0, 0, 1]]) @ lat) * BOHR_TO_ANGSTROM):
        other = BrillouinZone(cell=cell, npoints=60)
        x_other, ibz_other = electronsdb.get_bandstructure_path(other)
        assert np.array_equal(ibz_other, ibz)
        assert np.allclose(x_other / other.special_kpoints_distances(True)[-1],
                           x / bz.special_kpoints_distances(True)[-1])

    # The database is left unexpanded.
    assert electronsdb.expanded is False


def test_electrons_plot_bandstructure(electronsdb):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    bz = BrillouinZone(cell=np.asarray(electronsdb.lat) * BOHR_TO_ANGSTROM, npoints=60)
    fig, ax = plt.subplots()
    electronsdb.plot_bandstructure_ax(ax, bz, bandmin=0, bandmax=4)
    assert np.allclose(ax.get_xticks(), bz.special_kpoints_distances(True))
    assert [t.get_text() for t in ax.get_xticklabels()] == bz.path_labels_list(True, fmt="latex")
    assert len(ax.lines) >= 4
    plt.close(fig)


# -----------------------------------------------------------------------------
# KPathOptimizer tie-breaking
# -----------------------------------------------------------------------------

def test_turn_penalty_favours_straight_paths():
    from yambopy.bz.kpath_optimizer import KPathOptimizer

    k0, k1 = np.array([0.0, 0, 0]), np.array([1.0, 0, 0])
    assert KPathOptimizer._turn_penalty(k0, k1, np.array([2.0, 0, 0]), 1.0) == pytest.approx(0.0)
    assert KPathOptimizer._turn_penalty(k0, k1, np.array([1.0, 1, 0]), 1.0) == pytest.approx(0.05)
    assert KPathOptimizer._turn_penalty(k0, k1, np.array([0.0, 0, 0]), 1.0) == pytest.approx(0.1)


@pytest.mark.parametrize("method", ["greedy", "greedy_global"])
def test_greedy_breaks_ties_by_carrying_on(method):
    from yambopy.bz.kpath_optimizer import KPathOptimizer

    spglib = pytest.importorskip("spglib")
    cell = np.diag([3.0, 3.0, 10.0])  # square in-plane: (1/4, 0) has 4 equidistant copies
    sym = [w.T for w in spglib.get_symmetry((cell, [[0, 0, 0]], [1]))["rotations"]]
    opt = KPathOptimizer(np.linalg.inv(cell.T), sym)
    # X -> G -> P: every copy of P is 1/4 from G; carrying on means P = (-1/4, 0).
    result = opt.optimize_kpath([[0.5, 0, 0], [0, 0, 0], [0.25, 0, 0]], method=method)
    assert np.allclose(result["kpath_frac"][1:], [[0, 0, 0], [-0.25, 0, 0]])


def test_cell_to_ibrav_removed():
    import yambopy
    assert not hasattr(yambopy, "CellToIbrav")


# -----------------------------------------------------------------------------
# Path strings and labels
# -----------------------------------------------------------------------------

from yambopy.bz.labels import (  # noqa: E402
    format_label,
    label_to_latex,
    label_to_unicode,
    parse_label,
    parse_path,
)


@pytest.mark.parametrize("path, sections", [
    ("GMKG", [["G", "M", "K", "G"]]),
    ("GMKGALHA,LM,KH", [list("GMKGALHA"), ["L", "M"], ["K", "H"]]),
    ("G-M-K-G", [["G", "M", "K", "G"]]),
    ("Gamma M K Gamma", [["G", "M", "K", "G"]]),
    ("GammaMKGamma", [["G", "M", "K", "G"]]),
    ("GKK'G", [["G", "K", "K'", "G"]]),
    ("GXYSigmaGZSigma1NPY1Z,XP", [["G", "X", "Y", "Sigma", "G", "Z", "Sigma1", "N", "P", "Y1", "Z"], ["X", "P"]]),
    ("PhiPiPsiP", [["Phi", "Pi", "Psi", "P"]]),
    ("GM,", [["G", "M"]]),
])
def test_parse_path(path, sections):
    assert parse_path(path) == sections


@pytest.mark.parametrize("path", ["gMK", "GKpG", "G#M", "GΣM", "Gam", "Alpha"])
def test_parse_path_rejects_invalid(path):
    with pytest.raises(ValueError, match="Invalid path"):
        parse_path(path)


def test_parse_label():
    assert parse_label("Sigma1") == ("Sigma", "1", 0)
    assert parse_label("K''") == ("K", "", 2)
    with pytest.raises(ValueError, match="Invalid special-point label"):
        parse_label("Kp")


@pytest.mark.parametrize("label, latex, unicode", [
    ("G", r"$\Gamma$", "Γ"),
    ("M", "$M$", "M"),
    ("X1", "$X_{1}$", "X₁"),
    ("Sigma1", r"$\Sigma_{1}$", "Σ₁"),
    ("K'", r"$K^{\prime}$", "K′"),
    ("G2'", r"$\Gamma_{2}^{\prime}$", "Γ₂′"),
])
def test_label_rendering(label, latex, unicode):
    assert label_to_latex(label) == latex
    assert label_to_unicode(label) == unicode


def test_format_merged_and_legacy_labels():
    assert format_label("A — L", "plain") == "A — L"
    assert format_label("A — L", "latex") == "$A|L$"
    assert format_label("H — G", "unicode") == "H|Γ"
    assert format_label(r"$\Gamma$", "latex") == r"$\Gamma$"
    with pytest.raises(ValueError, match="Unknown label format"):
        format_label("G", "html")


def test_latex_labels_render_in_matplotlib():
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.mathtext import MathTextParser

    parser = MathTextParser("path")
    bz = BrillouinZone(ibrav=7, parameters={"a": 3, "c": 5}, npoints=50)
    for label in bz.path_labels_list(True, fmt="latex"):
        parser.parse(label)  # raises on invalid mathtext


def test_bct2_uses_sigma():
    bz = BrillouinZone(ibrav=7, parameters={"a": 3, "c": 5}, npoints=50)
    assert bz.variant == "BCT2"
    assert {"Sigma", "Sigma1"} <= set(bz.high_symmetry_points)
    assert "S" not in bz.high_symmetry_points
    assert bz.path_labels_list(True, fmt="unicode")[:5] == ["Γ", "X", "Y", "Σ", "Γ"]


def test_path_aliases_and_extra_points():
    bz = BrillouinZone(ibrav=4, parameters=HEX, path_string="Gamma-K-K'-Gamma",
                       extra_points={"K'": [-1 / 3, -1 / 3, 0]}, npoints=30)
    assert bz.path_sections == [["G", "K", "K'", "G"]]
    assert bz.path_string == "GKK'G"
    assert bz.path_labels_list(fmt="latex") == [r"$\Gamma$", "$K$", r"$K^{\prime}$", r"$\Gamma$"]


def test_undefined_label_is_reported():
    with pytest.raises(ValueError, match=r"undefined special points \['Q'\]"):
        BrillouinZone(ibrav=4, parameters=HEX, path_string="GMQ")


def test_single_point_sections_merge_consistently():
    bz = BrillouinZone(ibrav=4, parameters=HEX, path_string="GMK,A,LH", npoints=60)
    labels = bz.path_labels_list(True)
    assert labels == ["G", "M", "K — A — L", "H"]
    assert [label for _, label in bz.get_indices()] == labels
    assert len(bz.special_kpoints_distances(True)) == len(labels)
    assert [label for _, label in bz.get_indices(fmt="latex")] == bz.path_labels_list(True, fmt="latex")


def test_bandstructure_ticks_are_latex():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    bz = BrillouinZone(ibrav=4, parameters=HEX, npoints=60)
    nk = len(bz.kpoints())
    fig, ax = plt.subplots()
    YambopyBandStructure(np.linspace(0, 1, 2 * nk).reshape(nk, 2), bz.kpoints("car"), bz=bz).plot_ax(ax)
    assert [t.get_text() for t in ax.get_xticklabels()] == bz.path_labels_list(True, fmt="latex")
    plt.close(fig)


@pytest.mark.filterwarnings("ignore:FigureCanvasAgg is non-interactive")
def test_qepy_plots_use_latex_labels():
    import types
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from qepy.pwxml import PwXML
    from qepy.matdyn import Matdyn

    bz = BrillouinZone(ibrav=4, parameters=HEX, npoints=60)
    nk = len(bz.kpoints())
    expected = bz.path_labels_list(True, fmt="latex")

    pw = types.SimpleNamespace(kpoints=bz.kpoints("car"), lsda=False,
                               eigen1=np.zeros((nk, 2)), nbands=2)
    fig, ax = plt.subplots()
    PwXML.plot_eigen_ax(pw, ax, path_kpoints=bz.get_indices())
    assert [t.get_text() for t in ax.get_xticklabels()] == expected
    plt.close(fig)

    phonons = types.SimpleNamespace(nqpoints=nk, nmodes=3, eig=np.zeros((nk, 3)))
    Matdyn.plot_eigen(phonons, bz)
    ax = plt.gca()
    assert [t.get_text() for t in ax.get_xticklabels()] == expected
    assert np.allclose(ax.lines[-1].get_xdata(), bz.kpoints_distances())
    plt.close("all")
    with pytest.raises(ValueError, match="q-points"):
        Matdyn.plot_eigen(types.SimpleNamespace(nqpoints=nk - 1, nmodes=3, eig=np.zeros((nk - 1, 3))), bz)


# -----------------------------------------------------------------------------
# Cells identical to a standard one keep its special points
# -----------------------------------------------------------------------------

@pytest.mark.parametrize("ibrav, parameters, variant", LATTICES, ids=LATTICE_IDS)
def test_standard_cell_keeps_tabulated_points(ibrav, parameters, variant):
    from yambopy.bz.lattice_data import lattice_data_for_cell

    cell, _, points, _ = get_lattice_data(ibrav, parameters)
    _, var, points_cell, _ = lattice_data_for_cell(cell)
    assert var == variant
    for label, k in points.items():
        assert np.allclose(points_cell[label], k), label


# -----------------------------------------------------------------------------
# Path optimisation options
# -----------------------------------------------------------------------------

FAR_K = {"K": [-1 / 3, -1 / 3, 0], "H": [-1 / 3, -1 / 3, 0.5]}  # corners not adjacent to M, L


@pytest.fixture(scope="module")
def hex_cell_and_sym():
    cell = get_lattice_data(4, HEX)[0]
    return cell, hexagonal_sym_red(cell)


def test_optimise_options_normalisation():
    from yambopy.bz.brillouin import OPTIMISE_DEFAULTS, _optimise_options

    assert _optimise_options(False) is None
    assert _optimise_options(None) is None
    assert _optimise_options(True) == OPTIMISE_DEFAULTS
    assert _optimise_options({"method": "brute"}) == {**OPTIMISE_DEFAULTS, "method": "brute"}
    assert _optimise_options({"max_combinations": np.int64(10)})["max_combinations"] == 10


@pytest.mark.parametrize("optimise, error, match", [
    ({"methd": "brute"}, ValueError, "Unknown optimise options"),
    ({"method": "exact"}, ValueError, "Unknown optimise method"),
    ({"max_combinations": 0}, ValueError, "max_combinations"),
    ({"max_combinations": True}, ValueError, "max_combinations"),
    ({"parallel_segments": "yes"}, ValueError, "parallel_segments"),
    ("brute", TypeError, "bool or a dict"),
])
def test_optimise_options_rejected(hex_cell_and_sym, optimise, error, match):
    cell, sym = hex_cell_and_sym
    with pytest.raises(error, match=match):
        BrillouinZone(cell=cell, npoints=20, optimise=optimise, sym_red=sym)


def test_optimise_method_is_used(hex_cell_and_sym):
    cell, sym = hex_cell_and_sym
    common = dict(cell=cell, path_string="GMKGALHA", extra_points=FAR_K, npoints=100, sym_red=sym)
    greedy = BrillouinZone(**common, optimise={"method": "greedy"})
    brute = BrillouinZone(**common, optimise={"method": "brute"})
    assert brute.optimisation_results[0]["improvement"] > greedy.optimisation_results[0]["improvement"] > 0
    assert brute.star_indices != greedy.star_indices
    assert brute.as_dict()["optimise"] == {"method": "brute", "max_combinations": 100000,
                                           "parallel_segments": True}


def test_optimise_brute_fallback_warns(hex_cell_and_sym):
    cell, sym = hex_cell_and_sym
    with pytest.warns(UserWarning, match="Too many combinations"):
        bz = BrillouinZone(cell=cell, path_string="GMKG", extra_points=FAR_K, npoints=50,
                           optimise={"method": "brute", "max_combinations": 1}, sym_red=sym)
    assert "warning" in bz.optimisation_results[0]


def test_not_optimised_has_no_results():
    bz = BrillouinZone(ibrav=4, parameters=HEX, npoints=20)
    assert bz.optimisation_results is None
    assert bz.as_dict()["optimise"] is False


@pytest.mark.parametrize("mode", [{"npoints": 80}, {"density": 50}, {"intervals": [10, 7, 12]}])
def test_optimise_path_after_construction(hex_cell_and_sym, mode):
    cell, sym = hex_cell_and_sym
    common = dict(cell=cell, path_string="GMKG", extra_points=FAR_K, **mode)
    at_init = BrillouinZone(**common, optimise={"method": "brute"}, sym_red=sym)

    later = BrillouinZone(**common)
    before = later.kpoints("car").copy()
    results = later.optimise_path(sym, method="brute")

    assert results is later.optimisation_results and results[0]["improvement"] > 0
    assert later.star_indices == at_init.star_indices
    assert np.allclose(later.kpoints("car"), at_init.kpoints("car"))
    assert not np.allclose(later.kpoints("car")[-3:], before[-3:]) or len(before) != len(later.kpoints())
    if "intervals" in mode:
        assert later.intervals == mode["intervals"]


def test_optimise_path_survives_serialisation(hex_cell_and_sym):
    import json

    cell, sym = hex_cell_and_sym
    bz = BrillouinZone(cell=cell, path_string="GMKG", extra_points=FAR_K, npoints=60)
    bz.optimise_path(sym, method="greedy", parallel_segments=False)
    clone = BrillouinZone.from_dict(json.loads(json.dumps(bz.as_dict())))
    assert clone.as_dict()["optimise"] == {"method": "greedy", "max_combinations": 100000,
                                           "parallel_segments": False}
    assert clone.star_indices == bz.star_indices
    assert np.allclose(clone.kpoints("car"), bz.kpoints("car"))


def test_optimise_path_reoptimises_from_original_points(hex_cell_and_sym):
    cell, sym = hex_cell_and_sym
    bz = BrillouinZone(cell=cell, path_string="GMKGALHA", extra_points=FAR_K, npoints=100)
    bz.optimise_path(sym, method="brute")
    brute_indices = [list(section) for section in bz.star_indices]
    bz.optimise_path(sym, method="greedy")
    fresh = BrillouinZone(cell=cell, path_string="GMKGALHA", extra_points=FAR_K, npoints=100,
                          optimise={"method": "greedy"}, sym_red=sym)
    assert bz.star_indices == fresh.star_indices != brute_indices
    assert np.allclose(bz.kpoints("car"), fresh.kpoints("car"))


def test_optimise_path_needs_lattice(hex_cell_and_sym):
    _, sym = hex_cell_and_sym
    bz = BrillouinZone(path_string="GMKG", intervals=[5, 5, 5])
    with pytest.raises(ValueError, match="needs a lattice"):
        bz.optimise_path(sym)


# -----------------------------------------------------------------------------
# Library methods taking a BrillouinZone instead of the old Path
# -----------------------------------------------------------------------------

import types  # noqa: E402

FCC_PATH = "GXWKGLUWLK,UX"  # two sections, to exercise the merged tick


def _agg():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def test_bandstructure_on_grid_points_uses_path_positions(latticedb):
    bz = BrillouinZone(cell=latticedb.lat, path_string=FCC_PATH, npoints=60)
    kpts, _, d = bz.find_collinear_grid(latticedb, expand_kpoints=True)
    bands = YambopyBandStructure(np.zeros((len(kpts), 2)), kpts, bz=bz)
    assert np.allclose(bands.distances, d)
    # The first k-point (Gamma) sits on the first tick, and no position
    # goes beyond the last tick (section jumps are not counted).
    ticks = bz.special_kpoints_distances(True)
    assert bands.distances[0] == pytest.approx(ticks[0])
    assert max(bands.distances) <= ticks[-1] + 1e-9


def test_exciton_dispersion_on_bz(latticedb):
    from copy import deepcopy
    from yambopy.bse.bse_dispersion import ExcitonDispersion

    lattice = deepcopy(latticedb)
    lattice.expand_kpoints(verbose=0)
    exc_energies = np.arange(lattice.ibz_nkpoints, dtype=float)[:, None] * [1.0, 2.0]
    disp = types.SimpleNamespace(red_qpoints=lattice.red_kpoints, lattice=lattice, exc_energies=exc_energies)

    bz = BrillouinZone(cell=np.asarray(lattice.lat) * BOHR_TO_ANGSTROM, path_string=FCC_PATH, npoints=60)
    ybs = ExcitonDispersion.get_dispersion(disp, bz)
    _, ibz_ref, d_ref = bz.find_collinear_grid(latticedb, expand_kpoints=True)

    assert ybs.bz is bz
    assert np.allclose(ybs.bands[:, 0], ibz_ref)   # energy rows follow the IBZ index of each q
    assert np.allclose(ybs.distances, d_ref)


def test_exciton_weight_plot_on_bz(latticedb):
    from yambopy.bse.excitonweight import YamboExcitonWeight

    plt = _agg()
    lat = np.asarray(latticedb.lat)
    kpts = np.asarray(latticedb.car_kpoints)  # IBZ points, yambopy convention
    nk = len(kpts)
    weight = types.SimpleNamespace(
        lat=lat,
        calc_kpts_transitions=lambda **kw: (kpts, [np.ones(nk)], np.arange(1, nk + 1)),
        eigenvalues=np.tile(np.arange(4.0), (nk, 1)),
        nbands=4,
        transitions_v_to_c={(1, 3): 1.0},
    )
    bz = BrillouinZone(cell=lat * BOHR_TO_ANGSTROM, path_string=FCC_PATH, npoints=60)
    fig, ax = plt.subplots()
    YamboExcitonWeight.plot_exciton_bs(weight, ax, bz)
    _, _, d = bz.find_collinear_lattice_kpoints(lat, kpts)
    assert np.allclose(ax.lines[-1].get_xdata(), d)
    assert [t.get_text() for t in ax.get_xticklabels()] == bz.path_labels_list(True, fmt="latex")
    plt.close(fig)


def test_unfolding_plot_uses_bz_indices():
    from qepy.unfolding import Unfolding

    plt = _agg()
    bz = BrillouinZone(path_string="MKG", intervals=[10, 12])  # labels only, as for a supercell
    nk = len(bz.kpoints_distances()) if bz.rcell is not None else sum(bz.intervals) + 1
    unf = types.SimpleNamespace(kpoints=np.c_[np.arange(nk), np.zeros(nk), np.zeros(nk)],
                                nbands_pc=1, eigen_pc=np.zeros((nk, 1)), band_min=0,
                                nbands_sc=1, eigen_sc=np.zeros((nk, 1)), projection=np.ones((nk, 1)))
    fig, ax = plt.subplots()
    Unfolding.plot_eigen_ax(unf, ax, path=bz.get_indices())
    assert [t.get_text() for t in ax.get_xticklabels()] == bz.path_labels_list(True, fmt="latex")
    plt.close(fig)


def test_spin_polarised_exciton_bs_passes_bz():
    from yambopy.dbs.excitondb import YamboExcitonDB

    bz = BrillouinZone(ibrav=4, parameters=HEX, npoints=30)
    nk = 5
    calls = []

    def fake_bs(energies, bz_arg, excitons, debug):
        calls.append(bz_arg)
        e = np.zeros((nk, 4))
        return bz.kpoints("car")[:nk], e, e, e + 1, e + 1

    exc = types.SimpleNamespace(spin_pol="pol", start_band=0, mband=4, exciton_bs_spin_pol=fake_bs)
    up, dw = YamboExcitonDB.get_exciton_bs_spin_pol(exc, None, bz, (0,))
    assert calls == [bz] and up.bz is bz and dw.bz is bz


def test_exciton_bs_plot_wrappers_pass_bz():
    from yambopy.dbs.excitondb import YamboExcitonDB

    plt = _agg()
    bz = BrillouinZone(ibrav=4, parameters=HEX, npoints=30)
    nk = len(bz.kpoints())
    received = {}

    def fake_get_exciton_bs(energies_db, bz_arg, excitons, size=1, f=None, debug=False):
        received.update(bz=bz_arg, size=size)
        return YambopyBandStructure(np.linspace(0, 1, 2 * nk).reshape(nk, 2), bz.kpoints("car"), bz=bz_arg)

    exc = types.SimpleNamespace(get_exciton_bs=fake_get_exciton_bs)
    exc.plot_exciton_bs_ax = lambda *a, **k: YamboExcitonDB.plot_exciton_bs_ax(exc, *a, **k)
    fig, ax = plt.subplots()
    YamboExcitonDB.plot_exciton_bs_ax(exc, ax, None, bz, (0,), size=2)
    assert received == {"bz": bz, "size": 2}
    fig = YamboExcitonDB.plot_exciton_bs(exc, None, bz, (0,), show=False)
    assert received["bz"] is bz
    plt.close("all")


def test_path_distances():
    bz = BrillouinZone(ibrav=4, parameters=HEX, path_string="GMKG,AL", npoints=60)
    G, M, K, A, L = (np.asarray(bz.high_symmetry_points[x][0]) @ bz.rcell for x in "GMKAL")
    ticks = bz.special_kpoints_distances(True)  # G, M, K, G—A, L

    # Gamma visited twice, then the next section's start on the merged tick
    x = bz.path_distances(np.array([G, M, K, G, A, L]))
    assert np.allclose(x, ticks[[0, 1, 2, 3, 3, 4]])
    # Equivalent copies (shifted by a reciprocal lattice vector) are found too
    shift = bz.rcell[0] + bz.rcell[1]
    assert np.allclose(bz.path_distances(np.array([G + shift, M])), ticks[:2])
    # A point off the path, or out of path order
    assert bz.path_distances(np.array([0.3 * bz.rcell[2]])) is None
    assert bz.path_distances(np.array([K, M])) is None
