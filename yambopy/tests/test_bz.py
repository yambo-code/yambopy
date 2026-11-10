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
    results = []
    for cell in (lat, lat * BOHR_TO_ANGSTROM, rotated * BOHR_TO_ANGSTROM):
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

@pytest.mark.filterwarnings("ignore::DeprecationWarning:spglib")
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

@pytest.mark.filterwarnings("ignore::DeprecationWarning:spglib")
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
@pytest.mark.filterwarnings("ignore::DeprecationWarning:spglib")
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
