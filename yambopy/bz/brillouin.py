"""
brillouin.py
============
``BrillouinZone`` — high-symmetry k-paths in the Brillouin zones of the
QE Bravais lattices (see :mod:`yambopy.bz.lattice_data`):

* build the path from ``ibrav`` + parameters or from any cell describing a
  supported lattice, optionally optimised over symmetry-equivalent points;
* interpolate it (``npoints``, ``density`` or ``intervals``) and provide the
  band-structure x-axis, tick positions, labels and indices;
* find the k-points of a grid (e.g. a yambo database) lying on the path,
  whatever the units or cell of the grid's lattice.
"""

import re
from itertools import product
from typing import Any, Sequence
from copy import deepcopy

import numpy as np

from yambopy.bz.lattice_data import find_basis_change, get_lattice_data, lattice_data_for_cell
from yambopy.bz.kpath_optimizer import KPathOptimizer, get_star_of_k
from yambopy.dbs.latticedb import YamboLatticeDB


# =============================================================================
# Path-label regex (compiled once)
# =============================================================================

#: Matches path labels like G, M, K, A1, B1, P2, X1, etc.
_LABEL_RE = re.compile(r"[A-Z][a-z0-9]*")


def _to_builtin(obj: Any) -> Any:
    """Recursively convert numpy arrays and scalars to plain Python types."""
    if isinstance(obj, dict):
        return {key: _to_builtin(value) for key, value in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_to_builtin(value) for value in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.generic):
        return obj.item()
    return obj


# =============================================================================
# BrillouinZone class
# =============================================================================

class BrillouinZone:
    """
    Construct interpolated k-paths in any of the supported Brillouin zones.

    Follows Quantum ESPRESSO's Bravais-lattice classification (ibrav [1]_)
    and the Setyawan-Curtarolo special points and default paths.  Not all
    ibrav values are supported (no monoclinic or triclinic lattices).

    The lattice is given either by *ibrav* + *parameters* or by *cell*.
    Without either, only *path_string* + *intervals* are used: labels,
    ticks and indices are available, but no k-point coordinates.

    Parameters
    ----------
    ibrav : int or None
        QE Bravais-lattice index (1-11, plus -3, -5, -9).
    parameters : dict or None
        Lattice parameters for *ibrav* (e.g. ``{'a': 3.615}``); see
        :func:`~yambopy.bz.lattice_data.get_lattice_data`.
    cell : array-like, shape (3, 3), or None
        Direct lattice vectors as rows, in any orientation, choice of
        primitive vectors and length unit.  The Bravais lattice is
        identified and the special points are given in reduced
        coordinates of this cell (see
        :func:`~yambopy.bz.lattice_data.lattice_data_for_cell`).
    path_string : str or None
        Piecewise special k-point path (comma-separated sections, e.g.
        ``'GMKGALHA,LM,KH'``).  When *None*, uses the default path for the
        lattice type.
    extra_points : dict or None
        Additional or overriding special points, in reduced coordinates of
        the cell (e.g. ``{'M': [0, 0.5, 0], 'K': [1/3, 1/3, 0]}``).
    npoints : int or None
        Total number of k-points along the whole path.  Incompatible with
        *density* and *intervals*.
    density : float or None
        Linear density of k-points, per unit of reciprocal length of the
        cell (e.g. 1/Å or 1/bohr, without 2π).  Incompatible with
        *npoints* and *intervals*.  Default: 5 when none of the three is
        given.
    intervals : list[int] or None
        Number of k-points on each path segment, its end point excluded
        (see :meth:`interpolate_intervals`).  Incompatible with *npoints*
        and *density*.
    optimise : bool
        Replace each special point by the symmetry-equivalent copy giving
        the shortest path (see :meth:`optimise_path`).  Requires a lattice
        and *sym_red*.
    sym_red : list of ndarray, shape (3, 3), or None
        Point-group operations in yambopy's ``sym_red`` convention, used
        when *optimise* is True.

    Attributes
    ----------
    cell, rcell : ndarray, shape (3, 3)
        Direct and reciprocal lattice vectors as rows (``rcell`` without
        2π, so Cartesian k-points are ``k_red @ rcell``).
    variant : str
        Setyawan-Curtarolo lattice variant (e.g. ``'HEX'``, ``'ORCF1'``).
    high_symmetry_points : dict[str, list]
        For each label, the star of the special point (a single entry
        unless the path was optimised); the first member is the original.
    path_sections : list[list[str]]
        Labels of each comma-separated section of the path.
    star_indices : list[list[int]]
        Member of the star used for each point of each section.
    path_string, default_path : str
        Normalised path, and the lattice's default path.
    intervals : list[int]
        Number of k-points on each segment, its end point excluded.
    kpts_red, kpts_car : list[ndarray]
        Interpolated k-points of each section (see :meth:`kpoints`).

    References
    ----------
    .. [1] https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm226
    """

    def __init__(
        self,
        ibrav: int = None,
        parameters: dict[str, float] = None,
        cell: Sequence = None,
        path_string: str = None,
        extra_points: dict[str, Sequence[float]] = None,
        npoints: int = None,
        density: float = None,
        intervals: list[int] = None,
        optimise: bool = False,
        sym_red: list[np.ndarray] = None
    ):
        self.cell = None
        self.high_symmetry_points = {}

        has_lattice = (ibrav is not None and parameters is not None) or cell is not None
        if optimise and (not has_lattice or sym_red is None):
            raise ValueError(
                "To optimise the path, supply a cell or ibrav+parameters, "
                "and the symmetry matrices (sym_red)."
            )
        if sym_red is not None:
            sym_red = [np.asarray(s) for s in sym_red]

        if ibrav is not None and parameters is not None:
            self.cell, self.variant, self.high_symmetry_points, self.default_path = (
                get_lattice_data(ibrav, parameters)
            )
        elif cell is not None:
            self.cell, self.variant, self.high_symmetry_points, self.default_path = (
                lattice_data_for_cell(cell)
            )
        elif path_string is None and intervals is None:
            raise ValueError("Either specify:\n 1: ibrav and parameters\n 2: cell\n 3: path_string and intervals")

        # Reciprocal cell: rows are reciprocal primitive vectors (no 2π),
        # so k_car = k_red @ rcell.
        self.rcell = None
        if self.cell is not None:
            self.rcell = np.linalg.inv(self.cell.T)

        # Merge user-supplied extra points (may overwrite defaults).
        if extra_points is not None and isinstance(extra_points, dict):
            for label, coords in extra_points.items():
                self.high_symmetry_points[label] = np.array(coords)

        for label, coords in self.high_symmetry_points.items():
            self.high_symmetry_points[label] = [coords]

        # Use default path when none is provided.
        if path_string is None:
            path_string = self.default_path

        # Parse path sections.
        self.path_sections: list[list[str]] = []
        for section in path_string.split(","):
            labels = _LABEL_RE.findall(section)
            if labels:
                self.path_sections.append(labels)

        # Validate all labels are recognised.
        if self.high_symmetry_points:
            all_labels_recognized = all(
                label in self.high_symmetry_points
                for section in self.path_sections
                for label in section
            )
            if not all_labels_recognized:
                raise ValueError("Given path contains unrecognised labels.")

        self.star_indices = [[0]*len(section) for section in self.path_sections]

        if optimise and sym_red is not None:
            self.optimise_path(sym_red)

        # Reconstruct normalised path string.
        section_strings = ["".join(s) for s in self.path_sections]
        self.path_string = ",".join(section_strings)

        # Keep intervals
        if intervals is not None:
            self.intervals = intervals

        # Interpolate k-points along the path.
        if self.rcell is not None:
            self.kpts_red: list[np.ndarray] = []
            self.kpts_car: list[np.ndarray] = []
            if intervals is not None:
                self.interpolate_intervals(intervals)
            else:
                self.interpolate(npoints, density)

        # Persist arguments for serialisation.
        self.arguments: dict[str, Any] = {
            "ibrav": ibrav,
            "parameters": parameters,
            "cell": cell,
            "path_string": path_string,
            "extra_points": extra_points,
            "npoints": npoints,
            "density": density,
            "intervals": intervals,
            "optimise": optimise,
            "sym_red": sym_red
        }

    # ------------------------------------------------------------------
    # Serialisation
    # ------------------------------------------------------------------

    def as_dict(self) -> dict[str, Any]:
        """
        Return all arguments needed to reconstruct this instance, as plain
        Python types (numpy arrays become lists) so it can be stored as JSON.
        """
        return _to_builtin(self.arguments)

    @classmethod
    def from_dict(cls, args: dict[str, Any]) -> "BrillouinZone":
        """
        Construct a new instance from a dictionary returned by
        :meth:`as_dict`.
        """
        return cls(
            ibrav=args["ibrav"],
            parameters=args["parameters"],
            cell=args["cell"],
            path_string=args["path_string"],
            extra_points=args["extra_points"],
            npoints=args["npoints"],
            density=args["density"],
            intervals=args["intervals"],
            optimise=args["optimise"],
            sym_red=args["sym_red"]
        )

    # ------------------------------------------------------------------
    # Optimise k-path
    # ------------------------------------------------------------------

    def optimise_path(self,
        sym_red: list[np.ndarray],
        rcell: np.ndarray = None,
        method: str = "greedy_global",
        max_combinations: int = 100000,
        parallel_segments: bool = True):
        """
        Replace each special point of the path by the member of its star
        (symmetry-equivalent copies) that gives the shortest path, using
        :class:`~yambopy.bz.kpath_optimizer.KPathOptimizer`.

        Each comma-separated section is optimised independently.  The
        choices are stored in :attr:`star_indices` and each special point
        in :attr:`high_symmetry_points` becomes its star (first member: the
        original point).  Only the special points are optimised; call
        :meth:`interpolate` (or :meth:`interpolate_intervals`) afterwards
        to rebuild the k-points (the constructor does this when
        ``optimise=True``).

        Parameters
        ----------
        sym_red : list of ndarray, shape (3, 3)
            Point-group operations in yambopy's ``sym_red`` convention,
            referred to the same cell as this ``BrillouinZone`` (e.g. a
            yambo database's ``sym_red`` when the zone was built from its
            ``lat``).
        rcell : ndarray or None
            Reciprocal cell for path lengths; defaults to :attr:`rcell`.
        method, max_combinations, parallel_segments
            Passed to :meth:`~yambopy.bz.kpath_optimizer.KPathOptimizer.optimize_kpath`.
        """

        if rcell is None:
            rcell = self.rcell

        opt = KPathOptimizer(rcell, sym_red)
        for i, section in enumerate(self.path_sections):
            kpath_frac = [self.high_symmetry_points[label][0] for label in section]
            result = opt.optimize_kpath(kpath_frac, method, max_combinations, parallel_segments)
            self.star_indices[i] = result['best_indices']

        for section in self.path_sections:
            for label in section:
                self.high_symmetry_points[label] = get_star_of_k(self.high_symmetry_points[label][0], opt.sym)


    def print_optimisation_info(self):
        """
        Print, for each special point of the path, its coordinates and,
        when the optimisation chose another member of its star, the
        original and chosen coordinates (reduced) and the star index.
        """
        print('Path:')
        for i, section in enumerate(self.path_sections):
            print('---')
            for j, label in enumerate(section):
                rcoords = np.array(self.high_symmetry_points[label][0])
                index = self.star_indices[i][j]
                rcoords_idx = np.array(self.high_symmetry_points[label][index])
                print(f"Point name: {label}")
                if index != 0:
                    print(f"Original coordinates (red): {np.round(rcoords, 4)}")
                    print(f"Optimized coordinates (red): {np.round(rcoords_idx, 4)}")
                    print(f"Corresponding star: {index}")
                else:
                    print(f"Coordinates (red): {np.round(rcoords, 4)}")


    # ------------------------------------------------------------------
    # Interpolation
    # ------------------------------------------------------------------

    def interpolate(self,
        npoints: int = None,
        density: float = None
    ):
        """
        Interpolate evenly-spaced k-points along the path.

        Points are shared between segments in proportion to their lengths,
        with at least 2 per segment; :attr:`intervals` is set accordingly.

        Parameters
        ----------
        npoints : int or None
            Exact total number of k-points along the whole path (unless it
            is too small for 2 points per segment).
        density : float or None
            Linear density of k-points, per unit of reciprocal length of
            the cell (no 2π).  When both *npoints* and *density* are
            *None*, defaults to *density=5*.
        """
        if npoints is not None and density is not None:
            raise ValueError("You may define npoints or density, but not both.")

        # Cartesian length of every segment, in path order.
        seg_lengths = [
            float(np.linalg.norm(section[i + 1] - section[i]))
            for section in self.special_kpoints("car")
            for i in range(len(section) - 1)
        ]
        length = sum(seg_lengths)
        nsegments = len(seg_lengths)

        if npoints is None:
            if density is None:
                density = 5.0
            npoints = int(round(length * density))

        # Each segment contributes `n` points (its end point excluded) and
        # each section adds its final point, so the segments share
        # npoints - nsections points, distributed in proportion to their
        # lengths.  Every non-degenerate segment gets at least 2 points.
        budget = max(npoints - len(self.path_sections), 2 * nsegments)
        intervals: list[int] = []
        remaining_length = length
        remaining_points = budget
        for seg_len in seg_lengths:
            if seg_len < 1e-12 or remaining_length < 1e-12:
                n = 1
            else:
                n = max(2, int(round(seg_len * remaining_points / remaining_length)))
            intervals.append(n)
            remaining_length -= seg_len
            remaining_points = max(remaining_points - n, 0)

        self._build_kpoints(intervals)
        self.intervals = intervals

    def interpolate_intervals(self, intervals: list[int]):
        """
        Interpolate with a specific number of points on each segment.

        Parameters
        ----------
        intervals : list[int]
            Number of points on each path segment (in order along the
            path), the segment's start point included and its end point
            excluded; each section is closed with its last special point.
            Its length must equal the number of segments in the path.
        """
        num_intervals = sum(len(s) - 1 for s in self.path_sections)
        if num_intervals != len(intervals):
            raise ValueError(
                f"Number of intervals ({len(intervals)}) does not match the "
                f"number of path segments ({num_intervals})."
            )

        self._build_kpoints(intervals)
        self.intervals = list(intervals)

    def _build_kpoints(self, intervals: list[int]):
        """
        Fill :attr:`kpts_red` and :attr:`kpts_car` (one array per path
        section) placing ``intervals[i]`` evenly spaced points on segment
        *i* (its end point excluded) and closing each section with its
        last special point.
        """
        self.kpts_red = []
        self.kpts_car = []
        i = 0

        for spoints in self.special_kpoints("red", merge_sections=False):
            kpoints = np.array(spoints, dtype=float)
            section_kpts: list[np.ndarray] = []

            for kpt, diff in zip(kpoints[:-1], kpoints[1:] - kpoints[:-1]):
                for t in np.linspace(0, 1, intervals[i], endpoint=False):
                    section_kpts.append(kpt + t * diff)
                i += 1

            if len(kpoints) > 0:
                section_kpts.append(kpoints[-1])

            red = np.array(section_kpts) if section_kpts else np.empty((0, 3))
            self.kpts_red.append(red)
            self.kpts_car.append(red @ self.rcell)

    # ------------------------------------------------------------------
    # K-point access
    # ------------------------------------------------------------------

    def kpoints(self, coords: str = "red", qe: bool = False) -> np.ndarray:
        """
        Return the coordinates of all interpolated k-points on the path.

        Parameters
        ----------
        coords : ``'red'`` or ``'car'``
            ``'red'`` for reduced reciprocal coordinates, ``'car'`` for
            Cartesian (``k_red @ self.rcell``, no 2π).
        qe : bool
            If *True*, append a trailing column of ones (QE format).

        Returns
        -------
        ndarray, shape (N, 3) or (N, 4)
        """
        if coords == "red":
            kpts = np.concatenate(self.kpts_red) if self.kpts_red else np.empty((0, 3))
        elif coords == "car":
            kpts = np.concatenate(self.kpts_car) if self.kpts_car else np.empty((0, 3))
        else:
            raise ValueError(f"coords: {coords} not supported. Use 'red' or 'car'.")

        if qe:
            return np.pad(kpts, [(0, 0), (0, 1)], constant_values=1)
        return kpts

    # ------------------------------------------------------------------
    # Special k-points
    # ------------------------------------------------------------------

    def special_kpoints(self, coords: str = "red", merge_sections: bool = False):
        """
        Return the coordinates of special (high-symmetry) k-points defining
        the path.

        Parameters
        ----------
        coords : ``'red'`` or ``'car'``
            Reduced, or Cartesian (``k_red @ self.rcell``, no 2π).
        merge_sections : bool
            If *True*, return a single flat array joining all path sections.

        Returns
        -------
        list[list[ndarray]] (one list per section) or ndarray, shape (N, 3)
        """
        if coords == "red":
            sub_lists = [
                [np.array(self.high_symmetry_points[label][self.star_indices[i][j]]) for j, label in enumerate(section)]
                for i, section in enumerate(self.path_sections)
            ]
        elif coords == "car":
            sub_lists = [
                [self.high_symmetry_points[label][self.star_indices[i][j]] @ self.rcell for j, label in enumerate(section)]
                for i, section in enumerate(self.path_sections)
            ]
        else:
            raise ValueError(f"coords: {coords} not supported. Use 'red' or 'car'.")

        if merge_sections:
            return np.concatenate(
                [np.array(sl) for sl in sub_lists]
            ) if sub_lists else np.empty((0, 3))
        return sub_lists

    # ------------------------------------------------------------------
    # Distances
    # ------------------------------------------------------------------

    def kpoints_distances(self) -> np.ndarray:
        """
        Cumulative Cartesian distances along the interpolated path.

        Returns one distance per k-point of :meth:`kpoints`, so it can be
        used directly as the band-structure x-axis.  The jump between the
        end of a section and the start of the next one is not counted:
        both points share the same x, where :meth:`special_kpoints_distances`
        (with ``merge_sections=True``) places the merged ``'A — B'`` tick.
        """
        kpt_dists: list[float] = []
        dist = 0.0
        for kpts in self.kpts_car:
            for nk in range(len(kpts)):
                if nk > 0:
                    dist += float(np.linalg.norm(kpts[nk] - kpts[nk - 1]))
                kpt_dists.append(dist)
        return np.array(kpt_dists)

    def special_kpoints_distances(self, merge_sections: bool = False) -> list[np.ndarray] | np.ndarray:
        """
        Cumulative Cartesian distances of the special k-points along the
        path.

        Parameters
        ----------
        merge_sections : bool
            If *False* (default), return one array per section, each
            starting at 0.  If *True*, return a single array accumulated
            across sections, where the end of a section and the start of
            the next share one entry: these are the band-structure tick
            positions, matching ``path_labels_list(True)``.

        Returns
        -------
        list[ndarray] or ndarray
        """
        spoints_piecewise = self.special_kpoints("car")
        sdistances: list[np.ndarray] = []

        for section in spoints_piecewise:
            distance = 0
            distances = [0]
            for nk in range(len(section) - 1):
                distance += float(np.linalg.norm(section[nk + 1] - section[nk]))
                distances.append(distance)
            sdistances.append(np.array(distances))

        if not merge_sections:
            return sdistances

        mdistances = sdistances[0].tolist()
        for dists in sdistances[1:]:
            accdistances = mdistances[-1] + dists
            for accdist in accdistances[1:]:
                mdistances.append(accdist)

        return np.array(mdistances)

    # ------------------------------------------------------------------
    # Labels
    # ------------------------------------------------------------------

    def path_labels_list(self, merge_sections: bool = False) -> list[str]:
        """
        Return a flat list of special-point labels for band-structure
        x-axis labelling.

        Parameters
        ----------
        merge_sections : bool
            If *True*, the last label of a section and the first of the
            next are merged into one ``'A — B'`` entry, as they share one
            tick on the band-structure axis.

        Returns
        -------
        list[str]
        """
        labels_list = [
            label for section in self.path_sections for label in section
        ]

        if merge_sections:
            section_lengths = [len(s) for s in self.path_sections]
            if len(section_lengths) > 0:
                # Last entry is never a boundary.
                section_lengths.pop()
                boundaries = np.cumsum(section_lengths)
                for k in boundaries:
                    labels_list[k - 1] = (
                        f"{labels_list[k - 1]} — {labels_list[k]}"
                    )
                for k in sorted(boundaries, reverse=True):
                    labels_list.pop(k)

        return labels_list

    # ------------------------------------------------------------------
    # Ticks
    # ------------------------------------------------------------------

    def get_indices(self) -> list[list]:
        """
        Return ``[index, label]`` pairs locating each special point in the
        array returned by :meth:`kpoints`.

        Labels are those of :meth:`path_labels_list` with
        ``merge_sections=True``; a merged ``'A — B'`` label points to *A*,
        the last k-point of its section (*B* is the next one).
        """
        indices = []
        index = 0
        n = 0  # running segment index into self.intervals
        for s, section in enumerate(self.path_sections):
            if s == 0:
                indices.append([index, section[0]])
            for label in section[1:]:
                index += self.intervals[n]
                n += 1
                indices.append([index, label])
            if s < len(self.path_sections) - 1:
                # Merge this section's end with the next section's start.
                indices[-1][1] = f"{indices[-1][1]} — {self.path_sections[s + 1][0]}"
                index += 1
        return indices

    # ------------------------------------------------------------------
    # Collinear k-points
    # ------------------------------------------------------------------

    def find_collinear_grid(
        self,
        latticedb: YamboLatticeDB,
        expand_kpoints: bool = True,
        tol: float = 1e-4
    ):
        """
        Find the k-points of a yambo database grid that lie on the path.

        The grid k-points are taken in reduced coordinates of the database
        lattice and re-expressed in this ``BrillouinZone``'s cell, so the
        two may differ in units (e.g. Å vs bohr), orientation and choice of
        primitive vectors, as long as they describe the same lattice.

        Parameters
        ----------
        latticedb : YamboLatticeDB
        expand_kpoints : bool
            If *True*, search the full-BZ grid (expanding the database copy
            if needed) and return indices into the irreducible k-points.
            If *False*, search the database k-points as they are.
        tol : float
            Tolerance passed to :meth:`find_collinear_kpoints`.

        Returns
        -------
        Same as :meth:`find_collinear_kpoints`; the k-points are Cartesian
        in this ``BrillouinZone``'s convention (``k_red @ self.rcell``).
        """
        latdb = latticedb
        if expand_kpoints:
            if latticedb.ibz_nkpoints == latticedb.nkpoints:
                # Expand a copy: expanding modifies the database in place.
                latdb = deepcopy(latticedb)
                latdb.expand_kpoints(atol=1e-6, verbose=0)
            kmap = latdb.BZ_to_IBZ_indexes

        col_kpts, col_idx, col_dsts = self.find_collinear_lattice_kpoints(
            latdb.lat, latdb.car_kpoints, tol=tol
        )

        if expand_kpoints:
            col_idx = np.asarray(kmap)[col_idx]
        return col_kpts, col_idx, col_dsts

    def find_collinear_lattice_kpoints(
        self,
        lat: np.ndarray,
        car_kpoints: np.ndarray,
        tol: float = 1e-4,
    ):
        """
        Find which k-points, given for another description of the same
        lattice, lie on the path.

        Parameters
        ----------
        lat : ndarray, shape (3, 3)
            Direct lattice vectors (rows) the k-points refer to, e.g. a yambo
            database's ``lat``.  It may differ from ``self.cell`` in units,
            orientation and choice of primitive vectors.
        car_kpoints : ndarray, shape (M, 3)
            Cartesian k-points in the convention of *lat*:
            ``k_red @ rec_lat(lat)``, without ``2π`` (yambopy's
            ``car_kpoints``).
        tol : float
            Tolerance passed to :meth:`find_collinear_kpoints`.

        Returns
        -------
        Same as :meth:`find_collinear_kpoints`, with indices into
        *car_kpoints* and k-points in this ``BrillouinZone``'s convention.
        """
        lat = np.asarray(lat, dtype=float)
        # rec_lat(lat) = inv(lat).T, so k_red = car_kpoints @ lat.T.
        red_kpts = np.asarray(car_kpoints, dtype=float) @ lat.T

        # Basis change between the two cells, compared at equal volume.
        def normalised(c):
            return c / abs(np.linalg.det(c)) ** (1 / 3)

        M = find_basis_change(normalised(self.cell), normalised(lat), rtol=1e-3)
        if M is None:
            raise ValueError(
                "The lattice of the database does not match the lattice of "
                "this BrillouinZone (up to units, rotation and choice of "
                "primitive vectors)."
            )
        # lat ~ M @ cell, so k_red(lat) = M @ k_red(cell).
        red_kpts = red_kpts @ np.linalg.inv(M).T

        return self.find_collinear_kpoints(red_kpts @ self.rcell, tol=tol)

    # ------------------------------------------------------------------
    # Collinear k-points on the optimised path
    # ------------------------------------------------------------------

    def find_collinear_kpoints(
        self,
        kpoints_car: np.ndarray,
        tol: float = 1e-4,
    ):
        """
        Find which k-points from a given set lie on the (possibly
        optimised) path, ordered by distance along that path.

        The input *kpoints_car* must be in the convention used by the rest
        of ``BrillouinZone`` (``k_red @ self.rcell``, *i.e.* **without** the
        ``2π`` factor).  Use :meth:`find_collinear_grid` for a yambo
        database grid.

        Input k-points are folded into ``[-1/2, 1/2)`` (reduced) and their
        G-shifted copies, over a range of shifts covering the whole path,
        are tested for collinearity with each segment.

        Results are **deduplicated**: a k-point (same original index and same
        Cartesian position, modulo a G-vector) is reported only once, at the
        earliest position where it is found along the path.  This avoids
        duplicate matches on adjacent segments that share a boundary.

        Parameters
        ----------
        kpoints_car : ndarray, shape (M, 3)
            Cartesian coordinates of k-points to test, in the ``BrillouinZone``
            convention (``k_frac @ self.rcell``, no ``2π`` factor).
        tol : float
            Maximum distance from the segment line (in the units of
            ``self.rcell``), and maximum overshoot past either segment end
            as a fraction of the segment length.

        Returns
        -------
        collinear_kpoints : ndarray, shape (N, 3)
            Cartesian coordinates (BZ convention) of the found k-points,
            **deduplicated** at segment boundaries.
        collinear_indices : ndarray, shape (N,)
            Index of each found k-point in the input *kpoints_car* array.
        collinear_distances : ndarray, shape (N,)
            Cumulative distance along the path for each found k-point.
        """
        kpoints_car = np.asarray(kpoints_car, dtype=float)
        if kpoints_car.ndim != 2 or kpoints_car.shape[1] != 3:
            raise ValueError(
                "kpoints_car must be a 2-D array of shape (M, 3)."
            )

        rcell = self.rcell

        # Fold the input into [-1/2, 1/2) in reduced coordinates.
        red = np.linalg.solve(rcell.T, kpoints_car.T).T
        red -= np.floor(red + 0.5)
        kpoints_car = red @ rcell

        # G-shifts: enough for the folded points to reach every point of
        # the path, whose (optimised) special points may lie outside the
        # first cell.
        path_red = self.special_kpoints("red", merge_sections=True)
        lo = np.floor(path_red.min(axis=0) - 0.5).astype(int)
        hi = np.ceil(path_red.max(axis=0) + 0.5).astype(int)
        gvecs = np.array(
            list(product(*(range(l, h + 1) for l, h in zip(lo, hi)))), dtype=float
        )
        shifts_cart = gvecs @ rcell

        collinear_kpoints: list[np.ndarray] = []
        collinear_distances: list[float] = []
        collinear_indices: list[int] = []
        prev_segment_keys: set[tuple] = set()
        path_distance = 0.0

        for seg_cart in self.special_kpoints('car'):
            # Sections are separated by commas in the path string; the
            # last point of one section is not adjacent to the first of
            # the next → reset boundary-dedup state.
            prev_segment_keys.clear()

            for i in range(len(seg_cart) - 1):
                A = seg_cart[i]
                B = seg_cart[i + 1]
                AB = B - A
                AB_len = np.linalg.norm(AB)

                if AB_len < 1e-15:
                    continue

                AB_hat = AB / AB_len
                current_segment_keys: set[tuple] = set()

                for idx, P in enumerate(kpoints_car):
                    for shift in shifts_cart:
                        P_shifted = P + shift
                        AP = P_shifted - A
                        proj = float(np.dot(AP, AB_hat))

                        if proj < -tol * AB_len or proj > AB_len * (1.0 + tol):
                            continue

                        cross = float(np.linalg.norm(np.cross(AP, AB)))
                        if cross > tol * AB_len:
                            continue

                        # Dedup only against the immediately preceding segment
                        # (catches boundary duplicates between consecutive segments).
                        key = (idx, tuple(np.round(P_shifted, decimals=8)))
                        if key in prev_segment_keys:
                            continue
                        current_segment_keys.add(key)

                        collinear_kpoints.append(P_shifted)
                        collinear_distances.append(path_distance + proj)
                        collinear_indices.append(idx)
                        break

                prev_segment_keys = current_segment_keys
                path_distance += AB_len

        if collinear_distances:
            order = np.argsort(collinear_distances, kind="stable")
            return (
                np.array(collinear_kpoints)[order],
                np.array(collinear_indices)[order],
                np.array(collinear_distances)[order]
            )

        return (
            np.empty((0, 3)),
            np.empty(0, dtype=int),
            np.empty(0)
        )
