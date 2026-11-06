"""
brillouin_v04.py
================
Extends ``brillouin_v03``: ``find_collinear_kpoints`` now **deduplicates**
k-points that are found on multiple adjacent segments (segment boundaries).

Provides:

* ``get_lattice_data`` — returns lattice vectors, Bravais variant, high-symmetry
  k-point coordinates and a default band-structure path for a given QE ibrav.

* ``BrillouinZone`` — constructs interpolated k-paths in any of the supported
  Brillouin zones.
  New in v04: deduplicated ``find_collinear_kpoints()``.

Dependencies
------------
* numpy
* yambopy (``yambopy.kpoints.expand_kpoints``, ``yambopy.lattice.red_car``,
  ``yambopy.lattice.isbetween``)
"""

# from __future__ import annotations

import re
from itertools import product
from typing import Any, Sequence
from copy import deepcopy

import numpy as np

from yambopy.bz.lattice_data import get_lattice_data
from yambopy.bz.kpath_optimizer import KPathOptimizer, get_star_of_k
from yambopy.bz.cell_to_ibrav import CellToIbrav
from yambopy.dbs.latticedb import YamboLatticeDB
from yambopy.kpoints import expand_kpoints
from yambopy.lattice import isbetween, red_car


# =============================================================================
# Path-label regex (compiled once)
# =============================================================================

#: Matches path labels like G, M, K, A1, B1, P2, X1, etc.
_LABEL_RE = re.compile(r"[A-Z][a-z0-9]*")


# =============================================================================
# BrillouinZone class
# =============================================================================

class BrillouinZone:
    """
    Construct interpolated k-paths in any of the supported Brillouin zones.

    Follows Quantum ESPRESSO's Bravais-lattice classification (ibrav [1]_).
    Not all ibrav values are supported.

    Parameters
    ----------
    ibrav : int
        QE Bravais-lattice index (1-11, plus -3, -5, -9).
    parameters : dict or None
        Lattice parameters (e.g. ``{'a': 3.615}``).  Required keys depend
        on *ibrav*.
    path_string : str or None
        Piecewise special k-point path (comma-separated sections, e.g.
        ``'GMKGALHA,LM,KH'``).  When *None*, uses the default path for the
        lattice type.
    extra_points : dict or None
        Additional or overriding special k-point definitions.
        (e.g. ``{'M': [0, 0.5, 0], 'K': [1/3, 1/3, 0]}``)
    npoints : int or None
        Total number of k-points along the whole path.  Incompatible with
        *density* and *intervals*.
    density : float or None
        Linear density of k-points (1/A).  Incompatible with *npoints* and
        *intervals*.  Default: 5 when none of the three is given.
    intervals : list[int] or None
        Number of k-points on each path segment.  Incompatible with *npoints*
        and *density*.

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

        if ibrav is not None and parameters is not None:
            self.cell, self.variant, self.high_symmetry_points, self.default_path = (
                get_lattice_data(ibrav, parameters)
            )
        elif optimise and (cell is None or sym_red is None):
                raise ValueError("To optimise the path, a cell or ibrav+parameters, and symmetry matrices must be supplied")
        elif cell is not None:
            cell2ibrav = CellToIbrav(cell)
            print(cell2ibrav.summary())

            self.cell, self.variant, self.high_symmetry_points, self.default_path = (
                get_lattice_data(cell2ibrav.ibrav, cell2ibrav.parameters_dict())
            )
        elif path_string is None and intervals is None:
            raise ValueError("Either specify:\n 1: ibrav and parameters\n 2: cell\n 3: path_string and intervals")

        # Reciprocal cell: columns are reciprocal primitive vectors.
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
        """Return all arguments needed to reconstruct this instance."""
        return self.arguments

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
        Piecewise optimised k-path using :class:`~kpath_optimizer.KPathOptimizer`.

        Each comma-separated section of the path is optimised independently.
        Returns a list of result dicts (one per section) as returned by
        :meth:`~kpath_optimizer.KPathOptimizer.optimize_kpath`.

        Notes
        -----
        * The optimisation runs on the *special k-points* defining each section;
          interpolated points are **not** re-optimised.
        * ``KSymmetry`` works in QE primitive fractional coordinates, while
          ``BrillouinZone.high_symmetry_points`` are given in the coordinate
          system of the cell returned by :func:`get_lattice_data`.  These
          coincide for primitive Bravais types (``aP``, ``mP``, ``oP``, ``tP``,
          ``hP``, ``cP``).  For centred lattices the coordinates are already in
          the primitive-cell basis for this implementation, so the integration
          is consistent.
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

        Parameters
        ----------
        npoints : int or None
            Total number of k-points along the whole path.
        density : float or None
            Linear density of k-points (1/A).  When both *npoints* and
            *density* are *None*, defaults to *density=5*.
        """
        if npoints is not None and density is not None:
            raise ValueError("You may define npoints or density, but not both.")

        length = sum(self.special_kpoints_distances(merge_sections=True))

        if npoints is None:
            if density is None:
                density = 5.0
            npoints = int(round(length * density))

        points = self.special_kpoints("red", merge_sections=False)
        lengths = self.special_kpoints_distances(merge_sections=True)

        self.kpts_red = []
        self.kpts_car = []
        accumulated: list[float] = []
        x0 = 0.0
        section_index = 1  # offset into the flat lengths array

        intervals = []

        for spoints in points:
            kpoints = np.array(spoints)
            diffs = kpoints[1:] - kpoints[:-1]
            segment_kpts_red: list[np.ndarray] = []

            for kpt, diff, seg_len in zip(
                kpoints[:-1],
                diffs,
                lengths[section_index: section_index + len(kpoints) - 1],
            ):
                remaining = length - x0
                if abs(remaining) < 1e-6:
                    n = 0
                else:
                    n = max(2, int(round(seg_len * (npoints - len(accumulated)) / remaining)))

                for t in np.linspace(0, 1, n, endpoint=False):
                    segment_kpts_red.append(kpt + t * diff)
                    accumulated.append(x0 + t * seg_len)

                x0 += seg_len

                intervals.append(n)

            section_index += len(kpoints) - 1

            if len(kpoints) > 0:
                segment_kpts_red.append(kpoints[-1])

            if not segment_kpts_red:
                segment_kpts_red = []

            self.kpts_car.append(
                np.array([np.dot(kpt, self.rcell) for kpt in segment_kpts_red])
            )
            self.kpts_red.append(np.array(segment_kpts_red) if segment_kpts_red else np.empty((0, 3)))

        self.intervals = intervals

    def interpolate_intervals(self, intervals: list[int]):
        """
        Interpolate with a specific number of points on each segment.

        Parameters
        ----------
        intervals : list[int]
            Number of points on each path segment (in order along the path).
            The total length must equal the number of segments defined by
            the path string.
        """
        num_intervals = sum(len(s) - 1 for s in self.path_sections)
        if num_intervals != len(intervals):
            raise ValueError(
                f"Number of intervals ({len(intervals)}) does not match the "
                f"number of path segments ({num_intervals})."
            )

        points = self.special_kpoints("red", merge_sections=False)

        self.kpts_red = []
        self.kpts_car = []
        i = 0

        for spoints in points:
            kpoints = np.array(spoints)
            segment_kpts_red: list[np.ndarray] = []

            for kpt, _diff in zip(kpoints, kpoints[1:] - kpoints[:-1]):
                for t in np.linspace(0, 1, intervals[i], endpoint=False):
                    segment_kpts_red.append(kpt + t * _diff)
                i += 1

            if len(kpoints) > 0:
                segment_kpts_red.append(kpoints[-1])

            self.kpts_car.append(
                np.array([np.dot(kpt, self.rcell) for kpt in segment_kpts_red])
            )
            self.kpts_red.append(np.array(segment_kpts_red) if segment_kpts_red else np.empty((0, 3)))

    # ------------------------------------------------------------------
    # K-point access
    # ------------------------------------------------------------------

    def kpoints(self, coords: str = "red", qe: bool = False) -> np.ndarray:
        """
        Return the coordinates of all interpolated k-points on the path.

        Parameters
        ----------
        coords : ``'red'`` or ``'car'``
            ``'red'`` for fractional reciprocal, ``'car'`` for Cartesian.
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
        merge_sections : bool
            If *True*, return a single flat array joining all path sections.

        Returns
        -------
        list[list[ndarray]] or ndarray
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

        Useful for band-structure x-axis plotting.
        """
        kpt_dists = [0.0]
        dist = 0.0
        for kpts in self.kpts_car:
            for nk in range(len(kpts) - 1):
                dist += float(np.linalg.norm(kpts[nk + 1] - kpts[nk]))
                kpt_dists.append(dist)
        return np.array(kpt_dists)

    def special_kpoints_distances(self, merge_sections: bool = False) -> list[np.ndarray] | np.ndarray:
        """
        Cumulative Cartesian distances between consecutive special
        k-points on the path.

        Parameters
        ----------
        merge_sections : bool
            If *False* (default), each section starts again at the
            cumulative distance of the previous sections.  If *True*,
            distances are accumulated monotonically across sections.

        Returns
        -------
        ndarray
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
            If *True*, consecutive sections are joined with a ``' — '``
            separator at the shared boundary.

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

    def get_indices(self) -> list[tuple]:
        indices = []
        index = 0
        labels = self.path_labels_list(True)
        for n, label in enumerate(labels[:-1]):
            indices.append([index, label])
            index += self.intervals[n]
        indices.append([index, labels[-1]])
        return indices

    # ------------------------------------------------------------------
    # Collinear k-points
    # ------------------------------------------------------------------

    def get_collinear_kpoints(self, kpoints_car, sym_car=None, debug=False):
        """
        Given a set of k-points, return those which are collinear with
        the path.

        Parameters
        ----------
        kpoints_car : ndarray
            Cartesian coordinates of k-points (full or irreducible BZ).
        sym_car : list or None
            Symmetry operations to expand IBZ k-points, or *None* if
            *kpoints_car* already covers the full BZ.
        debug : bool
            Print debug information.

        Returns
        -------
        collinear_kpoints : ndarray
            Cartesian coordinates of collinear k-points.
        collinear_indices : ndarray
            Indices of the matching k-points.
        collinear_distances : ndarray
            Distances of the projected k-points along the path.
        """
        rlat = self.rcell

        if sym_car is None:
            kpoints_indices = list(range(len(kpoints_car)))
        else:
            _, kpoints_indices, _, kpoints_car_expanded = expand_kpoints(
                kpoints_car, sym_car, rlat
            )
            kpoints_car = kpoints_car_expanded

        spoints_car = self.special_kpoints("car")
        spoints_distances = self.special_kpoints_distances(True)

        collinear_kpoints = []
        collinear_indices = []
        collinear_distances = []
        prev_segment_keys: set[tuple] = set()
        kdist = 0

        for section in spoints_car:

            prev_segment_keys.clear()

            for k in range(len(section) - 1):
                data = {}
                start_kpoint = section[k]
                end_kpoint = section[k + 1]
                current_segment_keys: set[tuple] = set()

                for x, y, z in product(range(-1, 2), repeat=3):
                    shift = red_car(np.array([[x, y, z]]), rlat)[0]

                    for index, kpoint in zip(kpoints_indices, kpoints_car):
                        kpoint_shift = kpoint + shift

                        if isbetween(start_kpoint, end_kpoint, kpoint_shift):
                            # Dedup only against the immediately preceding segment
                            # (catches boundary duplicates between consecutive segments).
                            dup_key = (index, tuple(np.round(kpoint_shift, decimals=8)))
                            if dup_key in prev_segment_keys:
                                continue
                            current_segment_keys.add(dup_key)

                            key = tuple(
                                round(kpt, 4) for kpt in kpoint_shift
                            )
                            distance = float(
                                np.linalg.norm(start_kpoint - kpoint_shift)
                            )
                            distance_within_path = (
                                spoints_distances[kdist] + distance
                            )
                            data[key] = [
                                index, distance, kpoint_shift, distance_within_path,
                            ]
                prev_segment_keys = current_segment_keys
                kdist += 1

                # Sort by distance from start of segment (nearest first).
                sorted_data = sorted(data.values(), key=lambda i: i[1])

                for index, distance, kpoint_shift, distance_within_path in sorted_data:
                    collinear_indices.append(index)
                    collinear_kpoints.append(kpoint_shift)
                    collinear_distances.append(distance_within_path)

                    if debug:
                        print(
                            f"{kpoint_shift[0]:12.8f} {kpoint_shift[1]:12.8f} "
                            f"{kpoint_shift[2]:12.8f}  {index}  {distance}  "
                            f"{distance_within_path}"
                        )

        return (
            np.array(collinear_kpoints),
            np.array(collinear_indices),
            np.array(collinear_distances),
        )

    def find_collinear_grid(
        self,
        latticedb: YamboLatticeDB,
        expand_kpoints: bool = True,
        tol: float = 1e-4
    ):
        # Safe copy (deleted below)
        latdb = deepcopy(latticedb)
        rlat = latdb.rlat
        # Expand if `lattice` unexpanded
        if expand_kpoints:
            if latdb.ibz_nkpoints == latdb.nkpoints:
                latdb.expand_kpoints(atol=1e-6, verbose=1)
            kmap = latdb.BZ_to_IBZ_indexes
        # red_kpts = latdb.red_kpoints
        car_kpts = latdb.car_kpoints
        # ibz_red_kpts = latdb.get_ibz_kpoints(units='red')
        # ibz_car_kpts = latdb.get_ibz_kpoints(units='car')

        # inv_kmap = latdb.IBZ_to_BZ_indexes
        # sym_red = latdb.sym_red

        del latdb

        col_kpts, col_idx, col_dsts = self.find_collinear_kpoints(car_kpts, rcell=rlat, tol=tol)

        if expand_kpoints:
            col_idx = kmap[col_idx]
        return col_kpts, col_idx, col_dsts

    # ------------------------------------------------------------------
    # Collinear k-points on the optimised path
    # ------------------------------------------------------------------

    def find_collinear_kpoints(
        self,
        kpoints_car: np.ndarray,
        # sym_red: list[np.ndarray],
        rcell: np.ndarray | None = None,
        tol: float = 1e-4,
    ):
        """
        Find which k-points from a given set lie on the piecewise optimised
        path, ordered by distance along that path.

        The optimised path is obtained from :attr:`optimized_path` and its
        Cartesian coordinates are internally converted to the same convention
        used by the rest of ``BrillouinZone`` (``k_frac @ self.rcell``, *i.e.*
        **without** the ``2π`` factor).  The input *kpoints_car* must be in the
        **same** convention.

        For each segment, G-shifted copies of each input k-point (shifts in
        ``[-1, 0, 1]`` per reciprocal-lattice direction) are tested for
        collinearity.

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
            Relative tolerance for collinearity (fraction of segment length).

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

        if rcell is None:
            rcell = self.rcell

        # self.opt_path = self.compute_optimized_path(sym_red, rcell,'brute', parallel_segments=False)
        # print(self.opt_path)

        collinear_kpoints: list[np.ndarray] = []
        collinear_distances: list[float] = []
        collinear_indices: list[int] = []
        prev_segment_keys: set[tuple] = set()
        path_distance = 0.0

        # Pre-compute G-shift Cartesian vectors (BZ convention).
        gvecs = np.array(list(product(range(-1, 2), repeat=3)), dtype=float)
        shifts_cart = gvecs @ rcell  # (27, 3)
        # shifts_cart = red_car(gvecs, rcell)

        for seg_cart in self.special_kpoints('car'):
            # Sections are separated by commas in the path string; the
            # last point of one section is not adjacent to the first of
            # the next → reset boundary-dedup state.
            prev_segment_keys.clear()

            # Convert optimised kpath_cart from KSymmetry convention
            # (which includes 2π) to BZ convention (k_frac @ rcell).
            # kpath_frac = section["kpath_frac"]
            # seg_cart = kpath_frac @ rcell  # (Nk, 3), BZ convention
            # seg_cart = red_car(kpath_frac, rcell)
            # seg_lengths = section["segment_lengths"]

            for i in range(len(seg_cart) - 1):
                A = seg_cart[i]
                B = seg_cart[i + 1]
                AB = B - A
                # AB_len = seg_lengths[i] / (2 * np.pi)
                # AB_len = seg_lengths[i]
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
            order = np.argsort(collinear_distances)
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
