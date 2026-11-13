"""
kpath_optimizer.py
==================
``KPathOptimizer`` — optimise a k-path by finding the shortest combination
of symmetry-equivalent k-points.

For each k-point in the input path, the class considers all symmetry-
equivalent k-points (the star of k, from the crystal point-group
operations) and searches for the combination that minimises the total
Cartesian path length and individual segment lengths.

Conventions: k-points are in reduced reciprocal coordinates; Cartesian
coordinates are ``k_red @ rlat`` (rows of *rlat* are the reciprocal
vectors, without ``2π``).  Symmetry operations ``W`` are yambopy's
``sym_red`` (acting on reduced direct coordinates as ``x' = x @ W``), so
they act on reduced k-points as ``k' = W^-1 k``.
"""

from __future__ import annotations
import itertools
from typing import Sequence

import numpy as np

def apply_W_to_k(W: np.ndarray, k: np.ndarray) -> np.ndarray:
    """Apply a ``sym_red`` operation to a reduced k-point: ``W^-1 k``."""
    return np.linalg.solve(W, k)

def k_equiv(k1: np.ndarray, k2: np.ndarray, atol: float = 1e-5) -> bool:
    """True if two reduced k-points are equal modulo a reciprocal lattice vector."""
    d = k1 - k2
    return np.allclose(d - np.round(d), 0.0, atol=atol)

def k_in_list(k: np.ndarray, lst: list[np.ndarray], atol: float = 1e-5) -> bool:
    """True if *k* is equivalent (see :func:`k_equiv`) to any k-point in *lst*."""
    return any(k_equiv(k, kr, atol) for kr in lst)

def frac_to_cart_k(k_frac: Sequence, B: np.ndarray) -> np.ndarray:
    """Convert k from reduced to Cartesian coordinates (``k @ B``, no 2π)."""
    return B.T @ np.asarray(k_frac, dtype=float)

def get_star_of_k(k_frac: Sequence, crystal_pg_ops: list[np.ndarray], atol: float = 1e-5) -> np.ndarray:
    """
    Return all symmetry-equivalent k-points (the star of k).

    Parameters
    ----------
    k_frac : array-like, shape (3,)
        k-point in fractional reciprocal coordinates.
    crystal_pg_ops : list of ndarray, shape (3, 3)
        Point-group operations in yambopy's ``sym_red`` convention.
    atol : float
        Tolerance to consider two k-points equal modulo a G-vector.

    Returns
    -------
    star : ndarray, shape (N, 3)
        Distinct k-points in the star.  The input k-point is always
        the first row.
    """
    k0 = np.asarray(k_frac, dtype=float)
    star: list[np.ndarray] = [k0]
    for W in crystal_pg_ops:
        k_new = apply_W_to_k(W, k0)
        if not k_in_list(k_new, star, atol):
            star.append(k_new)
    return np.array(star)

def kpath_segment_length(k1_frac: Sequence, k2_frac: Sequence, B: np.ndarray) -> float:
    """Cartesian distance between two k-points in fractional reciprocal coords."""
    k1_cart = frac_to_cart_k(np.asarray(k1_frac, dtype=float), B)
    k2_cart = frac_to_cart_k(np.asarray(k2_frac, dtype=float), B)
    return float(np.linalg.norm(k2_cart - k1_cart))

class KPathOptimizer:
    """
    Optimise a k-path using symmetry-equivalent k-points.

    Parameters
    ----------
    rlat : ndarray, shape (3, 3)
        Reciprocal lattice vectors as rows (no ``2π``).
    sym : list of ndarray, shape (3, 3)
        Point-group operations in yambopy's ``sym_red`` convention.  Null
        matrices (yambo pads its symmetry arrays with them) are ignored.

    Methods
    -------
    optimize_kpath(kpath_frac, method='greedy_global', ...)
        Return the shortest combination of equivalent k-points.
    """

    def __init__(self, rlat: np.ndarray, sym: list[np.ndarray]):
        self._rlat = rlat
        self.sym = [s for s in sym if not np.all(s == np.zeros((3, 3)))]

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def optimize_kpath(
        self,
        kpath_frac: Sequence[Sequence],
        method: str = "greedy_global",
        max_combinations: int = 100000,
        parallel_segments: bool = True
    ) -> dict:
        """
        Optimize a k-path by finding the shortest combination of
        symmetry-equivalent k-points.

        Parameters
        ----------
        kpath_frac : sequence of array-like, shape (N, 3)
            Input k-path in fractional reciprocal coordinates.
        method : str
            ``'brute'`` (exhaustive, exact), ``'greedy'`` (fast sequential),
            or ``'greedy_global'`` (greedy with multiple starting points).
        max_combinations : int
            Fall back to greedy when brute-force combinations exceed this
            threshold (default 100000).
        parallel_segments : bool
            If True, break ties between (nearly) equally short choices in
            favour of smoother paths: ``'brute'`` adds a penalty for unequal
            segment lengths (up to 10% of their spread over their mean);
            the greedy methods add up to 10% to a segment's length for
            turning back on the previous one (none for carrying straight on).

        Returns
        -------
        dict with keys:
            kpath_cart, kpath_frac, segment_lengths, total_length,
            original_length, improvement (%), n_combinations, stars,
            best_indices (index of the chosen member of each point's star;
            0 is the input point), and ``warning`` when brute force fell
            back to greedy.  ``kpath_frac`` holds the chosen coordinates:
            each star member is used at its copy (shifted by a reciprocal
            lattice vector) nearest to the previous point, so the path stays
            continuous across the zone boundary.
        """
        kpath = [np.asarray(k, dtype=float) for k in kpath_frac]

        if len(kpath) < 2:
            # Nothing to optimise: a single point (or none) keeps its
            # original position, i.e. the first member of its star.
            kpath_frac = np.array(kpath).reshape(-1, 3)
            return {
                "kpath_cart": kpath_frac @ self._rlat,
                "kpath_frac": kpath_frac,
                "segment_lengths": np.array([]),
                "total_length": 0.0,
                "original_length": 0.0,
                "improvement": 0.0,
                "n_combinations": 1,
                "stars": [get_star_of_k(k, self.sym) for k in kpath],
                "best_indices": [0] * len(kpath),
            }

        stars = [get_star_of_k(k, self.sym) for k in kpath]
        star_sizes = [len(s) for s in stars]
        total_combinations = np.prod(star_sizes)

        original_length = self._path_length(kpath)

        if method == "brute":
            if total_combinations > max_combinations:
                result = self._optimize_kpath_greedy(kpath, stars, parallel_segments)
                result["warning"] = (
                    f"Too many combinations ({total_combinations} > "
                    f"{max_combinations}). Used greedy method instead."
                )
            else:
                result = self._optimize_kpath_brute(kpath, stars, parallel_segments)
        elif method == "greedy":
            result = self._optimize_kpath_greedy(kpath, stars, parallel_segments)
        elif method == "greedy_global":
            result = self._optimize_kpath_greedy_global(kpath, stars, parallel_segments)
        else:
            raise ValueError(
                f"Unknown method '{method}'. Use 'brute', 'greedy', or "
                "'greedy_global'."
            )

        result["original_length"] = original_length
        result["improvement"] = (
            (1 - result["total_length"] / original_length) * 100
            if original_length > 0 else 0.0
        )
        return result

    # ------------------------------------------------------------------
    # Internal: choosing copies and comparing lengths
    # ------------------------------------------------------------------

    #: Shifts tried around the rounded difference, (0, 0, 0) first.
    _SHIFTS = np.array(sorted(itertools.product((-1, 0, 1), repeat=3),
                              key=lambda v: (sum(abs(x) for x in v), v)), dtype=float)

    def _nearest_image(self, k: np.ndarray, ref: np.ndarray) -> np.ndarray:
        """
        Copy of the reduced k-point *k* (shifted by a reciprocal lattice
        vector) nearest, in Cartesian distance, to the reduced point *ref*.
        Among equally near copies, the least shifted one (*k* itself if it
        is one of them).
        """
        shifts = np.vstack([np.zeros(3), np.round(ref - k) + self._SHIFTS])
        candidates = k + shifts
        distances = np.linalg.norm((candidates - ref) @ self._rlat, axis=1)
        nearest = np.flatnonzero(~np.array([self._shorter(d, distances.min()) or
                                            self._shorter(distances.min(), d) for d in distances]))
        return candidates[nearest[np.argmin(np.abs(shifts[nearest]).sum(axis=1))]]

    @staticmethod
    def _is_shifted(k: np.ndarray, star_member: np.ndarray) -> bool:
        """True if *k* is a shifted copy of *star_member*, not the member itself."""
        return not np.allclose(k, star_member, atol=1e-8)

    def _chain(self, points: list[np.ndarray]) -> list[np.ndarray]:
        """
        Keep the first point and replace each next one by its copy nearest
        to the previous: since each shift only matters relative to the
        previous point, this minimises every segment, hence the length.
        """
        chained = [np.asarray(points[0], dtype=float)]
        for k in points[1:]:
            chained.append(self._nearest_image(np.asarray(k, dtype=float), chained[-1]))
        return chained

    @staticmethod
    def _shorter(a: float, b: float) -> bool:
        """
        True if *a* is shorter than *b* beyond rounding noise, so that
        (nearly) equal choices keep the first one found.
        """
        if not np.isfinite(b):
            return a < b
        return a < b - 1e-9 * abs(b) - 1e-15

    def _result(self, chosen: list[np.ndarray], indices: list[int],
                stars: list[np.ndarray], n_combinations: int) -> dict:
        """Result dict for the chosen reduced coordinates."""
        chosen = np.array(chosen, dtype=float)
        segment_lengths = [
            kpath_segment_length(chosen[i], chosen[i + 1], self._rlat)
            for i in range(len(chosen) - 1)
        ]
        return {
            "kpath_frac": chosen,
            "kpath_cart": chosen @ self._rlat,
            "segment_lengths": np.array(segment_lengths),
            "total_length": sum(segment_lengths),
            "n_combinations": n_combinations,
            "stars": stars,
            "best_indices": list(indices),
        }

    # ------------------------------------------------------------------
    # Internal: brute force search
    # ------------------------------------------------------------------

    def _optimize_kpath_brute(
        self,
        kpath: list[np.ndarray],
        stars: list[np.ndarray],
        parallel_segments: bool = True,
    ) -> dict:
        """
        Exhaustive search over all combinations of star members, each
        placed at its copy nearest to the previous point (see
        :meth:`_chain`), scoring path length plus, with
        *parallel_segments*, the segment-length uniformity penalty.  Ties
        keep the combination with the fewest shifted copies.
        """
        best_score = float("inf")
        best_indices, best_path = None, None
        n_combinations = 0

        best_shifts = None
        for indices in itertools.product(*[range(len(s)) for s in stars]):
            n_combinations += 1
            members = [stars[i][idx] for i, idx in enumerate(indices)]
            current_path = self._chain(members)
            current_score = self._path_length(current_path)
            if parallel_segments:
                current_score += self._segment_uniformity_penalty(current_path)
            n_shifts = sum(self._is_shifted(k, m) for k, m in zip(current_path, members))
            tie = best_path is not None and not self._shorter(current_score, best_score) \
                and not self._shorter(best_score, current_score)
            if (best_path is None or self._shorter(current_score, best_score)
                    or (tie and n_shifts < best_shifts)):
                best_score, best_shifts = current_score, n_shifts
                best_indices, best_path = list(indices), current_path

        return self._result(best_path, best_indices, stars, n_combinations)

    # ------------------------------------------------------------------
    # Internal: greedy sequential search
    # ------------------------------------------------------------------

    def _extend_greedily(
        self,
        chosen: list[np.ndarray],
        indices: list[int],
        stars: list[np.ndarray],
        parallel_segments: bool,
    ):
        """
        Extend a path (lists *chosen* / *indices*, modified in place) point
        by point: each next point is the star member whose nearest copy is
        closest to the previous choice (plus, with *parallel_segments*, the
        turn penalty); ties prefer an unshifted member, then the lowest index.
        """
        for i in range(len(chosen), len(stars)):
            prev = chosen[-1]
            prev_cart = prev @ self._rlat
            best_dist, best_idx, best_k, best_shifted = float("inf"), 0, None, True
            for j, k_star in enumerate(stars[i]):
                k = self._nearest_image(np.asarray(k_star, dtype=float), prev)
                k_cart = k @ self._rlat
                dist = float(np.linalg.norm(k_cart - prev_cart))
                if parallel_segments and len(chosen) > 1:
                    dist += self._turn_penalty(chosen[-2] @ self._rlat, prev_cart, k_cart, dist)
                shifted = self._is_shifted(k, k_star)
                tie = best_k is not None and not self._shorter(dist, best_dist) \
                    and not self._shorter(best_dist, dist)
                if (best_k is None or self._shorter(dist, best_dist)
                        or (tie and best_shifted and not shifted)):
                    best_dist, best_idx, best_k, best_shifted = dist, j, k, shifted
            chosen.append(best_k)
            indices.append(best_idx)

    def _optimize_kpath_greedy(
        self,
        kpath: list[np.ndarray],
        stars: list[np.ndarray],
        parallel_segments: bool = True,
    ) -> dict:
        """
        Greedy sequential optimisation: keep the first point and pick each
        next one as in :meth:`_extend_greedily`.
        """
        chosen, indices = [np.asarray(stars[0][0], dtype=float)], [0]
        self._extend_greedily(chosen, indices, stars, parallel_segments)
        return self._result(chosen, indices, stars, sum(len(s) for s in stars))

    # ------------------------------------------------------------------
    # Internal: greedy global search
    # ------------------------------------------------------------------

    def _optimize_kpath_greedy_global(
        self,
        kpath: list[np.ndarray],
        stars: list[np.ndarray],
        parallel_segments: bool = True,
    ) -> dict:
        """
        Greedy global: start from every member of the first point's star,
        extend greedily as in :meth:`_extend_greedily`, and keep the
        shortest resulting path (ties keep the first start).
        """
        best = None
        for start_idx in range(len(stars[0])):
            chosen, indices = [np.asarray(stars[0][start_idx], dtype=float)], [start_idx]
            self._extend_greedily(chosen, indices, stars, parallel_segments)
            length = self._path_length(chosen)
            if best is None or self._shorter(length, best[0]):
                best = (length, chosen, indices)

        n_combinations = len(stars[0]) * sum(len(s) for s in stars[1:])
        return self._result(best[1], best[2], stars, n_combinations)

    # ------------------------------------------------------------------
    # Internal: path and segment helpers
    # ------------------------------------------------------------------

    def _path_length(self, kpath: list[np.ndarray]) -> float:
        """Total Cartesian path length."""
        return sum(
            kpath_segment_length(kpath[i], kpath[i + 1], self._rlat)
            for i in range(len(kpath) - 1)
        )

    def _segment_uniformity_penalty(self, kpath: list[np.ndarray]) -> float:
        """Penalty for non-uniform segment lengths (avoids very short + very
        long segments)."""
        if len(kpath) < 3:
            return 0.0
        lengths = [
            kpath_segment_length(kpath[i], kpath[i + 1], self._rlat)
            for i in range(len(kpath) - 1)
        ]
        if not lengths:
            return 0.0
        mean_length = float(np.mean(lengths))
        if mean_length == 0.0:
            return 0.0
        return float(np.std(lengths) / mean_length * 0.1)

    @classmethod
    def _turn_penalty(cls, k0: np.ndarray, k1: np.ndarray, k2: np.ndarray, dist: float) -> float:
        """
        Penalty for the turn at *k1* on the way k0 -> k1 -> k2 (Cartesian):
        0 when the path carries straight on, up to 10% of *dist* when it
        turns straight back.  Used by the greedy methods to choose between
        star members giving (nearly) equal segment lengths.
        """
        angle = cls._angle_between_vectors(k1 - k0, k2 - k1)
        return dist * (angle / np.pi) * 0.1

    @staticmethod
    def _angle_between_vectors(v1: np.ndarray, v2: np.ndarray) -> float:
        """Angle between two vectors in radians."""
        n1 = float(np.linalg.norm(v1))
        n2 = float(np.linalg.norm(v2))
        if n1 == 0.0 or n2 == 0.0:
            return 0.0
        cos_a = float(np.dot(v1, v2) / (n1 * n2))
        return float(np.arccos(np.clip(cos_a, -1.0, 1.0)))
