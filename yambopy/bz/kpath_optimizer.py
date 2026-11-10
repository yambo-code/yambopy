"""
kpath_optimizer.py
==================
``KPathOptimizer`` — optimise a k-path by finding the shortest combination
of symmetry-equivalent k-points.

For each k-point in the input path, the class considers all symmetry-
equivalent k-points (the star of k, obtained from a :class:`KSymmetry`
instance) and searches for the combination that minimises the total
Cartesian path length and individual segment lengths.

Dependencies: numpy, ksymmetry.
"""

from __future__ import annotations
import itertools
from typing import Sequence

import numpy as np

def apply_W_to_k(W: np.ndarray, k: np.ndarray) -> np.ndarray:
    return np.linalg.solve(W, k)

def k_equiv(k1: np.ndarray, k2: np.ndarray, atol: float = 1e-5) -> bool:
    d = k1 - k2
    return np.allclose(d - np.round(d), 0.0, atol=atol)

def k_in_list(k: np.ndarray, lst: list[np.ndarray], atol: float = 1e-5) -> bool:
    return any(k_equiv(k, kr, atol) for kr in lst)

def frac_to_cart_k(k_frac: Sequence, B: np.ndarray) -> np.ndarray:
    """Convert k from fractional reciprocal to Cartesian (2*pi/Ang)."""
    return B.T @ np.asarray(k_frac, dtype=float)

def cart_to_frac_k(k_cart: Sequence, B: np.ndarray) -> np.ndarray:
    """Convert k from Cartesian (2*pi/Ang) to fractional reciprocal."""
    return np.linalg.solve(B.T, np.asarray(k_cart, dtype=float))

def get_star_of_k(k_frac: Sequence, crystal_pg_ops: list[np.ndarray], atol: float = 1e-5) -> np.ndarray:
    """
    Return all symmetry-equivalent k-points (the star of k).

    Parameters
    ----------
    k_frac : array-like, shape (3,)
        k-point in fractional reciprocal coordinates.

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
        # k_new = cart_to_frac_k(np.dot(W, frac_to_cart_k(k0, B)), B)
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
    symmetry : KSymmetry
        A :class:`~ksymmetry.KSymmetry` instance providing the crystal
        symmetry (point-group operations, star computation, coordinate
        conversions).

    Methods
    -------
    optimize_kpath(kpath_frac, method='brute', ...)
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
            If True, also penalise non-uniform segment directions.

        Returns
        -------
        dict with keys:
            kpath_cart, kpath_frac, segment_lengths, total_length,
            original_length, improvement (%), n_combinations, stars.
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

        # original_length = sum(
        #     kpath_segment_length(kpath[i], kpath[i + 1], self._rlat)
        #     for i in range(len(kpath) - 1)
        # )
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
    # Internal: brute force search
    # ------------------------------------------------------------------

    def _optimize_kpath_brute(
        self,
        kpath: list[np.ndarray],
        stars: list[np.ndarray],
        parallel_segments: bool = True,
    ) -> dict:
        """Exhaustive search over all combinations of equivalent k-points."""
        n_points = len(kpath)
        initial_path = [np.asarray(kpath[i], dtype=float) for i in range(n_points)]
        best_score = self._path_length(initial_path)
        if parallel_segments:
            best_score += self._segment_uniformity_penalty(initial_path)

        best_indices = [0] * n_points
        n_combinations = 0

        for indices in itertools.product(*[range(len(s)) for s in stars]):
            n_combinations += 1
            current_path = [stars[i][idx] for i, idx in enumerate(indices)]
            current_length = self._path_length(current_path)
            if parallel_segments:
                penalty = self._segment_uniformity_penalty(current_path)
                current_score = current_length + penalty
            else:
                current_score = current_length
            if current_score < best_score:
                best_score = current_score
                best_indices = list(indices)

        optimized_frac = [stars[i][idx] for i, idx in enumerate(best_indices)]
        optimized_cart = [frac_to_cart_k(k, self._rlat) for k in optimized_frac]
        segment_lengths = [
            kpath_segment_length(optimized_frac[i], optimized_frac[i + 1], self._rlat)
            for i in range(n_points - 1)
        ]
        return {
            "kpath_frac": np.array(optimized_frac),
            "kpath_cart": np.array(optimized_cart),
            "segment_lengths": np.array(segment_lengths),
            "total_length": sum(segment_lengths),
            "n_combinations": n_combinations,
            "stars": stars,
            "best_indices": best_indices,
        }

    # ------------------------------------------------------------------
    # Internal: greedy sequential search
    # ------------------------------------------------------------------

    def _optimize_kpath_greedy(
        self,
        kpath: list[np.ndarray],
        stars: list[np.ndarray],
        parallel_segments: bool = True,
    ) -> dict:
        """Greedy sequential optimisation: pick each k-point to minimise the
        segment to the previously chosen point."""
        n_points = len(kpath)
        best_indices = [0] * n_points
        for i in range(n_points):
            if i == 0:
                best_indices[0] = 0
                continue
            prev_k_cart = frac_to_cart_k(stars[i - 1][best_indices[i - 1]], self._rlat)
            best_dist = float("inf")
            best_idx = 0
            for j, k_star in enumerate(stars[i]):
                k_cart = frac_to_cart_k(k_star, self._rlat)
                dist = float(np.linalg.norm(k_cart - prev_k_cart))
                if parallel_segments and i > 1:
                    prev_prev_k_cart = frac_to_cart_k(
                        stars[i - 2][best_indices[i - 2]],
                        self._rlat
                    )
                    seg1 = prev_k_cart - prev_prev_k_cart
                    seg2 = k_cart - prev_k_cart
                    angle = self._angle_between_vectors(seg1, seg2)
                    dist += dist * (1.0 - angle / np.pi) * 0.1
                if j == 0 or dist < best_dist:
                    best_dist = dist
                    best_idx = j
            best_indices[i] = best_idx

        optimized_frac = [stars[i][idx] for i, idx in enumerate(best_indices)]
        optimized_cart = [frac_to_cart_k(k, self._rlat) for k in optimized_frac]
        segment_lengths = [
            kpath_segment_length(optimized_frac[i], optimized_frac[i + 1], self._rlat)
            for i in range(n_points - 1)
        ]
        return {
            "kpath_frac": np.array(optimized_frac),
            "kpath_cart": np.array(optimized_cart),
            "segment_lengths": np.array(segment_lengths),
            "total_length": sum(segment_lengths),
            "n_combinations": sum(len(s) for s in stars),
            "stars": stars,
            "best_indices": best_indices,
        }

    # ------------------------------------------------------------------
    # Internal: greedy global search
    # ------------------------------------------------------------------

    def _optimize_kpath_greedy_global(
        self,
        kpath: list[np.ndarray],
        stars: list[np.ndarray],
        parallel_segments: bool = True,
    ) -> dict:
        """Greedy global: try every possible starting k-point, then greedily
        extend, and keep the overall shortest path."""
        n_points = len(kpath)
        if n_points == 2:
            return self._optimize_kpath_greedy(kpath, stars, parallel_segments)

        best_indices = None
        best_length = float("inf")
        for start_idx in range(len(stars[0])):
            current_indices = [start_idx]
            current_path = [stars[0][start_idx]]
            for i in range(1, n_points):
                prev_k_cart = frac_to_cart_k(current_path[-1], self._rlat)
                best_local_idx = 0
                best_local_dist = float("inf")
                for j in range(len(stars[i])):
                    k_cart = frac_to_cart_k(stars[i][j], self._rlat)
                    dist = float(np.linalg.norm(k_cart - prev_k_cart))
                    if parallel_segments and len(current_path) > 1:
                        prev_prev_k_cart = frac_to_cart_k(current_path[-2], self._rlat)
                        seg1 = prev_k_cart - prev_prev_k_cart
                        seg2 = k_cart - prev_k_cart
                        angle = self._angle_between_vectors(seg1, seg2)
                        dist += dist * (1.0 - angle / np.pi) * 0.1
                    if dist < best_local_dist:
                        best_local_dist = dist
                        best_local_idx = j
                current_indices.append(best_local_idx)
                current_path.append(stars[i][best_local_idx])
            current_length = self._path_length(current_path)
            if current_length < best_length:
                best_length = current_length
                best_indices = list(current_indices)

        optimized_frac = [stars[i][idx] for i, idx in enumerate(best_indices)]
        optimized_cart = [frac_to_cart_k(k, self._rlat) for k in optimized_frac]
        segment_lengths = [
            kpath_segment_length(optimized_frac[i], optimized_frac[i + 1], self._rlat)
            for i in range(n_points - 1)
        ]
        return {
            "kpath_frac": np.array(optimized_frac),
            "kpath_cart": np.array(optimized_cart),
            "segment_lengths": np.array(segment_lengths),
            "total_length": sum(segment_lengths),
            "n_combinations": len(stars[0]) * sum(len(s) for s in stars[1:]),
            "stars": stars,
            "best_indices": best_indices,
        }

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

    @staticmethod
    def _angle_between_vectors(v1: np.ndarray, v2: np.ndarray) -> float:
        """Angle between two vectors in radians."""
        n1 = float(np.linalg.norm(v1))
        n2 = float(np.linalg.norm(v2))
        if n1 == 0.0 or n2 == 0.0:
            return 0.0
        cos_a = float(np.dot(v1, v2) / (n1 * n2))
        return float(np.arccos(np.clip(cos_a, -1.0, 1.0)))

    def _fold_k_to_first_bz(self, k_frac: Sequence) -> np.ndarray:
        """Fold a k-point to the first Brillouin zone [-0.5, 0.5)."""
        k = np.asarray(k_frac, dtype=float)
        return k - np.floor(k + 0.5)
