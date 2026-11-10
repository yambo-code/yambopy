"""
lattice_data.py
===============
Lattice and high-symmetry k-point data for the QE Bravais lattices
(ibrav 1-11, -3, -5, -9), with the Setyawan-Curtarolo special points and
default band-structure paths:

* :func:`get_lattice_data` — cell, variant, special points and default
  path for an ibrav and its lattice parameters;
* :func:`lattice_data_for_cell` — the same for an arbitrary cell, keeping
  the cell and expressing the special points in its reduced coordinates;
* :func:`find_basis_change` — integer change of basis between two cells
  of the same lattice.

Special points are in reduced reciprocal coordinates of the returned cell,
named following :mod:`yambopy.bz.labels` (``G`` for Γ, ``Sigma1`` for Σ₁).
"""

from itertools import product
from math import cos, isclose, radians, sqrt, tan
from typing import Sequence

import numpy as np


#: Required lattice parameters for each supported ibrav.
IBRAV_REQUIRED_PARAMETERS: dict[int, list[str]] = {
    1: ["a"],
    2: ["a"],
    3: ["a"],
    -3: ["a"],
    4: ["a", "c"],
    5: ["a", "gamma"],
    -5: ["a", "gamma"],
    6: ["a", "c"],
    7: ["a", "c"],
    8: ["a", "b", "c"],
    9: ["a", "b", "c"],
    -9: ["a", "b", "c"],
    10: ["a", "b", "c"],
    11: ["a", "b", "c"],
}


# =============================================================================
# Lattice data
# =============================================================================

#: ibrav values whose Setyawan-Curtarolo special points assume ordered
#: lattice constants, and the parameters that must be sorted.  (Simple
#: orthorhombic, ibrav=8, is valid for any order: sorting would only
#: rename X, Y, Z, so it is left as given.)
_ORDERED_PARAMETERS: dict[int, tuple[str, ...]] = {
    9: ("a", "b"),
    -9: ("a", "b"),
    10: ("a", "b", "c"),
    11: ("a", "b", "c"),
}


def get_lattice_data(
    ibrav: int,
    parameters: dict[str, float],
) -> tuple[np.ndarray, str, dict, str]:
    """
    Return lattice vectors and high-symmetry k-point information for a
    given QE Bravais-lattice index.

    For the centred orthorhombic lattices (ibrav 9, -9, 10, 11) the
    Setyawan-Curtarolo formulas assume ordered lattice constants
    (a < b, or a < b < c).  Other orderings are handled by computing the
    points for the sorted constants and expressing them in the reduced
    coordinates of the requested cell, which is returned as QE defines it.

    Parameters
    ----------
    ibrav : int
        QE Bravais-lattice index (1-11, plus -3, -5, -9).
    parameters : dict
        Lattice parameters dict (e.g. ``{'a': 3.615}``).  Required keys
        depend on *ibrav* (see :data:`IBRAV_REQUIRED_PARAMETERS`).  Lengths
        are in any unit (the cell uses the same one); the rhombohedral
        angle ``gamma`` is in degrees.

    Returns
    -------
    cell : ndarray, shape (3, 3)
        Direct-lattice basis vectors as rows.
    variant : str
        Bravais-lattice variant label (e.g. ``'CUB'``, ``'FCC'``).
    high_symmetry_points : dict[str, list[float] or ndarray]
        Reduced coordinates of the high-symmetry k-points, in the basis of
        the reciprocal of *cell*.
    default_path : str
        Default piecewise band-structure path (comma-separated sections).
    """
    data = _standard_lattice_data(ibrav, parameters)

    keys = _ORDERED_PARAMETERS.get(ibrav)
    if keys is None:
        return data
    values = [parameters[k] for k in keys]
    if values == sorted(values):
        return data

    cell = data[0]
    sorted_parameters = dict(parameters, **dict(zip(keys, sorted(values))))
    std_cell, variant, points, default_path = _standard_lattice_data(ibrav, sorted_parameters)

    M = find_basis_change(std_cell, cell)
    if M is None:  # cannot happen: both cells describe the same lattice
        raise RuntimeError(f"No basis change found for ibrav={ibrav}, {parameters}.")

    points = {label: M @ np.asarray(k, dtype=float) for label, k in points.items()}
    return cell, variant, points, default_path


def _standard_lattice_data(
    ibrav: int,
    parameters: dict[str, float],
) -> tuple[np.ndarray, str, dict, str]:
    """
    :func:`get_lattice_data` without the handling of unordered lattice
    constants: the Setyawan-Curtarolo special points are only correct when
    the constants are in the order they assume.

    Parameters
    ----------
    ibrav : int
        QE Bravais-lattice index (1-11, plus -3, -5, -9).
    parameters : dict
        Lattice parameters dict (e.g. ``{'a': 3.615}``).  Required keys
        depend on *ibrav* (see :data:`IBRAV_REQUIRED_PARAMETERS`).  Lengths
        are in any unit (the cell uses the same one); the rhombohedral
        angle ``gamma`` is in degrees.

    Returns
    -------
    cell : ndarray, shape (3, 3)
        Direct-lattice basis vectors as rows.
    variant : str
        Bravais-lattice variant label (e.g. ``'CUB'``, ``'FCC'``).
    high_symmetry_points : dict[str, list[float]]
        Fractional coordinates of high-symmetry k-points.
    default_path : str
        Default piecewise band-structure path (comma-separated sections).
    """
    if not isinstance(parameters, dict):
        raise TypeError(f"{parameters} is not a dictionary")

    if ibrav not in IBRAV_REQUIRED_PARAMETERS:
        raise ValueError(f"ibrav: {ibrav} not supported.")

    for p in IBRAV_REQUIRED_PARAMETERS[ibrav]:
        if p not in parameters:
            raise ValueError(
                f"ibrav: {ibrav} lattice needs parameter: '{p}'"
            )

    for param, value in parameters.items():
        if param in ("a", "b", "c") and value <= 0:
            raise ValueError(f"{param} parameter must be positive.")

    a = parameters.get("a", 1.0)
    b = parameters.get("b", 1.0)
    c = parameters.get("c", 1.0)
    gamma = parameters.get("gamma", 90.0)  # degrees; only ibrav ±5 uses it

    match ibrav:
        # --- CUB ---
        case 1:
            v1 = [a, 0, 0]
            v2 = [0, a, 0]
            v3 = [0, 0, a]
            cell = np.array([v1, v2, v3])
            variant = "CUB"
            high_symmetry_points = {
                "G": [0, 0, 0],
                "M": [1 / 2, 1 / 2, 0],
                "R": [1 / 2, 1 / 2, 1 / 2],
                "X": [0, 1 / 2, 0],
            }
            default_path = "GXMGRX,MR"

        # --- FCC ---
        case 2:
            v1 = [-a / 2, 0, a / 2]
            v2 = [0, a / 2, a / 2]
            v3 = [-a / 2, a / 2, 0]
            cell = np.array([v1, v2, v3])
            variant = "FCC"
            high_symmetry_points = {
                "G": [0, 0, 0],
                "K": [-3 / 8, 3 / 8, 0],
                "L": [0, 1 / 2, 0],
                "U": [0, 5 / 8, 3 / 8],
                "W": [-1 / 4, 1 / 2, 1 / 4],
                "X": [0, 1 / 2, 1 / 2],
            }
            default_path = "GXWKGLUWLK,UX"

        # --- BCC ---
        case 3:
            v1 = [a / 2, a / 2, a / 2]
            v2 = [-a / 2, a / 2, a / 2]
            v3 = [-a / 2, -a / 2, a / 2]
            cell = np.array([v1, v2, v3])
            variant = "BCC"
            high_symmetry_points = {
                "G": [0, 0, 0],
                "H": [1 / 2, 1 / 2, -1 / 2],
                "P": [3 / 4, 1 / 4, -1 / 4],
                "N": [1 / 2, 0, -1 / 2],
            }
            default_path = "GHNGPH,PN"

        # --- BCC (alternative) ---
        case -3:
            v1 = [-a / 2, a / 2, a / 2]
            v2 = [a / 2, -a / 2, a / 2]
            v3 = [a / 2, a / 2, -a / 2]
            cell = np.array([v1, v2, v3])
            variant = "BCC"
            high_symmetry_points = {
                "G": [0, 0, 0],
                "H": [1 / 2, -1 / 2, 1 / 2],
                "P": [1 / 4, 1 / 4, 1 / 4],
                "N": [0, 0, 1 / 2],
            }
            default_path = "GHNGPH,PN"

        # --- HEX ---
        case 4:
            v1 = [a, 0, 0]
            v2 = [-a / 2, a * sqrt(3) / 2, 0]
            v3 = [0, 0, c]
            cell = np.array([v1, v2, v3])
            variant = "HEX"
            high_symmetry_points = {
                "G": [0, 0, 0],
                "A": [0, 0, 1 / 2],
                # M = (1/2, 0, 0) is the midpoint of the BZ edge joining the
                # corners (2/3, -1/3, 0) and (1/3, 1/3, 0); either works as K.
                "H": [2 / 3, -1 / 3, 1 / 2],
                "K": [2 / 3, -1 / 3, 0],
                "L": [1 / 2, 0, 1 / 2],
                "M": [1 / 2, 0, 0],
            }
            default_path = "GMKGALHA,LM,KH"

        # --- RHL (rhombohedral) ---
        case 5:
            cs = cos(radians(gamma))
            tx = sqrt((1 - cs) / 2)
            ty = sqrt((1 - cs) / 6)
            tz = sqrt((1 + 2 * cs) / 3)
            v1 = [a * tx, -a * ty, a * tz]
            v2 = [0, a * 2 * ty, a * tz]
            v3 = [-a * tx, -a * ty, a * tz]
            cell = np.array([v1, v2, v3])

            if 0 < gamma < 90:
                variant = "RHL1"
                eta = (1 + 4 * cos(radians(gamma))) / (2 + 4 * cos(radians(gamma)))
                nu = 3 / 4 - eta / 2
                high_symmetry_points = {
                    "G": [0, 0, 0],
                    "B": [eta, 1 / 2, 1 - eta],
                    "B1": [1 / 2, 1 - eta, eta - 1],
                    "F": [1 / 2, 1 / 2, 0],
                    "L": [1 / 2, 0, 0],
                    "L1": [0, 0, -1 / 2],
                    "P": [eta, nu, nu],
                    "P1": [1 - nu, 1 - nu, 1 - eta],
                    "P2": [nu, nu, eta - 1],
                    "Q": [1 - nu, nu, 0],
                    "X": [nu, 0, -nu],
                    "Z": [1 / 2, 1 / 2, 1 / 2],
                }
                default_path = "GLB1,BZGX,QFP1Z,LP"
            elif 90 < gamma < 120:
                variant = "RHL2"
                eta = 1 / (2 * tan(radians(gamma / 2)) ** 2)
                nu = 3 / 4 - eta / 2
                high_symmetry_points = {
                    "G": [0, 0, 0],
                    "F": [0, 1 / 2, -1 / 2],
                    "L": [0, 1 / 2, 0],
                    "P": [1 - nu, 1 - nu, -nu],
                    "P1": [nu - 1, nu, nu - 1],
                    "Q": [eta, eta, eta],
                    "Q1": [-eta, 1 - eta, -eta],
                    "Z": [1 / 2, 1 / 2, -1 / 2],
                }
                default_path = "GPZQGFP1Q1LZ"
            else:
                raise ValueError("Invalid gamma value for ibrav=5. Must be in (0,90) or (90,120).")

        # --- RHL (alternative) ---
        case -5:
            a /= sqrt(3)
            cs = cos(radians(gamma))
            ty = sqrt((1 - cs) / 6)
            tz = sqrt((1 + 2 * cs) / 3)
            u_ = tz - 2 * sqrt(2) * ty
            v_ = tz + sqrt(2) * ty
            v1 = [a * u_, a * v_, a * v_]
            v2 = [a * v_, a * u_, a * v_]
            v3 = [a * v_, a * v_, a * u_]
            cell = np.array([v1, v2, v3])

            if 0 < gamma < 90:
                variant = "RHL1"
                eta = (1 + 4 * cos(radians(gamma))) / (2 + 4 * cos(radians(gamma)))
                nu = 3 / 4 - eta / 2
                high_symmetry_points = {
                    "G": [0, 0, 0],
                    "B": [eta, 1 / 2, 1 - eta],
                    "B1": [1 / 2, 1 - eta, eta - 1],
                    "F": [1 / 2, 1 / 2, 0],
                    "L": [1 / 2, 0, 0],
                    "L1": [0, 0, -1 / 2],
                    "P": [eta, nu, nu],
                    "P1": [1 - nu, 1 - nu, 1 - eta],
                    "P2": [nu, nu, eta - 1],
                    "Q": [1 - nu, nu, 0],
                    "X": [nu, 0, -nu],
                    "Z": [1 / 2, 1 / 2, 1 / 2],
                }
                default_path = "GLB1,BZGX,QFP1Z,LP"
            elif 90 < gamma < 120:
                variant = "RHL2"
                eta = 1 / (2 * tan(radians(gamma / 2)) ** 2)
                nu = 3 / 4 - eta / 2
                high_symmetry_points = {
                    "G": [0, 0, 0],
                    "F": [0, 1 / 2, -1 / 2],
                    "L": [0, 1 / 2, 0],
                    "P": [1 - nu, 1 - nu, -nu],
                    "P1": [nu - 1, nu, nu - 1],
                    "Q": [eta, eta, eta],
                    "Q1": [-eta, 1 - eta, -eta],
                    "Z": [1 / 2, 1 / 2, -1 / 2],
                }
                default_path = "GPZQGFP1Q1LZ"
            else:
                raise ValueError("Invalid gamma value for ibrav=-5. Must be in (0,90) or (90,120).")

        # --- TET ---
        case 6:
            v1 = [a, 0, 0]
            v2 = [0, a, 0]
            v3 = [0, 0, c]
            cell = np.array([v1, v2, v3])
            variant = "TET"
            high_symmetry_points = {
                "G": [0, 0, 0],
                "A": [1 / 2, 1 / 2, 1 / 2],
                "M": [1 / 2, 1 / 2, 0],
                "R": [0, 1 / 2, 1 / 2],
                "X": [0, 1 / 2, 0],
                "Z": [0, 0, 1 / 2],
            }
            default_path = "GXMGZRAZ,XR,MA"

        # --- BCT (body-centred tetragonal) ---
        case 7:
            v1 = [a / 2, -a / 2, c / 2]
            v2 = [a / 2, a / 2, c / 2]
            v3 = [-a / 2, -a / 2, c / 2]
            cell = np.array([v1, v2, v3])

            if isclose(a, c, rel_tol=1e-5):
                raise ValueError("For ibrav=7, a and c must differ. Choose either c < a or c > a.")
            elif c < a:
                variant = "BCT1"
                eta = (1 + c ** 2 / a ** 2) / 4
                high_symmetry_points = {
                    "G": [0, 0, 0],
                    "M": [1 / 2, 1 / 2, -1 / 2],
                    "N": [1 / 2, 1 / 2, 0],
                    "P": [1 / 4, 3 / 4, -1 / 4],
                    "X": [0, 1 / 2, -1 / 2],
                    "Z": [eta, eta, eta],
                    "Z1": [1 - eta, 1 - eta, -eta],
                }
                default_path = "GXMGZPNZ1M,XP"
            else:
                variant = "BCT2"
                eta = (1 + a ** 2 / c ** 2) / 4
                zeta = a ** 2 / (2 * c ** 2)
                high_symmetry_points = {
                    "G": [0, 0, 0],
                    "N": [1 / 2, 1 / 2, 0],
                    "P": [1 / 4, 3 / 4, -1 / 4],
                    "Sigma": [eta, eta, -eta],
                    "Sigma1": [1 - eta, 1 - eta, eta],
                    "X": [0, 1 / 2, -1 / 2],
                    "Y": [zeta, 1 / 2, -1 / 2],
                    "Y1": [1 / 2, 1 - zeta, zeta],
                    "Z": [1 / 2, 1 / 2, 1 / 2],
                }
                default_path = "GXYSigmaGZSigma1NPY1Z,XP"

        # --- ORC (simple orthorhombic) ---
        case 8:
            v1 = [a, 0, 0]
            v2 = [0, b, 0]
            v3 = [0, 0, c]
            cell = np.array([v1, v2, v3])
            variant = "ORC"
            high_symmetry_points = {
                "G": [0, 0, 0],
                "R": [1 / 2, 1 / 2, 1 / 2],
                "S": [1 / 2, 1 / 2, 0],
                "T": [0, 1 / 2, 1 / 2],
                "U": [1 / 2, 0, 1 / 2],
                "X": [1 / 2, 0, 0],
                "Y": [0, 1 / 2, 0],
                "Z": [0, 0, 1 / 2],
            }
            default_path = "GXSYGZURTZ,YT,UX,SR"

        # --- ORCC (C-face-centred orthorhombic) ---
        case 9:
            v1 = [a / 2, b / 2, 0]
            v2 = [-a / 2, b / 2, 0]
            v3 = [0, 0, c]
            cell = np.array([v1, v2, v3])
            variant = "ORCC"
            zeta = (1 + a ** 2 / b ** 2) / 4
            high_symmetry_points = {
                "G": [0, 0, 0],
                "A": [zeta, -zeta, 1 / 2],
                "A1": [1 - zeta, zeta, 1 / 2],
                "R": [1 / 2, 0, 1 / 2],
                "S": [1 / 2, 0, 0],
                "T": [1 / 2, 1 / 2, 1 / 2],
                "X": [zeta, -zeta, 0],
                "X1": [1 - zeta, zeta, 0],
                "Y": [1 / 2, 1 / 2, 0],
                "Z": [0, 0, 1 / 2],
            }
            default_path = "GXSRAZGYX1A1TY,ZT"

        # --- ORCC (alternative, A-face) ---
        case -9:
            v1 = [a / 2, -b / 2, 0]
            v2 = [a / 2, b / 2, 0]
            v3 = [0, 0, c]
            cell = np.array([v1, v2, v3])
            variant = "ORCC"
            zeta = (1 + a ** 2 / b ** 2) / 4
            high_symmetry_points = {
                "G": [0, 0, 0],
                "A": [zeta, zeta, 1 / 2],
                "A1": [-zeta, 1 - zeta, 1 / 2],
                "R": [0, 1 / 2, 1 / 2],
                "S": [0, 1 / 2, 0],
                "T": [-1 / 2, 1 / 2, 1 / 2],
                "X": [zeta, zeta, 0],
                "X1": [-zeta, 1 - zeta, 0],
                "Y": [-1 / 2, 1 / 2, 0],
                "Z": [0, 0, 1 / 2],
            }
            default_path = "GXSRAZGYX1A1TY,ZT"

        # --- ORCF (F-face-centred orthorhombic) ---
        case 10:
            v1 = [a / 2, 0, c / 2]
            v2 = [a / 2, b / 2, 0]
            v3 = [0, b / 2, c / 2]
            cell = np.array([v1, v2, v3])

            ia = 1 / a ** 2
            ibc = 1 / b ** 2 + 1 / c ** 2
            if isclose(ia, ibc, rel_tol=1e-5):
                variant = "ORCF3"
            elif ia < ibc:
                variant = "ORCF2"
            else:
                variant = "ORCF1"

            if variant in ("ORCF1", "ORCF3"):
                zeta = (1 + a ** 2 / b ** 2 - a ** 2 / c ** 2) / 4
                eta = (1 + a ** 2 / b ** 2 + a ** 2 / c ** 2) / 4
                high_symmetry_points = {
                    "G": [0, 0, 0],
                    "A": [1 / 2 + zeta, zeta, 1 / 2],
                    "A1": [1 / 2 - zeta, 1 - zeta, 1 / 2],
                    "L": [1 / 2, 1 / 2, 1 / 2],
                    "T": [1 / 2, 1 / 2, 1],
                    "X": [eta, eta, 0],
                    "X1": [1 - eta, 1 - eta, 1],
                    "Y": [0, 1 / 2, 1 / 2],
                    "Z": [1 / 2, 0, 1 / 2],
                }
                if variant == "ORCF1":
                    default_path = "GYTZGXA1Y,TX1,XAZ,LG"
                else:
                    default_path = "GYTZGXA1Y,XAZ,LG"
            else:
                eta = (1 + a ** 2 / b ** 2 - a ** 2 / c ** 2) / 4
                delta = (1 + b ** 2 / a ** 2 - b ** 2 / c ** 2) / 4
                phi = (1 + c ** 2 / b ** 2 - c ** 2 / a ** 2) / 4
                high_symmetry_points = {
                    "G": [0, 0, 0],
                    "C": [1 / 2 - eta, 1 - eta, 1 / 2],
                    "C1": [1 / 2 + eta, eta, 1 / 2],
                    "D": [1 / 2, 1 - delta, 1 / 2 - delta],
                    "D1": [1 / 2, delta, 1 / 2 + delta],
                    "L": [1 / 2, 1 / 2, 1 / 2],
                    "H": [1 / 2 - phi, 1 / 2, 1 - phi],
                    "H1": [1 / 2 + phi, 1 / 2, phi],
                    "X": [1 / 2, 1 / 2, 0],
                    "Y": [0, 1 / 2, 1 / 2],
                    "Z": [1 / 2, 0, 1 / 2],
                }
                default_path = "GYCDXGZD1HC,C1Z,XH1,HY,LG"

        # --- ORCI (I-body-centred orthorhombic) ---
        case 11:
            v1 = [a / 2, b / 2, c / 2]
            v2 = [-a / 2, b / 2, c / 2]
            v3 = [-a / 2, -b / 2, c / 2]
            cell = np.array([v1, v2, v3])
            variant = "ORCI"
            mu = (a ** 2 + b ** 2) / (4 * c ** 2)
            delta = (b ** 2 - a ** 2) / (4 * c ** 2)
            zeta = (1 + a ** 2 / c ** 2) / 4
            eta = (1 + b ** 2 / c ** 2) / 4
            high_symmetry_points = {
                "G": [0, 0, 0],
                "L": [1 / 2 - delta, -mu, delta - 1 / 2],
                "L1": [1 / 2 + delta, mu, -1 / 2 - delta],
                "L2": [1 - mu, 1 / 2 - delta, mu],
                "R": [1 / 2, 0, 0],
                "S": [1 / 2, 1 / 2, 0],
                "T": [1 / 2, 0, -1 / 2],
                "W": [3 / 4, 1 / 4, -1 / 4],
                "X": [zeta, -zeta, -zeta],
                "X1": [1 - zeta, zeta, zeta],
                "Y": [eta, eta, -eta],
                "Y1": [1 - eta, 1 - eta, eta],
                "Z": [1 / 2, 1 / 2, 1 / 2],
            }
            default_path = "GXLTWRX1ZGYSW,L1Y,Y1Z"

        case _:
            raise ValueError(f"ibrav: {ibrav} not supported.")

    return cell, variant, high_symmetry_points, default_path


# =============================================================================
# Special points for an arbitrary cell
# =============================================================================

#: Space group of a one-atom lattice (its holohedry) -> QE ibrav.
_HOLOHEDRY_TO_IBRAV = {
    221: 1,   # Pm-3m   cubic P
    225: 2,   # Fm-3m   cubic F
    229: 3,   # Im-3m   cubic I
    191: 4,   # P6/mmm  hexagonal
    166: 5,   # R-3m    rhombohedral
    123: 6,   # P4/mmm  tetragonal P
    139: 7,   # I4/mmm  tetragonal I
    47: 8,    # Pmmm    orthorhombic P
    65: 9,    # Cmmm    orthorhombic C
    69: 10,   # Fmmm    orthorhombic F
    71: 11,   # Immm    orthorhombic I
}


def find_basis_change(
    std_cell: np.ndarray,
    cell: np.ndarray,
    rtol: float = 1e-4,
) -> np.ndarray | None:
    """
    Find an integer matrix ``M`` with ``|det M| = 1`` such that
    ``cell = M @ std_cell @ R.T`` for some orthogonal ``R``, i.e. both
    cells (basis vectors as rows) span the same lattice up to a rotation.

    The search matches metric tensors, ``cell @ cell.T == M G_std M.T``.
    A lattice vector ``v = n @ std_cell`` has ``|n_i| <= |v| |b_i|``, which
    bounds the coefficients, so the search is exhaustive.  Proper changes
    of basis (``det M = +1``) are preferred.  Returns *None* if none exists.

    Reduced reciprocal coordinates then transform as ``k_cell = M @ k_std``.

    Parameters
    ----------
    std_cell, cell : array-like, shape (3, 3)
        The two cells, basis vectors as rows, in the same length unit.
    rtol : float
        Tolerance on the metric tensor entries, relative to the largest
        squared length of *cell*'s basis vectors.
    """
    std_cell = np.asarray(std_cell, dtype=float)
    cell = np.asarray(cell, dtype=float)
    G_std = std_cell @ std_cell.T
    G = cell @ cell.T
    atol = rtol * float(np.max(np.abs(np.diag(G))))

    rcell_std = np.linalg.inv(std_cell.T)
    max_coeff = int(np.ceil(
        np.sqrt(np.max(np.diag(G))) * np.max(np.linalg.norm(rcell_std, axis=1)) + 1e-6
    ))
    rng = range(-max_coeff, max_coeff + 1)
    vecs = np.array([v for v in product(rng, repeat=3) if any(v)], dtype=float)
    norms = np.einsum("ij,jk,ik->i", vecs, G_std, vecs)
    cands = [vecs[np.abs(norms - G[i, i]) < atol] for i in range(3)]

    improper = None
    for v0 in cands[0]:
        for v1 in cands[1]:
            if abs(v0 @ G_std @ v1 - G[0, 1]) >= atol:
                continue
            for v2 in cands[2]:
                if (abs(v0 @ G_std @ v2 - G[0, 2]) >= atol
                        or abs(v1 @ G_std @ v2 - G[1, 2]) >= atol):
                    continue
                M = np.array([v0, v1, v2])
                det = round(np.linalg.det(M))
                if det == 1:
                    return M.astype(int)
                if det == -1 and improper is None:
                    improper = M.astype(int)
    return improper


def lattice_data_for_cell(
    cell: Sequence,
    symprec: float = 1e-5,
) -> tuple[np.ndarray, str, dict, str]:
    """
    Like :func:`get_lattice_data`, but for a user-given cell (basis vectors
    as rows), which is returned unchanged.

    The Bravais lattice is identified with spglib, from the symmetry of a
    one-atom crystal on *cell*, so the result does not depend on the
    cell's orientation or choice of primitive vectors (e.g. hexagonal with
    60 instead of 120 degrees, ibrav -5/-9 orientations, rotated cells).
    The special points of the matching standard cell are then expressed
    in the reduced coordinates of *cell*.

    Parameters
    ----------
    cell : array-like, shape (3, 3)
        Direct-lattice basis vectors as rows.
    symprec : float
        Relative tolerance; spglib gets ``symprec`` times the longest
        basis vector as its (absolute) distance tolerance.
    """
    import warnings

    import spglib

    cell = np.asarray(cell, dtype=float)
    scale = float(np.max(np.linalg.norm(cell, axis=1)))
    with warnings.catch_warnings():
        # spglib >= 2.5 warns about its legacy error handling on every call;
        # failures are still reported by returning None, handled below.
        warnings.simplefilter("ignore", DeprecationWarning)
        dataset = spglib.get_symmetry_dataset(
            (cell, [[0.0, 0.0, 0.0]], [1]), symprec=symprec * scale
        )
    if dataset is None:
        raise ValueError("spglib could not determine the lattice symmetry of the cell.")

    ibrav = _HOLOHEDRY_TO_IBRAV.get(dataset.number)
    if ibrav is None:
        raise ValueError(
            f"Lattice with holohedry {dataset.international} "
            f"(space group {dataset.number}) is not supported "
            "(monoclinic and triclinic lattices have no ibrav data here)."
        )

    # Conventional lattice constants from spglib's standardised cell.
    a, b, c = np.linalg.norm(dataset.std_lattice, axis=1)
    params: dict[str, float] = {"a": a, "b": b, "c": c}
    if ibrav == 5:
        # spglib gives the hexagonal setting; convert to rhombohedral.
        params = {
            "a": float(np.sqrt(a ** 2 / 3 + c ** 2 / 9)),
            "gamma": float(np.degrees(np.arccos(
                (2 * c ** 2 - 3 * a ** 2) / (2 * c ** 2 + 6 * a ** 2)
            ))),
        }
    # Any ordering of a, b, c is fine here: get_lattice_data handles the
    # orderings the Setyawan-Curtarolo formulas do not assume.

    std_cell, variant, points, default_path = get_lattice_data(ibrav, params)

    M = find_basis_change(std_cell, cell, rtol=max(10 * symprec, 1e-6))
    if M is None:
        raise ValueError(
            f"Could not relate the given cell to the standard ibrav={ibrav} "
            f"cell {params}. Pass ibrav and parameters instead."
        )

    points = {label: M @ np.asarray(k, dtype=float) for label, k in points.items()}
    return cell.copy(), variant, points, default_path
