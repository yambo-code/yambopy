#
# License-Identifier: GPL
#
# Copyright (C) 2026 The Yambo Team
#
# This file is part of the yambopy project
#
"""Electronic Berry curvature from Yambo velocity dipoles."""

import numpy as np
import warnings

from yambopy.units import bohr2ang, ha2ev


def compute_electronic_berry_curvature(dipolesdb, electronsdb,
                                       bands_range=None, degen_tol=1.e-5,
                                       convention='faria'):
    """Compute the isolated-band Berry curvature in the Yambo IBZ.

    The default convention follows Eq. (4) of Faria Junior *et al.*,
    New J. Phys. 24, 083004 (2022),

    ``Omega_gamma = 2 Im[v_i(n,m) v_j(m,n)]/(E_n-E_m)^2``.

    Yambo stores ``DIP_v`` in the velocity-gauge unit ``Hartree*bohr``
    (atomic velocity, since ``hbar=1``). Yambopy exposes the electronic
    energies in eV, so they are converted back to Hartree before the sum.
    The native result is consequently in ``bohr^2``.

    Parameters
    ----------
    dipolesdb : YamboDipolesDB
        Velocity database loaded with ``dip_type='v'``, ``expand=False`` and
        ``project=False``. The Yambo run must contain ``DipBandsAll``.
    electronsdb : YamboElectronsDB
        Electronic energies on the same irreducible k-point grid.
    bands_range : sequence of two int, optional
        One-based inclusive range of output bands. All stored dipole bands
        remain in the intermediate-state sum.
    degen_tol : float, optional
        Degeneracy tolerance in eV. Affected scalar band values are invalid.
    convention : str, optional
        ``'faria'`` follows the cited paper. ``'wannier90'`` uses the opposite
        sign associated with ``A=<u|i grad_k|u>`` and ``Omega=curl_k A``.

    Returns
    -------
    result : dict
        ``Omega_bohr2`` and ``Omega_ang2`` have layout
        ``(3,nk,nspin,ntarget)``. ``valid`` marks isolated target bands.
        Values with ``valid=False`` omit near-degenerate partners and must
        only be used as diagnostics. Converge the intermediate-band sum and
        use energies and velocities from the same Hamiltonian.
    """
    # Check that the database contains the complete velocity matrix
    if dipolesdb.dip_type != 'v':
        raise ValueError("Berry curvature requires dip_type='v'")
    if dipolesdb.dip_bands_ordered or not hasattr(dipolesdb, 'dipoles_full'):
        raise ValueError("Berry curvature requires a DipBandsAll database")
    if not np.isfinite(degen_tol) or degen_tol < 0:
        raise ValueError("degen_tol must be finite and non-negative")

    signs = {'faria': 1., 'wannier90': -1.}
    if convention not in signs:
        raise ValueError("convention must be 'faria' or 'wannier90'")

    # Restore the spin axis suppressed by YamboDipolesDB when spin == 1
    velocity = np.asarray(dipolesdb.dipoles_full)
    if dipolesdb.spin == 1:
        velocity = velocity[None, ...]
    if (velocity.ndim != 5 or velocity.shape[2] != 3 or
            velocity.shape[-2] != velocity.shape[-1] or
            velocity.shape[0] != dipolesdb.spin):
        raise ValueError("unexpected velocity layout; expected (spin,k,3,band,band)")

    # Match the electronic energies to the stored dipole window
    energies_ev = np.asarray(electronsdb.eigenvalues_ibz)
    if energies_ev.ndim == 2:
        energies_ev = energies_ev[None, ...]
    if energies_ev.ndim != 3 or energies_ev.shape[:2] != velocity.shape[:2]:
        raise ValueError("electron and dipole spin/k-point grids do not match")
    if hasattr(dipolesdb, 'lattice'):
        for name in ('lat', 'alat', 'iku_kpoints'):
            dip_grid = np.asarray(getattr(dipolesdb.lattice, name))
            if name == 'iku_kpoints' and hasattr(dipolesdb.lattice, 'get_ibz_kpoints'):
                dip_grid = np.asarray(dipolesdb.lattice.get_ibz_kpoints())
            elec_grid = np.asarray(getattr(electronsdb, name))
            if dip_grid.shape != elec_grid.shape or not np.allclose(
                    dip_grid, elec_grid, atol=1.e-7, rtol=1.e-6):
                raise ValueError("electron and dipole lattice/k-point order do not match")

    first = dipolesdb.bands_range[0]
    nstored = velocity.shape[-1]
    energies_ha = energies_ev[..., first-1:first-1+nstored] / ha2ev
    if energies_ha.shape[-1] != nstored:
        raise ValueError("electron database does not cover the dipole band window")

    # Select output bands without truncating the sum over stored bands
    if bands_range is None:
        target_first, target_last = first, first+nstored-1
    else:
        if len(bands_range) != 2 or any(int(b) != b for b in bands_range):
            raise ValueError("bands_range must contain two integer band numbers")
        target_first, target_last = map(int, bands_range)
    if not (first <= target_first <= target_last < first+nstored):
        raise ValueError("bands_range lies outside the DipBandsAll window")
    target = np.arange(target_first-first, target_last-first+1)

    # Yambopy stores [alpha,n,m] = <m|v_alpha|n>
    v_mn = velocity[:, :, :, target, :].transpose(0, 1, 3, 4, 2)
    if not np.all(np.isfinite(v_mn)) or not np.all(np.isfinite(energies_ha)):
        raise ValueError("velocity matrix and energies must be finite")
    delta = energies_ha[:, :, target, None]-energies_ha[:, :, None, :]
    keep = np.abs(delta) > degen_tol/ha2ev
    keep &= (target[:, None] != np.arange(nstored)[None, :])[None, None]

    # Sum over states: (i,j,gamma) = (y,z,x), (z,x,y), (x,y,z)
    denominator = np.where(keep, delta**2, np.inf)
    omega = np.empty((velocity.shape[0], velocity.shape[1], len(target), 3))
    for gamma, (i, j) in enumerate(((1, 2), (2, 0), (0, 1))):
        omega[..., gamma] = 2*np.sum(
            np.imag(np.conj(v_mn[..., i])*v_mn[..., j])/denominator,
            axis=-1)

    # A scalar band curvature is not assigned across a degeneracy
    valid = np.sum(~keep, axis=-1) == 1
    if not np.all(valid):
        warnings.warn("Berry sums omit near-degenerate partners; use only "
                      "bands where valid=True", RuntimeWarning, stacklevel=2)
    omega = signs[convention]*omega.transpose(3, 1, 0, 2)

    return {
        'Omega_bohr2': omega,
        'Omega_ang2': omega*bohr2ang**2,
        'valid': valid.transpose(1, 0, 2),
        'bands_range': (target_first, target_last),
        'convention': convention,
        'formula_sign': int(signs[convention]),
        'velocity_unit': 'Hartree*bohr (atomic velocity, hbar=1)',
        'energy_input_unit': 'eV (YamboElectronsDB)',
        'energy_internal_unit': 'Hartree',
    }


def expand_berry_curvature(lattice, berry):
    """Expand an isolated-band Berry curvature from the IBZ to the full BZ.

    ``lattice`` must be a :class:`YamboLatticeDB` loaded with ``Expand=True``.
    Berry curvature is a time-reversal-odd axial vector: a unitary operation
    ``R`` gives ``det(R) R Omega``, while an antiunitary operation adds a
    minus sign. The band-validity mask is copied to every symmetry image.
    """
    omega_ibz = np.asarray(berry['Omega_bohr2'])
    valid_ibz = np.asarray(berry['valid'])
    if omega_ibz.ndim != 4 or omega_ibz.shape[0] != 3:
        raise ValueError("Omega_bohr2 must have layout (3,nk,nspin,nband)")
    if valid_ibz.shape != omega_ibz.shape[1:]:
        raise ValueError("Berry-curvature values and valid mask do not match")
    if omega_ibz.shape[2] != 1:
        raise NotImplementedError("full-BZ Berry expansion requires one storage-spin "
                                  "channel; collinear spin-channel mapping is not implemented")
    if not hasattr(lattice, 'kpoints_indexes'):
        raise ValueError("lattice must be loaded with Expand=True")
    if omega_ibz.shape[1] != lattice.ibz_nkpoints:
        raise ValueError("Berry-curvature and lattice IBZ grids do not match")

    # Use Yambopy's established BZ-to-IBZ and symmetry maps
    bz_to_ibz = np.asarray(lattice.kpoints_indexes, dtype=int)
    symmetry = np.asarray(lattice.symmetry_indexes, dtype=int)
    rotations = np.asarray(lattice.sym_car)[symmetry]
    time_reversal = np.asarray(lattice.time_rev_list, dtype=bool)[symmetry]

    axial = np.linalg.det(rotations)[:, None, None]*rotations
    axial[time_reversal] *= -1
    omega_bz = np.einsum(
        'kij,jksb->iksb', axial, omega_ibz[:, bz_to_ibz], optimize=True)

    result = dict(berry)
    result['Omega_bohr2_BZ'] = omega_bz
    result['Omega_ang2_BZ'] = omega_bz*bohr2ang**2
    result['valid_BZ'] = valid_ibz[bz_to_ibz]
    result['kpoints_BZ'] = np.asarray(lattice.car_kpoints)
    return result


def compute_electronic_berry_curvature_BZ(lattice, dipolesdb, electronsdb,
                                          bands_range=None, degen_tol=1.e-5,
                                          convention='faria'):
    """Compute the band Berry curvature in the IBZ and expand it to the BZ."""
    berry = compute_electronic_berry_curvature(
        dipolesdb, electronsdb, bands_range=bands_range,
        degen_tol=degen_tol, convention=convention)
    return expand_berry_curvature(lattice, berry)
