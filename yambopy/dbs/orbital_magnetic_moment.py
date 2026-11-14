#
# License-Identifier: GPL
#
# Copyright (C) 2026 The Yambo Team
#
# This file is part of the yambopy project
#
"""Electronic magnetic moments from Yambo velocity and spin dipoles."""

import numpy as np
import warnings
from yambopy.dbs.spindb import YamboSpinDB
from yambopy.units import ha2ev


FREE_ELECTRON_G = 2.00231930436092


def compute_electronic_magnetic_moment(dipolesdb, electronsdb,
                                       filename='ndb.dipoles',
                                       bands_range=None, degen_tol=1.e-5,
                                       g_factor=FREE_ELECTRON_G, spindb=None):
    """Return orbital, spin and total band magnetic operators.

    All returned operator arrays have the Yambopy generic-operator layout
    ``(3,nk,nspin,nbands,nbands)``.  The physical moment is
    ``mu/mu_B = -(L/hbar + g*sigma/2)``; the Zeeman operator is its opposite.
    Here ``zeeman_g_operator`` is the coefficient of ``mu_B B``, not a
    two-branch splitting g factor. Only the spin part is a full band matrix.
    """
    result = compute_electronic_orbital_moment(
        dipolesdb, electronsdb, bands_range=bands_range,
        degen_tol=degen_tol, return_operator=True)
    if spindb is None:
        spindb = YamboSpinDB.from_db_file(filename, result['bands_range'])
    if tuple(spindb.bands_range) != result['bands_range']:
        raise ValueError("DIP_spin and orbital operator address different bands")
    if not np.isfinite(g_factor):
        raise ValueError("g_factor must be finite")
    spin_sigma = spindb.as_operator()
    if spin_sigma.shape != result['operator'].shape:
        raise ValueError("DIP_spin and orbital-operator grids do not match")
    spin_g = (g_factor/2.)*spin_sigma
    spin_moment = -spin_g
    orbital_g = -result['operator']
    zeeman_g = orbital_g + spin_g
    result['spin_sigma'] = spin_sigma
    result['spin_g_operator'] = spin_g
    result['orbital_g_operator'] = orbital_g
    result['zeeman_g_operator'] = zeeman_g
    result['spin_moment_over_muB'] = spin_moment
    result['total_moment_over_muB'] = -zeeman_g
    result['total_operator'] = result['total_moment_over_muB']
    result['free_electron_g'] = g_factor
    return result


def compute_electronic_orbital_moment(dipolesdb, electronsdb,
                                      bands_range=None, degen_tol=1.e-5,
                                      return_operator=False):
    """Compute diagonal Bloch-band orbital moments from ``DipBandsAll``.

    Yambo stores ``DIP_v`` in atomic velocity units (``Hartree*bohr`` with
    ``hbar=1``), including the nonlocal correction when enabled. Energies from
    ``YamboElectronsDB`` are converted from eV to Hartree, giving
    dimensionless ``L/hbar``.  ``mu_orb/mu_B = -L/hbar``.

    Parameters
    ----------
    dipolesdb : YamboDipolesDB
        Full velocity database generated with ``DipBandsAll`` and loaded with
        ``dip_type='v'``, ``expand=False`` and ``project=False``.
    electronsdb : YamboElectronsDB
        Electronic energies on the same irreducible k-point grid.  Yambopy
        stores these energies in eV.
    bands_range : sequence of two int, optional
        One-based inclusive target-band range, for example ``[23, 30]``.
        All bands in ``dipolesdb`` are still used as intermediate states.
    degen_tol : float, optional
        Energy tolerance in eV.  A target with another band closer than this
        value is marked invalid for the isolated-band formula.
    return_operator : bool, optional
        If True, embed the diagonal moments in the layout consumed by
        ``exciton_X_matelem``. Off-diagonal orbital elements are not computed.

    Returns
    -------
    result : dict
        Contains ``L_over_hbar``, ``mu_orb_over_muB``, ``valid`` and
        ``bands_range``.  It also contains ``operator`` when requested.

    Notes
    -----
    ``L_over_hbar`` is the self-rotation quantity defined by the velocity
    sum, not an atom-projected canonical ``r cross p`` for a nonlocal
    Hamiltonian. Use energies and velocities from the same Hamiltonian and
    converge the intermediate-band sum. Where ``valid=False``, the returned
    finite sum omits near-degenerate partners and is diagnostic only.
    """
    if dipolesdb.dip_type != 'v':
        raise ValueError("orbital magnetic moments require dip_type='v'")
    if dipolesdb.dip_bands_ordered or not hasattr(dipolesdb, 'dipoles_full'):
        raise ValueError("orbital magnetic moments require a DipBandsAll database")
    if not np.isfinite(degen_tol) or degen_tol < 0:
        raise ValueError("degen_tol must be finite and non-negative")

    velocity = np.asarray(dipolesdb.dipoles_full)
    if dipolesdb.spin == 1:
        velocity = velocity[None, ...]
    if (velocity.ndim != 5 or velocity.shape[2] != 3 or
            velocity.shape[-2] != velocity.shape[-1] or
            velocity.shape[0] != dipolesdb.spin):
        raise ValueError("unexpected velocity layout; expected (spin,k,3,band,band)")

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
    energies_ev = energies_ev[..., first-1:first-1+nstored]
    if energies_ev.shape[-1] != nstored:
        raise ValueError("electron database does not cover the dipole band window")

    if bands_range is None:
        target_first, target_last = first, first + nstored - 1
    else:
        if len(bands_range) != 2 or any(int(b) != b for b in bands_range):
            raise ValueError("bands_range must contain two integer band numbers")
        target_first, target_last = map(int, bands_range)
    if not (first <= target_first <= target_last < first + nstored):
        raise ValueError("bands_range lies outside the DipBandsAll window")
    target = np.arange(target_first-first, target_last-first+1)

    # YamboDipolesDB stores [alpha,n,l] = <l|v_alpha|n>.
    v_ln = velocity[:, :, :, target, :].transpose(0, 1, 3, 4, 2)
    if not np.all(np.isfinite(v_ln)) or not np.all(np.isfinite(energies_ev)):
        raise ValueError("velocity matrix and energies must be finite")
    energies_ha = energies_ev / ha2ev
    delta = energies_ha[:, :, target, None] - energies_ha[:, :, None, :]
    keep = np.abs(delta) > degen_tol / ha2ev
    external = target[:, None] == np.arange(nstored)[None, :]
    keep &= ~external[None, None, :, :]

    # L_gamma/hbar = 2 Im[v_i* v_j]/(E_n-E_l).
    angular = np.empty((velocity.shape[0], velocity.shape[1], len(target), 3))
    safe_delta = np.where(keep, delta, np.inf)
    for gamma, (i, j) in enumerate(((1, 2), (2, 0), (0, 1))):
        angular[..., gamma] = 2*np.sum(
            np.imag(np.conj(v_ln[..., i]) * v_ln[..., j]) / safe_delta,
            axis=-1,
        )

    valid = np.sum(~keep, axis=-1) == 1
    if not np.all(valid):
        warnings.warn("Orbital sums omit near-degenerate partners; use only "
                      "bands where valid=True", RuntimeWarning, stacklevel=2)
    angular = angular.transpose(3, 1, 0, 2)
    magnetic = -angular
    result = {
        'L_over_hbar': angular,
        'mu_orb_over_muB': magnetic,
        'valid': valid.transpose(1, 0, 2),
        'bands_range': (target_first, target_last),
        'is_diagonal_orbital_approximation': True,
    }
    if return_operator:
        operator = np.zeros((3, velocity.shape[1], velocity.shape[0],
                             len(target), len(target)), dtype=magnetic.dtype)
        diagonal = np.arange(len(target))
        operator[..., diagonal, diagonal] = magnetic
        result['operator'] = operator
    return result


def compute_electronic_orbital_moment_BZ(wfdb, dipolesdb, electronsdb,
                                         bands_range=None, degen_tol=1.e-5,
                                         dmats=None):
    """Expand the diagonal electronic orbital moment to the full BZ.

    Orbital moment is an axial, time-reversal-odd operator, so the expansion
    reuses :meth:`YamboSpinDB.expand_operator_fullBZ`.
    """
    ibz = compute_electronic_orbital_moment(
        dipolesdb, electronsdb, bands_range=bands_range,
        degen_tol=degen_tol, return_operator=True)
    operator_ibz = ibz['operator']   # (3, nk_ibz, nspin, ntarget, ntarget)
    valid_ibz = ibz['valid']         # (nk_ibz, nspin, ntarget)

    lattice = wfdb.ydb
    nk_ibz, nspin, ntarget = operator_ibz.shape[1:4]
    if nk_ibz != lattice.ibz_nkpoints:
        raise ValueError("orbital-moment and wavefunction IBZ grids do not match")

    first, last = ibz['bands_range']
    helper = YamboSpinDB(nk_ibz, nspin, first, last,
                         [first, last], ntarget, None)
    operator_bz = helper.expand_operator_fullBZ(operator_ibz, wfdb, dmats=dmats)
    valid_bz = valid_ibz[np.asarray(lattice.kpoints_indexes, dtype=int)]

    return {
        'operator_BZ': operator_bz,
        'valid_BZ': valid_bz,
        'bands_range': ibz['bands_range'],
        'is_diagonal_orbital_approximation': True,
    }


def compute_electronic_orbital_moment_BZ_from_save(path='.', save='SAVE',
                                                    dipoles_filename='ndb.dipoles',
                                                    bands_range=None, degen_tol=1.e-5):
    """Load standard Yambo databases and return the full-BZ orbital moment."""
    import os

    from yambopy.dbs.electronsdb import YamboElectronsDB
    from yambopy.dbs.dipolesdb import YamboDipolesDB
    from yambopy.dbs.latticedb import YamboLatticeDB
    from yambopy.dbs.wfdb import YamboWFDB

    lattice = YamboLatticeDB.from_db_file(os.path.join(path, save, 'ns.db1'))
    electronsdb = YamboElectronsDB.from_db_file(folder=os.path.join(path, save), Expand=False)
    dipolesdb = YamboDipolesDB.from_db_file(
        lattice, filename=os.path.join(path, dipoles_filename),
        dip_type='v', expand=False, project=False)

    first, last = (bands_range if bands_range is not None
                   else dipolesdb.bands_range)
    wfdb = YamboWFDB(path=path, save=save, latdb=lattice,
                     bands_range=[first - 1, last])

    result = compute_electronic_orbital_moment_BZ(
        wfdb, dipolesdb, electronsdb, bands_range=bands_range, degen_tol=degen_tol)
    result['lattice'] = lattice
    result['wfdb'] = wfdb
    return result
