"""Read spin matrix elements from Yambo ``ndb.dipoles`` databases."""

import os
import re

import numpy as np
from netCDF4 import Dataset

from yambopy.tools.citations import citation


class YamboSpinDB:
    """Read ``DIP_spin`` from a Yambo ``ndb.dipoles`` database.

    ``spin_full`` follows the ``YamboDipolesDB.dipoles_full`` convention:
    ``spin_full[k,alpha,n,m] = <m,k|sigma_alpha|n,k>`` for one storage-spin
    channel. Values are dimensionless Pauli matrix elements; divide by two
    for ``S/hbar``. Two-channel collinear spin is not supported.
    """

    def __init__(self, nk_ibz, spin, min_band, max_band, bands_range,
                 nbands, spin_full):
        """Initialize the YamboSpinDB."""
        self.nk_ibz = nk_ibz
        self.spin = spin
        self.min_band = min_band
        self.max_band = max_band
        self.bands_range = bands_range
        self.nbands = nbands
        self.spin_full = spin_full

    @classmethod
    def from_db_file(cls, filename='ndb.dipoles', bands_range=None):
        """
        Initialize the class from a Yambo dipole database.

        Parameters
        ----------
        filename : str, optional
            Path to ``ndb.dipoles`` or to the directory containing it.
        bands_range : list, optional
            One-based inclusive band range. An empty list reads all bands.
        """
        bands_range = [] if bands_range is None else list(bands_range)
        if bands_range and (len(bands_range) != 2 or
                            any(int(b) != b for b in bands_range)):
            raise ValueError("bands_range must contain two integer band numbers")
        if os.path.isdir(filename):
            filename = os.path.join(filename, 'ndb.dipoles')
        if not os.path.isfile(filename):
            raise FileNotFoundError("error opening %s in YamboSpinDB" % filename)

        with Dataset(filename) as database:
            # General parameters
            nk_ibz = int(database.variables['HEAD_R_LATT'][2])
            if 'DIP_spin' in database.variables:
                variable = database.variables['DIP_spin']
                if (variable.ndim != 6 or variable.shape[1] != nk_ibz or
                        variable.shape[2] != variable.shape[3]):
                    raise ValueError("DIP_spin must contain square band matrices on the IBZ")
                spin = variable.shape[0]
                if spin != 1:
                    raise ValueError(
                        "two-channel collinear DIP_spin is not supported")
                min_band, max_band = cls._stored_band_range(
                    database, variable.shape[2])

                if len(bands_range) == 0:
                    bands_range = [min_band, max_band]
                first, last = map(int, bands_range)
                if not min_band <= first <= last <= max_band:
                    raise ValueError(
                        "invalid bands_range, db contains [%d,%d]" %
                        (min_band, max_band))
                nbands = last-first+1
                start, stop = first-min_band, last-min_band+1

                raw = variable[:, :, start:stop, start:stop, :, :]
                spin_full = cls._as_complex(raw)
                if spin == 1:
                    spin_full = np.squeeze(spin_full, axis=0)
                spin_full = np.swapaxes(spin_full, spin, spin+2)
            else:
                spin = 1
                min_band = (int(np.rint(database.variables['PARS'][0]).item())
                            if 'PARS' in database.variables else 1)
                spin_full, bands_range, max_band = cls._read_fragments(
                    filename, nk_ibz, min_band, bands_range)
                nbands = spin_full.shape[-1]

        return cls(nk_ibz, spin, min_band, max_band, list(bands_range),
                   nbands, spin_full)

    @staticmethod
    def _as_complex(raw):
        """Convert the final Yambo real-imaginary axis to complex values."""
        if np.ma.isMaskedArray(raw) and np.any(np.ma.getmaskarray(raw)):
            raise ValueError("DIP_spin contains missing values")
        raw = np.asarray(raw)
        if raw.shape[-2:] == (2, 3):
            raw = np.swapaxes(raw, -2, -1)
        if raw.shape[-2:] != (3, 2):
            raise ValueError("DIP_spin must end with (xyz,re_im) axes")
        values = raw[..., 0] + 1j*raw[..., 1]
        if not np.all(np.isfinite(values)):
            raise ValueError("DIP_spin contains non-finite values")
        return values

    @staticmethod
    def _stored_band_range(database, axis_size):
        if 'PARS' in database.variables:
            first, last = np.rint(database.variables['PARS'][:2]).astype(int)
            if last-first+1 == axis_size:
                return first, last
            raise ValueError("DIP_spin band axes do not match PARS")
        return 1, axis_size

    @classmethod
    def _read_fragments(cls, filename, nk_ibz, min_band, bands_range):
        fragments = sorted(
            p for p in cls._fragment_files(filename))
        if not fragments:
            raise ValueError(
                'DIP_spin is absent; run Yambo with DipComputed="R P V Spin"')

        pattern = re.compile(r'DIP_spin_k_(\d{4})_spin_(\d{4})$')
        records = {}
        selected = None
        stored_bands = None
        for fragment in fragments:
            with Dataset(fragment) as database:
                for name, variable in database.variables.items():
                    match = pattern.fullmatch(name)
                    if match is None:
                        continue
                    ik, spin = map(int, match.groups())
                    if spin != 1:
                        raise ValueError(
                            "two-channel collinear DIP_spin is not supported")
                    if variable.ndim != 4 or variable.shape[0] != variable.shape[1]:
                        raise ValueError("fragmented DIP_spin must contain square band matrices")
                    if stored_bands is not None and stored_bands != variable.shape[0]:
                        raise ValueError("DIP_spin fragments have different band windows")
                    stored_bands = variable.shape[0]
                    selected = ([min_band, min_band+stored_bands-1]
                                if len(bands_range) == 0
                                else list(map(int, bands_range)))
                    first, last = selected
                    start, stop = first-min_band, last-min_band+1
                    if start < 0 or stop > stored_bands or first > last:
                        raise ValueError(
                            "bands_range lies outside fragmented DIP_spin")
                    if ik in records or not 1 <= ik <= nk_ibz:
                        raise ValueError("duplicate or out-of-range DIP_spin k point")
                    records[ik] = cls._as_complex(
                        variable[start:stop, start:stop])

        missing = sorted(set(range(1, nk_ibz+1)).difference(records))
        if missing:
            raise ValueError("fragmented DIP_spin is missing k points: %s" % missing)
        spin_full = np.stack([records[ik] for ik in range(1, nk_ibz+1)])
        return np.swapaxes(spin_full, 1, 3), selected, min_band+stored_bands-1

    @staticmethod
    def _fragment_files(filename):
        folder, name = os.path.split(filename)
        prefix = name + '_fragment_'
        return (os.path.join(folder, entry) for entry in os.listdir(folder or '.')
                if entry.startswith(prefix))

    def as_operator(self):
        """Return ``<m,k|sigma_alpha|n,k>`` in generic Yambopy layout."""
        spin = self.spin_full
        if self.spin == 1:
            spin = spin[None, ...]
        # spin_full follows dipoles_full, whose two band axes are interchanged
        return spin.transpose(2, 1, 0, 4, 3)

    def expand_fullBZ(self, wfdb, dmats=None):
        """Expand the spin operator from the IBZ to the full BZ.

        ``wfdb`` supplies the BZ-to-IBZ map, symmetry operations and the
        band-gauge matrices computed by :meth:`YamboWFDB.Dmat`.  The result
        has shape ``(3,nk_bz,nspin,nbands,nbands)``.
        """
        return self.expand_operator_fullBZ(self.as_operator(), wfdb, dmats)

    def expand_operator_fullBZ(self, operator_ibz, wfdb, dmats=None):
        """Expand an axial, time-reversal-odd IBZ operator with ``Dmat``.

        The operator is first rotated in band gauge, ``D O D^dagger``, and
        then as an axial vector, ``det(R) R``. Antiunitary images are
        conjugated and sign-flipped.  The result has layout
        ``(3,nk_bz,nspin,nbands,nbands)``.
        """
        lattice = wfdb.ydb
        if self.spin != 1:
            raise NotImplementedError("full-BZ spin-channel mapping is not implemented")
        if self.nk_ibz != lattice.ibz_nkpoints:
            raise ValueError("operator and wavefunctions use different IBZ grids")

        band_slice = self._band_slice_in(wfdb)
        dmats = self._checked_dmats(wfdb, dmats, band_slice)

        ibz_source = self._ibz_source_bz_indices(lattice)
        operator_ibz = np.asarray(operator_ibz)
        expected_operator = (3, self.nk_ibz, self.spin,
                             self.nbands, self.nbands)
        if operator_ibz.shape != expected_operator:
            raise ValueError("unexpected IBZ operator shape")

        nsym = len(lattice.sym_car)
        antiunitary_start = nsym/(1+int(np.rint(lattice.time_rev)))

        result = np.empty((3, wfdb.nkBZ, self.spin,
                           self.nbands, self.nbands), complex)
        for ik_bz, (ik_ibz, isym) in enumerate(zip(
                lattice.kpoints_indexes, lattice.symmetry_indexes)):
            is_antiunitary = isym >= antiunitary_start
            dmat = dmats[isym, ibz_source[ik_ibz]]
            if not np.allclose(dmat @ dmat.conj().transpose(0, 2, 1),
                               np.eye(self.nbands), atol=1.e-3, rtol=0.):
                raise ValueError("Dmat is not unitary in the selected band window; "
                                 "include the complete symmetry-connected subspace")
            result[:, ik_bz] = self._rotate_operator(
                operator_ibz[:, ik_ibz], dmat,
                lattice.sym_car[isym], is_antiunitary)

        return result

    @staticmethod
    def _rotate_operator(op_ibz, dmat, rotation, is_antiunitary):
        """Transform one Cartesian operator triplet from k to R*k.

        ``op_ibz`` has shape ``(3,nspin,nbands,nbands)``; ``dmat`` has shape
        ``(nspin,nbands,nbands)``; ``rotation`` is the 3x3 Cartesian
        symmetry matrix ``R``.
        """
        if is_antiunitary:
            op_ibz = op_ibz.conj()
        dmat_dagger = dmat.conj().transpose(0, 2, 1)
        gauge_rotated = np.einsum(
            'sab,vsbc,scd->vsad', dmat, op_ibz, dmat_dagger, optimize=True)
        cartesian_rotation = np.linalg.det(rotation)*rotation
        if is_antiunitary:
            cartesian_rotation = -cartesian_rotation
        return np.einsum(
            'uv,vsab->usab', cartesian_rotation, gauge_rotated, optimize=True)

    def _band_slice_in(self, wfdb):
        """Where the DIP_spin band window sits inside wfdb's band window."""
        first = self.bands_range[0]-1
        start = first-wfdb.min_bnd
        stop = start+self.nbands
        if start < 0 or stop > wfdb.nbands:
            raise ValueError("wavefunctions do not contain the DIP_spin bands")
        return slice(start, stop)

    def _checked_dmats(self, wfdb, dmats, band_slice):
        if dmats is None:
            dmats = wfdb.Dmat()
        dmats = np.asarray(dmats)
        expected = (len(wfdb.ydb.sym_car), wfdb.nkBZ, self.spin,
                    wfdb.nbands, wfdb.nbands)
        if dmats.shape != expected:
            raise ValueError("Dmat has incompatible symmetry, k-point, spin or band axes")
        return dmats[:, :, :, band_slice, band_slice]

    @staticmethod
    def _ibz_source_bz_indices(lattice):
        """Full-BZ index of each identity-mapped IBZ representative."""
        identity = np.argmin(np.linalg.norm(
            lattice.sym_car-np.eye(3)[None, :, :], axis=(1, 2)))
        nk_ibz = lattice.ibz_nkpoints
        source = np.empty(nk_ibz, dtype=int)
        for ik in range(nk_ibz):
            matches = np.flatnonzero(
                (lattice.kpoints_indexes == ik)
                & (lattice.symmetry_indexes == identity))
            if len(matches) != 1:
                raise ValueError("could not identify an identity-mapped IBZ point")
            source[ik] = matches[0]
        return source


def read_spin_m_e_iBZ(filename='ndb.dipoles', bands_range=None):
    """Read ``DIP_spin`` and return the generic Yambopy operator layout."""
    selected = [] if bands_range is None else bands_range
    return YamboSpinDB.from_db_file(filename, selected).as_operator()


@citation("A. R. Kshirsagar et al. Phys. Rev. B 112, 12 (2025)")
def compute_exciton_spin_from_dipoles(lattice, excdb, wfdb, spin_db, dmats=None,
                                       contribution='b', diagonal=False):
    """Compute exciton ``S_z`` using ``DIP_spin`` as electronic input.

    Parameters
    ----------
    lattice : YamboLatticeDB
        Lattice database (as used to build ``excdb`` and ``wfdb``).
    excdb : YamboExcitonDB
        BSE exciton database (``Load_WF=True``).
    wfdb : YamboWFDB
        Wavefunction database restricted to the same band window as
        ``spin_db`` and ``excdb.table``. Only its BZ map/``Dmat`` machinery
        is used here. Computing ``Dmat`` requires the wavefunctions unless
        the matrices have already been supplied.
    spin_db : YamboSpinDB
        The ``DIP_spin`` database, already restricted to the same one-based
        inclusive band range as ``wfdb``/``excdb``.
    dmats : ndarray, optional
        Precomputed :meth:`YamboWFDB.Dmat` array.
    contribution : str, optional
        'b' (both electron and hole, default), 'e', or 'h'.
    diagonal : bool, optional
        If True, only diagonal exciton spin elements are computed.

    Returns
    -------
    exe_Sz : ndarray
        Spin matrix elements for excitons.
    """
    from yambopy.bse.exciton_matrix_elements import exciton_X_matelem

    assert wfdb.nspinor == 2, "Makes sense only for nspinor = 2"
    assert np.min(excdb.table[:, 1]) - 1 == wfdb.min_bnd, \
        "wfdb and exciton db are inconsistant (Bands)"
    assert np.max(excdb.table[:, 2]) == wfdb.min_bnd + wfdb.nbands, \
        "wfdb and exciton db are inconsistant (Bands)"

    first = spin_db.bands_range[0] - 1
    assert first == wfdb.min_bnd and spin_db.nbands == wfdb.nbands, \
        "spin_db and wfdb address different bands"

    spin_bz = spin_db.expand_fullBZ(wfdb, dmats=dmats)   # (3,nkBZ,nspin,nb,nb)
    elec_sz = 0.5 * spin_bz[2, :, 0, :, :]                # sigma_z -> S_z
    assert elec_sz.shape == (wfdb.nkBZ, wfdb.nbands, wfdb.nbands)

    Akcv = excdb.get_Akcv()
    excQpt = lattice.lat @ excdb.car_qpoint

    exe_Sz = exciton_X_matelem(excQpt, np.array([0, 0, 0]), Akcv, Akcv,
                               elec_sz[None, :, None, ...], wfdb.kBZ,
                               diagonal_only=diagonal, contribution=contribution)
    return exe_Sz[0]
