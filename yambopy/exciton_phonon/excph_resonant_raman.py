##
# Authors: PMI, MN
##

import numpy as np
from tqdm import tqdm
from numba import njit, prange
from yambopy.units import ha2ev


#@citation("PHYSICAL REVIEW B 113, 085201 (2026)")
def ip_resonant_raman_tensor_oneph(laser_energies, ph_energies, el_energies,
                                   elec_dipoles, eph_g, cell_vol,
                                   broad=0.1, ph_freq_threshold=5.0):
    """
    ! 1-phonon Raman tensor at independent-particle (IP) level
    ! Phys. Rev. B 113, 085201 – Published 2 February, 2026
    ! DOI: 10.1103/ty8m-mgml
    Implements Eq. (D2) of the paper above.

    Mirrors `compute_Raman_oneph_ip` from
    https://github.com/muralidhar-nalabothula/PhdScripts/blob/main/exph/raman.py

    All band-window bookkeeping (which bands, IBZ -> BZ expansion, etc.)
    is the caller's responsibility. The function assumes that every input
    is already sliced and aligned onto the same nb-band window in the
    order (valence, then conduction) along the band axis. 

    Parameters
    ----------
    laser_energies : (nfreqs,) float ndarray
        Incoming laser energies in eV.
    ph_energies : (nmodes,) float ndarray
        Phonon energies at q=0 in eV.
    el_energies : (nk, nb) float ndarray
        Single-particle energies in eV, valence-first then conduction.
    elec_dipoles : (3, nk, nc, nv) complex ndarray
        Velocity-gauge dipoles <ck|v|vk> (yambo convention, a.u.).
        nv + nc must equal nb.
    eph_g : (nmodes, nk, nb, nb) complex ndarray
        Electron-phonon matrix elements at q=0 in Hartree, with index
        order (mode, k, final_band, initial_band) and the standard
        1/sqrt(2*omega_nu) normalization.
    cell_vol : float
        Unit-cell volume in bohr^3.
    broad : float, optional
        Lorentzian broadening Gamma in eV. Default: 0.1.
    ph_freq_threshold : float, optional
        Acoustic-mode cutoff in cm^-1. Default: 5.

    Returns
    -------
    raman_tensor : (nfreqs, nmodes, 3, 3) complex ndarray
    """
    cm1_to_Ha = 0.12398e-3 / ha2ev

    nk, nb     = el_energies.shape
    _, _, nc, nv = elec_dipoles.shape
    nmodes     = ph_energies.shape[0]
    if nv + nc != nb:
        raise ValueError("nv + nc != nb (%d + %d != %d)" % (nv, nc, nb))
    if elec_dipoles.shape != (3, nk, nc, nv):
        raise ValueError("elec_dipoles shape %s != (3, %d, %d, %d)"
                         % (elec_dipoles.shape, nk, nc, nv))
    if eph_g.shape != (nmodes, nk, nb, nb):
        raise ValueError("eph_g shape %s != (%d, %d, %d, %d)"
                         % (eph_g.shape, nmodes, nk, nb, nb))

    # Unit conversion (eV -> Ha)
    laser_Ha     = np.asarray(laser_energies) / ha2ev
    ph_Ha        = np.asarray(ph_energies)    / ha2ev
    El_Ha        = np.asarray(el_energies)    / ha2ev
    broad_Ha     = (broad / ha2ev) / 2.0
    ph_thresh_Ha = ph_freq_threshold * cm1_to_Ha

    # Vertical c-v gaps minus i*Gamma
    delta = (El_Ha[:, nv:, None] - El_Ha[:, None, :nv]) - 1j * broad_Ha   # (nk, nc, nv)

    d_abs = elec_dipoles.conj()    # incoming-photon vertex
    d_emi = elec_dipoles           # outgoing-photon vertex

    gcc = eph_g[:, :, nv:, nv:]    # (nmodes, nk, nc, nc)
    gvv = eph_g[:, :, :nv, :nv]    # (nmodes, nk, nv, nv)

    nfreqs  = laser_Ha.size
    Ram     = np.zeros((nfreqs, nmodes, 3, 3), dtype=np.complex128)
    ram_fac = 1.0 / nk / np.sqrt(cell_vol)

    for w_idx in tqdm(range(nfreqs), desc="IP one-phonon Raman"):
        wL = laser_Ha[w_idx]

        dipS_res  = d_abs / (wL - delta)[None, ...]
        dipS_ares = d_emi / (wL + delta)[None, ...]

        for m in range(nmodes):
            wph = ph_Ha[m]
            if abs(wph) <= ph_thresh_Ha:
                continue

            dipSp_res  = d_emi / (wL - delta - wph)[None, ...]
            dipSp_ares = d_abs / (wL + delta - wph)[None, ...]

            gcc_m = gcc[m][None, ...]
            gvv_m = gvv[m][None, ...]

            tmp_res  = gcc_m.conj() @ dipSp_res  - dipSp_res  @ gvv_m.conj()
            tmp_ares = gcc_m        @ dipSp_ares - dipSp_ares @ gvv_m

            Ram[w_idx, m] = (np.einsum('akcv,bkcv->ab', dipS_res,  tmp_res)
                           + np.einsum('akcv,bkcv->ab', dipS_ares, tmp_ares))
            Ram[w_idx, m] *= np.sqrt(abs(wL - wph) / wL) * ram_fac

    return Ram


@njit(cache=True, nogil=True, parallel=True)
def _exc_raman_oneph_kernel(laser_Ha, ph_Ha, BS_energies, exc_dip_absorp,
                            exc_ph, ram_fac, ph_thresh_Ha):
    """Numba kernel for `exc_resonant_raman_oneph`.

    Mirrors `compute_Raman_oneph_exc_numba` from Murali's PhdScripts.
    All inputs in atomic units.

    Parameters
    ----------
    laser_Ha       : (nfreqs,) float    – laser energies (Ha)
    ph_Ha          : (nmodes,) float    – phonon energies (Ha)
    BS_energies    : (nexc,)  complex   – exciton energies minus i*Gamma/2 (Ha)
    exc_dip_absorp : (npol, nexc) complex  – absorption-form excitonic dipoles
                                             (= conj of yambopy emission output)
    exc_ph         : (nmodes, nexc_in, nexc_out) complex
                     yambopy convention (phonon absorption, [mode, init, fin]).
                     The numba code conjugates internally to obtain the Stokes
                     (phonon emission) branch.
    ram_fac        : float              – 1 / N_k / sqrt(cell_vol)
    ph_thresh_Ha   : float              – acoustic-mode cutoff (Ha)
    """
    nfreqs = len(laser_Ha)
    nmodes, nexc, _ = exc_ph.shape
    npol = exc_dip_absorp.shape[0]
    Ram = np.zeros((nfreqs, nmodes, 3, 3), dtype=exc_dip_absorp.dtype)

    dipS_ares_base = np.conj(exc_dip_absorp)   # = emission-form dipole
    dipSp_res_base = np.conj(exc_dip_absorp)

    for i_ome in prange(nfreqs):
        wL = laser_Ha[i_ome]

        dipS_res      = exc_dip_absorp / (wL - BS_energies)
        dipS_res_conj = np.conj(dipS_res)
        dipS_ares     = dipS_ares_base / (wL + BS_energies)

        for i in range(nmodes):
            freq = ph_Ha[i]
            if np.abs(freq) <= ph_thresh_Ha:
                continue

            dipSp_res  = dipSp_res_base / (wL - BS_energies - freq)
            dipSp_ares = exc_dip_absorp / (wL + BS_energies - freq)

            ex_ph_T = exc_ph[i].T

            # Resonant term (Stokes): conj of full product gives g_emi
            dipSp_res_T_conj = np.conj(dipSp_res.T)
            term1 = dipS_res_conj @ ex_ph_T @ dipSp_res_T_conj
            term1 = np.conj(term1)

            # Anti-resonant term
            dipSp_ares_T = dipSp_ares.T
            term2 = dipS_ares @ ex_ph_T @ dipSp_ares_T

            scale = np.sqrt(np.abs(wL - freq) / wL) * ram_fac
            Ram[i_ome, i, :npol, :npol] = (term1 + term2) * scale

    return Ram


def exc_resonant_raman_oneph(laser_energies, ph_energies, exc_energies,
                             exc_dipoles, exc_ph_mat_el, n_kpts, cell_vol,
                             broad=0.1, ph_freq_threshold=5.0,
                             precision='d'):
    """
    1-phonon Raman tensor including excitonic effects (Stokes).

    Numba-accelerated. Mirrors `compute_Raman_oneph_exc` from
    https://github.com/muralidhar-nalabothula/PhdScripts/blob/main/exph/raman.py
    See also Sven Reichardt et al., Sci. Adv. 6, eabb5915 (2020).

    All band/exciton/k-point bookkeeping is the caller's responsibility.
    The function assumes the inputs are already sliced consistently
    (same nexc on every exciton-indexed axis) and that
    `exc_dipoles` is in yambo's emission convention,
    `exc_ph_mat_el` is in yambopy's raw (phonon-absorption) convention
    with shape (nmodes, nexc_in, nexc_out). The kernel handles the
    conjugation required to convert to phonon-emission (Stokes) Raman.

    Parameters
    ----------
    laser_energies : (nfreqs,) float ndarray
        Incoming laser energies in eV.
    ph_energies : (nmodes,) float ndarray
        Phonon energies at q=0 in eV.
    exc_energies : (nexc,) float ndarray
        Exciton energies at Q=0 in eV.
    exc_dipoles : (3, nexc) complex ndarray
        Velocity-gauge exciton dipoles in emission convention
        (a.u.), as returned by yambopy `exc_dipoles_pol`.
    exc_ph_mat_el : (nmodes, nexc, nexc) complex ndarray
        Exciton-phonon matrix elements at q=0 in Hartree, in yambopy's
        raw convention (phonon absorption, [mode, initial, final]).
        This is `exciton_phonon_matelem(...)[0]` (q=0 slice) after
        slicing nexc_in == nexc_out == nexc.
    n_kpts : int
        Number of k-points in the BZ (for the 1/N_k prefactor).
    cell_vol : float
        Unit-cell volume in bohr^3.
    broad : float, optional
        Lorentzian broadening Gamma in eV. Default: 0.1.
    ph_freq_threshold : float, optional
        Acoustic-mode cutoff in cm^-1. Default: 5.
    precision : {'d', 's'}, optional
        Floating-point precision used inside the numba kernel.
        'd' = double / complex128 (default), 's' = single / complex64.

    Returns
    -------
    raman_tensor : (nfreqs, nmodes, 3, 3) complex ndarray
    """
    cm1_to_Ha = 0.12398e-3 / ha2ev

    nfreqs = np.atleast_1d(laser_energies).size
    nmodes = ph_energies.shape[0]
    nexc   = exc_energies.shape[0]
    if exc_dipoles.shape[0] != 3 or exc_dipoles.shape[1] != nexc:
        raise ValueError("exc_dipoles shape %s != (3, %d)"
                         % (tuple(exc_dipoles.shape), nexc))
    if exc_ph_mat_el.shape != (nmodes, nexc, nexc):
        raise ValueError("exc_ph_mat_el shape %s != (%d, %d, %d)"
                         % (tuple(exc_ph_mat_el.shape), nmodes, nexc, nexc))
    if cell_vol <= 0.:
        raise ValueError("cell_vol must be positive")

    prec = str(precision).strip().lower()
    if prec in ('d', 'double'):
        f_type, c_type = np.float64, np.complex128
    elif prec in ('s', 'single'):
        f_type, c_type = np.float32, np.complex64
    else:
        raise ValueError("precision must be 'd' or 's'")

    # eV -> Ha
    laser_Ha     = np.atleast_1d(np.asarray(laser_energies)) / ha2ev
    ph_Ha        = np.asarray(ph_energies) / ha2ev
    broad_Ha     = (broad / ha2ev) / 2.0
    ph_thresh_Ha = ph_freq_threshold * cm1_to_Ha
    ram_fac      = 1.0 / float(n_kpts) / np.sqrt(cell_vol)

    BS_energies  = np.asarray(exc_energies) / ha2ev - 1j * broad_Ha
    dip_absorp   = np.conj(exc_dipoles)               # yambopy emission -> absorption

    # Contiguous arrays in the chosen precision (required by the @njit kernel).
    laser_c      = np.ascontiguousarray(laser_Ha,    dtype=f_type)
    ph_c         = np.ascontiguousarray(ph_Ha,       dtype=f_type)
    BS_c         = np.ascontiguousarray(BS_energies, dtype=c_type)
    dip_absorp_c = np.ascontiguousarray(dip_absorp,  dtype=c_type)
    exc_ph_c     = np.ascontiguousarray(exc_ph_mat_el, dtype=c_type)

    return _exc_raman_oneph_kernel(laser_c, ph_c, BS_c, dip_absorp_c,
                                   exc_ph_c, ram_fac, ph_thresh_Ha)


# 1 eV in cm^-1
_EV_TO_CM1 = 8065.54


def ip_raman_spectrum_oneph(laser_energy, ph_energies,
                            el_energies=None, elec_dipoles=None, eph_g=None,
                            cell_vol=None,
                            raman_tensor=None,
                            broad=0.1, ph_freq_threshold=5.0,
                            energy_grid=None, energy_units='cm-1',
                            spectrum_broad=0.1, broad_type='lorentzian',
                            pol_in=None, pol_out=None):
    """
    1-phonon IP Raman spectrum at a single laser energy.

    Implements Eq. (D1) of Phys. Rev. B 113, 085201 (2026):

        dσ/dΩ  ∝  Σ_μν |M^{μν}(ω_L, ω_λ)|² · δ(E − ℏω_λ)

    where M^{μν}(ω_L, ω_λ) is the per-mode 3×3 IP Raman tensor (D2),
    computed by `ip_resonant_raman_tensor_oneph`. The δ-function is
    replaced by a Lorentzian (default) or Gaussian of FWHM
    `spectrum_broad`.

    If a `raman_tensor` is given, it is used directly and no IP-Raman
    tensor calculation is performed; the function then only needs
    `laser_energy` (to set ω_L for the prefactor *already inside* the
    tensor — i.e., we just trust it) and `ph_energies` (to know where
    to put each peak).

    Parameters
    ----------
    laser_energy : float
        Single laser energy ω_L in eV.
    ph_energies : (nmodes,) float ndarray
        Phonon energies at q=0 in eV (used to place each peak on the
        spectrum and, if ``raman_tensor is None``, passed through to
        the tensor function).
    el_energies, elec_dipoles, eph_g, cell_vol : optional
        Required only when ``raman_tensor is None``. Same meaning as in
        `ip_resonant_raman_tensor_oneph`.
    raman_tensor : (nmodes, 3, 3) or (1, nmodes, 3, 3) complex ndarray, optional
        Pre-computed Raman tensor (as returned by
        `ip_resonant_raman_tensor_oneph` at a single laser energy).
        If given, the tensor is reused and the IP-Raman calculation is
        skipped.
    broad : float, optional
        Lorentzian broadening (eV) passed to the tensor function
        (electronic Γ). Default: 0.1.
    ph_freq_threshold : float, optional
        Acoustic-mode cutoff in cm^-1 passed to the tensor function.
        Default: 5. NOTE: acoustic modes are NOT removed from the
        output spectrum — they appear with whatever (typically zero or
        noisy) intensity comes out of the tensor calculation, so the
        user can inspect them.
    energy_grid : 1D ndarray, optional
        Energy grid on which the spectrum is evaluated, in
        `energy_units`. If None, an evenly-spaced grid is built
        covering the phonon energy range with a step ≈
        `spectrum_broad/4`.
    energy_units : {'cm-1', 'eV'}, optional
        Units of the output `energy_grid` and of `spectrum_broad`.
        Default: 'cm-1'.
    spectrum_broad : float, optional
        FWHM of the broadening function applied to each phonon
        δ-peak, in `energy_units`. Default: 0.1.
    broad_type : {'gaussian', 'lorentzian'}, optional
        Shape of the broadening function. Default: 'gaussian'.
    pol_in, pol_out : array-like of shape (3,), optional
        Incoming / outgoing photon polarisation 3-vectors (Cartesian,
        possibly complex). If BOTH are given, the per-mode intensity
        is |Σ_{μν} e_out[μ] R[λ,μ,ν] e_in[ν]|² (vectors normalised
        internally). If either is None (default), the unpolarised D1
        sum is used: I_λ = Σ_{μν} |R[λ,μ,ν]|².

    Returns
    -------
    energy_grid : (nE,) float ndarray
        Energy axis in `energy_units`.
    intensity   : (nE,) float ndarray
        Raman intensity at each energy.
    """
    # ---- Compute or accept the Raman tensor -------------------------------
    nmodes = int(ph_energies.shape[0])
    if raman_tensor is None:
        for name, arr in (('el_energies', el_energies),
                          ('elec_dipoles', elec_dipoles),
                          ('eph_g', eph_g),
                          ('cell_vol', cell_vol)):
            if arr is None:
                raise ValueError("`%s` must be provided when `raman_tensor` is None" % name)
        R = ip_resonant_raman_tensor_oneph(
                np.array([float(laser_energy)]),
                ph_energies, el_energies, elec_dipoles, eph_g, cell_vol,
                broad=broad, ph_freq_threshold=ph_freq_threshold)
        R = R[0]                                          # (nmodes, 3, 3)
    else:
        R = np.asarray(raman_tensor)
        if R.ndim == 4 and R.shape[0] == 1:
            R = R[0]
        if R.shape != (nmodes, 3, 3):
            raise ValueError("raman_tensor shape %s != (%d, 3, 3)" %
                             (tuple(R.shape), nmodes))

    # ---- Per-mode intensity ----------------------------------------------
    if (pol_in is None) or (pol_out is None):
        # Eq. (D1) sum over polarisation indices
        I_mode = np.sum(np.abs(R)**2, axis=(1, 2))        # (nmodes,)
    else:
        e_in  = np.asarray(pol_in,  dtype=complex)
        e_out = np.asarray(pol_out, dtype=complex)
        e_in  = e_in  / np.linalg.norm(e_in)
        e_out = e_out / np.linalg.norm(e_out)
        amp   = np.einsum('m, lmn, n -> l', e_out, R, e_in)
        I_mode = np.abs(amp)**2                            # (nmodes,)

    # ---- Energy axis and broadening unit conversions ---------------------
    units = str(energy_units).strip().lower()
    if units in ('cm-1', 'cm^-1', 'cm', 'wavenumber'):
        ph_axis  = np.asarray(ph_energies) * _EV_TO_CM1
        unit_lbl = 'cm-1'
    elif units in ('ev',):
        ph_axis  = np.asarray(ph_energies, dtype=float)
        unit_lbl = 'eV'
    else:
        raise ValueError("energy_units must be 'cm-1' or 'eV'")

    fwhm = float(spectrum_broad)
    if fwhm <= 0.:
        raise ValueError("spectrum_broad (FWHM) must be > 0")

    if energy_grid is None:
        emin = min(0.0, float(np.min(ph_axis)) - 5.0 * fwhm)
        emax = float(np.max(ph_axis)) + 5.0 * fwhm
        step = max(fwhm / 4.0, (emax - emin) / 5000.0)
        energy_grid = np.arange(emin, emax + step, step)
    energy_grid = np.asarray(energy_grid, dtype=float)

    # ---- Place each mode's intensity on the energy grid ------------------
    diff = energy_grid[:, None] - ph_axis[None, :]         # (nE, nmodes)
    btype = str(broad_type).strip().lower()
    if btype == 'gaussian':
        sigma   = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        profile = np.exp(-0.5 * (diff / sigma)**2) / (sigma * np.sqrt(2.0 * np.pi))
    elif btype == 'lorentzian':
        hwhm    = 0.5 * fwhm
        profile = (hwhm / np.pi) / (diff**2 + hwhm**2)
    else:
        raise ValueError("broad_type must be 'gaussian' or 'lorentzian'")

    intensity = np.einsum('em, m -> e', profile, I_mode)

    return energy_grid, intensity
