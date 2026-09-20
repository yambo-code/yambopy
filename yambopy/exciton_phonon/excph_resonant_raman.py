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
    """Numba kernel for `exc_resonant_raman_tensor_oneph`.

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


def exc_resonant_raman_tensor_oneph(laser_energies, ph_energies, exc_energies,
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
                            pol_in=None, pol_out=None,
                            n_probe=1):
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
    n_probe : int, optional
        Number of laser-energy probe points. Default 1 (single laser
        energy ω_L, identical to the previous behaviour). For
        ``n_probe > 1`` the per-mode intensity is computed at
        ``n_probe`` laser energies centred on ω_L and spaced by
        ``broad/4``, then integrated over that window (Simpson). Use an
        odd ``n_probe``. Ignored when ``raman_tensor`` is supplied.

    Returns
    -------
    energy_grid : (nE,) float ndarray
        Energy axis in `energy_units`.
    intensity   : (nE,) float ndarray
        Raman intensity at each energy.
    """
    # ---- Build laser window (n_probe pts, centred on ω_L, spaced broad/4) --
    if raman_tensor is None:
        for name, arr in (('el_energies', el_energies), ('elec_dipoles', elec_dipoles),
                          ('eph_g', eph_g), ('cell_vol', cell_vol)):
            if arr is None:
                raise ValueError("`%s` must be provided when `raman_tensor` is None" % name)
        offsets   = (np.arange(n_probe) - n_probe // 2) * (broad / 4.0)
        laser_arr = float(laser_energy) + offsets
        R = ip_resonant_raman_tensor_oneph(
                laser_arr, ph_energies, el_energies, elec_dipoles, eph_g, cell_vol,
                broad=broad, ph_freq_threshold=ph_freq_threshold)   # (P, nmodes, 3, 3)
    else:
        R = np.asarray(raman_tensor)
        if R.ndim == 3:
            R = R[None]                                             # (1, nmodes, 3, 3)
        laser_arr = np.array([float(laser_energy)])

    # ---- Per-mode intensity at each probe point ---------------------------
    if (pol_in is None) or (pol_out is None):
        I_probe = np.sum(np.abs(R)**2, axis=(2, 3))                 # (P, nmodes)
    else:
        e_in  = np.asarray(pol_in,  dtype=complex); e_in  /= np.linalg.norm(e_in)
        e_out = np.asarray(pol_out, dtype=complex); e_out /= np.linalg.norm(e_out)
        amp   = np.einsum('m, plmn, n -> pl', e_out, R, e_in)
        I_probe = np.abs(amp)**2                                    # (P, nmodes)

    # ---- Collapse window -> one intensity per mode (Simpson integral) -----
    if I_probe.shape[0] == 1:
        I_mode = I_probe[0]                                         # (nmodes,)
    else:
        from scipy.integrate import simpson
        I_mode = simpson(I_probe, laser_arr, axis=0)                # (nmodes,)

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


def exc_raman_spectrum_oneph(laser_energy, ph_energies,
                             exc_energies=None, exc_dipoles=None,
                             exc_ph_mat_el=None, n_kpts=None, cell_vol=None,
                             raman_tensor=None,
                             broad=0.1, ph_freq_threshold=5.0,
                             energy_grid=None, energy_units='cm-1',
                             spectrum_broad=0.1, broad_type='lorentzian',
                             pol_in=None, pol_out=None,
                             precision='d'):
    """
    1-phonon excitonic Raman spectrum at a single laser energy (Stokes).

    Excitonic analogue of `ip_raman_spectrum_oneph`. Implements the
    "place each phonon peak on the energy axis with a broadening" step
    on top of the excitonic Raman tensor returned by
    `exc_resonant_raman_tensor_oneph`:

        I(E) = Σ_λ |M^{μν}(ω_L, ω_λ)|² · δ(E − ℏω_λ)

    The δ-function is replaced by a Lorentzian (default) or Gaussian of
    FWHM `spectrum_broad`.

    If `raman_tensor` is given, it is used directly and the excitonic
    tensor calculation is skipped.

    Parameters
    ----------
    laser_energy : float
        Single laser energy ω_L in eV.
    ph_energies : (nmodes,) float ndarray
        Phonon energies at q=0 in eV (used to place each peak and,
        if ``raman_tensor is None``, forwarded to the tensor function).
    exc_energies, exc_dipoles, exc_ph_mat_el, n_kpts, cell_vol : optional
        Required only when ``raman_tensor is None``. Same meaning as in
        `exc_resonant_raman_tensor_oneph`.
    raman_tensor : (nmodes, 3, 3) or (1, nmodes, 3, 3) complex ndarray, optional
        Pre-computed excitonic Raman tensor at the single laser energy
        of interest (as returned by `exc_resonant_raman_tensor_oneph`).
    broad : float, optional
        Lorentzian broadening (eV) for the electronic / excitonic Γ.
        Default: 0.1.
    ph_freq_threshold : float, optional
        Acoustic-mode cutoff in cm^-1 passed to the tensor function.
        Acoustic modes are KEPT in the output spectrum so the user can
        inspect their (typically tiny / noise-level) contribution.
        Default: 5.
    energy_grid : 1D ndarray, optional
        Energy grid on which the spectrum is evaluated, in
        `energy_units`. Auto-built around the phonon range when None.
    energy_units : {'cm-1', 'eV'}, optional
        Units of the output `energy_grid` and of `spectrum_broad`.
        Default: 'cm-1'.
    spectrum_broad : float, optional
        FWHM of the broadening function applied to each phonon peak,
        in `energy_units`. Default: 0.1.
    broad_type : {'lorentzian', 'gaussian'}, optional
        Shape of the broadening function. Default: 'lorentzian'.
    pol_in, pol_out : array-like of shape (3,), optional
        Incoming / outgoing photon polarisation Cartesian 3-vectors
        (complex allowed). If both are given, the per-mode intensity is
        ``|Σ_{μν} e_out[μ] R[λ,μ,ν] e_in[ν]|²`` (vectors normalised
        internally). If either is None (default), the unpolarised D1
        sum is used: ``I_λ = Σ_{μν} |R[λ,μ,ν]|²``.
    precision : {'d', 's'}, optional
        Forwarded to `exc_resonant_raman_tensor_oneph` when the tensor
        is computed inside this function. Default: 'd'.

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
        for name, arr in (('exc_energies',  exc_energies),
                          ('exc_dipoles',   exc_dipoles),
                          ('exc_ph_mat_el', exc_ph_mat_el),
                          ('n_kpts',        n_kpts),
                          ('cell_vol',      cell_vol)):
            if arr is None:
                raise ValueError("`%s` must be provided when `raman_tensor` is None"
                                 % name)
        R = exc_resonant_raman_tensor_oneph(
                np.array([float(laser_energy)]),
                ph_energies, exc_energies, exc_dipoles, exc_ph_mat_el,
                n_kpts, cell_vol,
                broad=broad, ph_freq_threshold=ph_freq_threshold,
                precision=precision)
        R = R[0]                                              # (nmodes, 3, 3)
    else:
        R = np.asarray(raman_tensor)
        if R.ndim == 4 and R.shape[0] == 1:
            R = R[0]
        if R.shape != (nmodes, 3, 3):
            raise ValueError("raman_tensor shape %s != (%d, 3, 3)" %
                             (tuple(R.shape), nmodes))

    # ---- Per-mode intensity ----------------------------------------------
    if (pol_in is None) or (pol_out is None):
        I_mode = np.sum(np.abs(R)**2, axis=(1, 2))            # (nmodes,)
    else:
        e_in  = np.asarray(pol_in,  dtype=complex)
        e_out = np.asarray(pol_out, dtype=complex)
        e_in  = e_in  / np.linalg.norm(e_in)
        e_out = e_out / np.linalg.norm(e_out)
        amp   = np.einsum('m, lmn, n -> l', e_out, R, e_in)
        I_mode = np.abs(amp)**2                                # (nmodes,)

    # ---- Energy axis and broadening unit conversions ---------------------
    units = str(energy_units).strip().lower()
    if units in ('cm-1', 'cm^-1', 'cm', 'wavenumber'):
        ph_axis = np.asarray(ph_energies) * _EV_TO_CM1
    elif units in ('ev',):
        ph_axis = np.asarray(ph_energies, dtype=float)
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
    diff  = energy_grid[:, None] - ph_axis[None, :]            # (nE, nmodes)
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


# =============================================================================
# Diagnostic helpers — inspect every intermediate of the excitonic Raman
# tensor at a SINGLE laser energy. Pure numpy (no numba) so the code is easy
# to read; use for debugging / spectrum inspection.
# =============================================================================

def exc_raman_components_oneph(laser_energy, ph_energies, exc_energies,
                               exc_dipoles, exc_ph_mat_el, n_kpts, cell_vol,
                               broad=0.1, ph_freq_threshold=5.0,
                               save_per_pair=False, precision='s'):
    """
    Diagnostic version of `exc_resonant_raman_tensor_oneph`. Computes the
    Raman tensor at a SINGLE laser energy in pure numpy and returns every
    intermediate quantity that goes into it.

    The returned dict is intended to be fed to `save_raman_components`,
    which writes it both as a single .npz archive and as a folder of
    column-format .dat files suitable for gnuplot.

    Memory notes
    ------------
    The dominant arrays are the ``(nmodes, nexc, nexc)`` ones. To keep the
    footprint small this routine (a) stores only the absorption-form
    exciton-phonon matrix and applies its conjugate on the fly instead of
    keeping a second full copy, (b) frees the pair-intensity scratch as soon
    as it is consumed, (c) builds the optional per-pair amplitudes one mode
    at a time, and (d) works in single precision by default (see
    ``precision``). For ``nexc=1000, nmodes=6`` the base footprint drops from
    ~0.4 GB (double, old) to ~0.1 GB (single); with ``save_per_pair`` the
    per-pair peak drops roughly 6x.

    Parameters
    ----------
    laser_energy : float
        Single laser energy omega_L in eV.
    ph_energies, exc_energies, exc_dipoles, exc_ph_mat_el, n_kpts, cell_vol :
        Same as `exc_resonant_raman_tensor_oneph`.
    broad, ph_freq_threshold :
        Same as `exc_resonant_raman_tensor_oneph`.
    save_per_pair : bool, optional
        If True, include the full complex
        (nmodes, 3, 3, nexc, nexc) per-pair amplitudes in the dict.
        Memory ~ nmodes * 9 * nexc^2 * itemsize. Off by default.
    precision : {'s', 'd'}, optional
        Floating-point precision of the returned arrays. 's' = single
        (float32 / complex64, the default) halves the memory of every large
        array; 'd' = double (float64 / complex128) reproduces the original
        bit-for-bit-comparable values. Energies/denominators are always
        evaluated in double internally and only the stored arrays are cast.

    Returns
    -------
    dict
        See README produced by `save_raman_components` for the
        complete list of keys. (The redundant ``exc_ph_emi`` key is no longer
        included; it is simply ``conj(exc_ph_abs)``.)
    """
    cm1_to_Ha = 0.12398e-3 / ha2ev

    prec = str(precision).strip().lower()
    if prec in ('s', 'single'):
        f_type, c_type = np.float32, np.complex64
    elif prec in ('d', 'double'):
        f_type, c_type = np.float64, np.complex128
    else:
        raise ValueError("precision must be 's' (single) or 'd' (double)")

    # ----- Unit conversion (eV -> Ha) ------------------------------------
    wL           = float(laser_energy) / ha2ev
    ph_Ha        = np.asarray(ph_energies,  dtype=float) / ha2ev
    exc_Ha       = np.asarray(exc_energies, dtype=float) / ha2ev
    broad_Ha     = (broad / ha2ev) / 2.0
    ph_thresh_Ha = ph_freq_threshold * cm1_to_Ha
    ram_fac      = 1.0 / float(n_kpts) / np.sqrt(cell_vol)

    nmodes = ph_Ha.shape[0]
    nexc   = exc_Ha.shape[0]

    BS_energies = exc_Ha - 1j * broad_Ha   # E_lam - i*Gamma/2

    D_emi = np.asarray(exc_dipoles, dtype=complex)
    D_abs = np.conj(D_emi)

    inv_denom_res_v1  = 1.0 / (wL - BS_energies)                              # (nexc,)
    inv_denom_ares_v1 = 1.0 / (wL + BS_energies)
    inv_denom_res_v2  = 1.0 / (wL - BS_energies[None, :] - ph_Ha[:, None])    # (nmodes, nexc)
    inv_denom_ares_v2 = 1.0 / (wL + BS_energies[None, :] - ph_Ha[:, None])

    dipS_res   = (D_abs * inv_denom_res_v1[None, :]).astype(c_type)           # (3, nexc)
    dipS_ares  = (D_emi * inv_denom_ares_v1[None, :]).astype(c_type)
    dipSp_res  = (D_emi[None, :, :] * inv_denom_res_v2[:, None, :]).astype(c_type)   # (nmodes, 3, nexc)
    dipSp_ares = (D_abs[None, :, :] * inv_denom_ares_v2[:, None, :]).astype(c_type)

    # Big array: keep only the absorption form. The emission (Stokes) form is
    # just its conjugate and is applied on the fly, saving a full second
    # (nmodes, nexc, nexc) copy.
    exc_ph_abs = np.asarray(exc_ph_mat_el, dtype=c_type)     # yambopy raw (mode, init, fin)

    radiation_factor = np.sqrt(np.abs(wL - ph_Ha) / wL).astype(f_type)
    active_modes     = np.abs(ph_Ha) > ph_thresh_Ha
    prefactor        = (radiation_factor * ram_fac).astype(f_type)           # (nmodes,) real

    # Resonant term needs conj(exc_ph_abs). Using the identity
    #   conj( sum a*conj(g)*b ) = sum conj(a)*g*conj(b)
    # keeps the conjugations on the small dressed-dipole arrays and lets the
    # einsum read exc_ph_abs directly (no big (nmodes,nexc,nexc) conj temp).
    term_res  = np.conj(np.einsum('al, mLl, mbL -> mab',
                                  np.conj(dipS_res), exc_ph_abs, np.conj(dipSp_res),
                                  optimize=True))
    term_ares = np.einsum('al, mLl, mbL -> mab',
                          dipS_ares, exc_ph_abs, dipSp_ares, optimize=True)

    term_res  = (term_res  * prefactor[:, None, None]).astype(c_type)
    term_ares = (term_ares * prefactor[:, None, None]).astype(c_type)
    term_res[~active_modes]  = 0.0
    term_ares[~active_modes] = 0.0
    raman_tensor = term_res + term_ares

    # Per-pair factorized intensity (sum over polarizations); cheap & always returned.
    A_res     = np.sum(np.abs(dipS_res)**2,  axis=0)
    A_ares    = np.sum(np.abs(dipS_ares)**2, axis=0)
    B_res     = np.sum(np.abs(dipSp_res)**2,  axis=1)
    B_ares    = np.sum(np.abs(dipSp_ares)**2, axis=1)
    g_sq_pair = (np.abs(exc_ph_abs)**2).transpose(0, 2, 1).astype(f_type)
    pref_sq   = (prefactor**2)[:, None, None]

    pair_intensity_res  = (A_res[None, :, None]  * g_sq_pair * B_res[:, None, :] * pref_sq).astype(f_type)
    pair_intensity_ares = (A_ares[None, :, None] * g_sq_pair * B_ares[:, None, :] * pref_sq).astype(f_type)
    del g_sq_pair                                            # free the (nmodes, nexc, nexc) scratch
    pair_intensity_res[~active_modes]  = 0.0
    pair_intensity_ares[~active_modes] = 0.0

    components = {
        'laser_energy_eV'         : float(laser_energy),
        'broad_eV'                : float(broad),
        'ph_freq_threshold_cm-1'  : float(ph_freq_threshold),
        'n_kpts'                  : int(n_kpts),
        'cell_vol_bohr3'          : float(cell_vol),
        'ram_fac'                 : float(ram_fac),
        'exc_energies_eV'         : exc_Ha * ha2ev,
        'exc_energies_Ha'         : exc_Ha,
        'ph_energies_eV'          : ph_Ha * ha2ev,
        'ph_energies_Ha'          : ph_Ha,
        'active_modes'            : active_modes.astype(np.int64),
        'exc_dip_emi'             : D_emi,
        'exc_dip_absorp'          : D_abs,
        'exc_ph_abs'              : exc_ph_abs,
        'inv_denom_res_v1'        : inv_denom_res_v1,
        'inv_denom_ares_v1'       : inv_denom_ares_v1,
        'inv_denom_res_v2'        : inv_denom_res_v2,
        'inv_denom_ares_v2'       : inv_denom_ares_v2,
        'dipS_res'                : dipS_res,
        'dipS_ares'               : dipS_ares,
        'dipSp_res'               : dipSp_res,
        'dipSp_ares'              : dipSp_ares,
        'radiation_factor'        : radiation_factor,
        'pair_intensity_res'      : pair_intensity_res,
        'pair_intensity_ares'     : pair_intensity_ares,
        'term_res'                : term_res,
        'term_ares'               : term_ares,
        'raman_tensor'            : raman_tensor,
    }

    if save_per_pair:
        # Build one mode at a time into preallocated arrays: the einsum never
        # forms an all-modes (nmodes, 3, 3, nexc, nexc) temporary, so the peak
        # working set is a single mode (3, 3, nexc, nexc) plus one (nexc, nexc)
        # conjugate slice, rather than several full-size temporaries.
        per_pair_res  = np.zeros((nmodes, 3, 3, nexc, nexc), dtype=c_type)
        per_pair_ares = np.zeros((nmodes, 3, 3, nexc, nexc), dtype=c_type)
        for m in range(nmodes):
            if not active_modes[m]:
                continue
            g_m = exc_ph_abs[m]                              # (nexc, nexc) view
            per_pair_res[m]  = np.einsum('al, Ll, bL -> ablL',
                                         dipS_res, np.conj(g_m), dipSp_res[m],
                                         optimize=True) * prefactor[m]
            per_pair_ares[m] = np.einsum('al, Ll, bL -> ablL',
                                         dipS_ares, g_m, dipSp_ares[m],
                                         optimize=True) * prefactor[m]
        components['per_pair_res']  = per_pair_res
        components['per_pair_ares'] = per_pair_ares

    return components


def save_raman_components(components, out_dir='raman_components',
                          top_n_pairs=20, heatmap_modes='auto'):
    """
    Persist a diagnostic dict from `exc_raman_components_oneph` to disk.

    Writes one compressed .npz archive containing every array, plus a
    set of column-format .dat files for gnuplot, plus a README and a
    sample plot.gp.

    Parameters
    ----------
    components : dict
        As returned by `exc_raman_components_oneph`.
    out_dir : str
        Folder to create / overwrite.
    top_n_pairs : int
        Number of dominant (lambda1, lambda2) pairs listed per active mode.
    heatmap_modes : {'auto', 'all', None} or list of int
        Which modes get a full nexc x nexc pair-intensity heatmap file.
        'auto' = top 3 by |R|^2 (default). 'all' = every active mode.
        None = skip. Explicit list = those mode indices.
    """
    import os

    os.makedirs(out_dir, exist_ok=True)

    # 1. npz with everything
    array_components = {}
    for k, v in components.items():
        try:
            array_components[k] = np.asarray(v)
        except (TypeError, ValueError):
            pass
    np.savez_compressed(os.path.join(out_dir, 'components.npz'),
                        **array_components)

    wL_eV   = float(components['laser_energy_eV'])
    ph_eV   = np.asarray(components['ph_energies_eV'])
    exc_eV  = np.asarray(components['exc_energies_eV'])
    active  = np.asarray(components['active_modes']).astype(bool)
    rad_fac = np.asarray(components['radiation_factor'])
    nmodes  = ph_eV.shape[0]
    nexc    = exc_eV.shape[0]
    ph_cm   = ph_eV * _EV_TO_CM1

    term_res  = np.asarray(components['term_res'])
    term_ares = np.asarray(components['term_ares'])
    raman_t   = np.asarray(components['raman_tensor'])

    mod_res_sum  = np.sum(np.abs(term_res)**2,  axis=(1, 2))
    mod_ares_sum = np.sum(np.abs(term_ares)**2, axis=(1, 2))
    mod_total    = np.sum(np.abs(raman_t)**2,   axis=(1, 2))

    # 2. summary_modes.dat
    with open(os.path.join(out_dir, 'summary_modes.dat'), 'w') as f:
        f.write('# Excitonic Raman diagnostics at laser_energy = %.6f eV\n' % wL_eV)
        f.write('# 1:mode_idx  2:ph_freq_cm-1  3:ph_freq_eV  '
                '4:|R_res|^2_sum  5:|R_ares|^2_sum  6:|R_total|^2_sum  '
                '7:radiation_factor  8:active(1)/skipped(0)\n')
        for m in range(nmodes):
            f.write('%6d  %14.4f  %14.6e  %16.6e  %16.6e  %16.6e  %12.4e  %d\n'
                    % (m, ph_cm[m], ph_eV[m],
                       mod_res_sum[m], mod_ares_sum[m], mod_total[m],
                       rad_fac[m], int(active[m])))

    # 3. summary_excitons.dat
    D_emi         = np.asarray(components['exc_dip_emi'])
    inv_d_res_v1  = np.asarray(components['inv_denom_res_v1'])
    inv_d_ares_v1 = np.asarray(components['inv_denom_ares_v1'])
    Dx, Dy, Dz = np.abs(D_emi[0]), np.abs(D_emi[1]), np.abs(D_emi[2])
    D_sq = Dx**2 + Dy**2 + Dz**2

    with open(os.path.join(out_dir, 'summary_excitons.dat'), 'w') as f:
        f.write('# Excitonic Raman diagnostics at laser_energy = %.6f eV\n' % wL_eV)
        f.write('# 1:exc_idx  2:E_eV  3:|D_x|  4:|D_y|  5:|D_z|  6:|D|^2  '
                '7:|1/denom_res_v1|  8:|1/denom_ares_v1|\n')
        for l in range(nexc):
            f.write('%6d  %14.6f  %14.6e  %14.6e  %14.6e  %14.6e  %14.6e  %14.6e\n'
                    % (l, exc_eV[l], Dx[l], Dy[l], Dz[l], D_sq[l],
                       np.abs(inv_d_res_v1[l]), np.abs(inv_d_ares_v1[l])))

    # 4. dominant pairs per mode
    pair_res   = np.asarray(components['pair_intensity_res'])
    pair_ares  = np.asarray(components['pair_intensity_ares'])
    pair_total = pair_res + pair_ares
    exc_ph_abs = np.asarray(components['exc_ph_abs'])
    g_abs_sq   = np.abs(exc_ph_abs)**2

    pairs_dir = os.path.join(out_dir, 'pairs_per_mode')
    os.makedirs(pairs_dir, exist_ok=True)
    for m in range(nmodes):
        if not active[m]:
            continue
        flat = pair_total[m].ravel()
        if not np.any(flat > 0):
            continue
        n_top = min(top_n_pairs, flat.size)
        top_idx = np.argsort(flat)[::-1][:n_top]
        L1, L2 = np.unravel_index(top_idx, pair_total[m].shape)
        filename = os.path.join(pairs_dir, 'pairs_mode_%03d.dat' % m)
        with open(filename, 'w') as f:
            f.write('# Top %d pair contributions, mode %d '
                    '(omega = %.2f cm-1, %.4f eV)\n'
                    % (n_top, m, ph_cm[m], ph_eV[m]))
            f.write('# laser_energy = %.6f eV\n' % wL_eV)
            f.write('# lambda1 = first-vertex exciton; lambda2 = second-vertex exciton\n')
            f.write('# 1:rank  2:lambda1  3:lambda2  4:E_lambda1_eV  5:E_lambda2_eV  '
                    '6:|g(lambda2,lambda1)|^2  7:|c_res|^2  8:|c_ares|^2  9:|c_total|^2\n')
            for rank, (l1, l2) in enumerate(zip(L1, L2), start=1):
                f.write('%6d  %6d  %6d  %14.6f  %14.6f  %14.6e  '
                        '%14.6e  %14.6e  %14.6e\n'
                        % (rank, l1, l2,
                           exc_eV[l1], exc_eV[l2],
                           g_abs_sq[m, l2, l1],
                           pair_res[m, l1, l2],
                           pair_ares[m, l1, l2],
                           pair_total[m, l1, l2]))

    # 5. heatmap files
    if heatmap_modes is None:
        chosen = []
    elif isinstance(heatmap_modes, str):
        opt = heatmap_modes.strip().lower()
        if opt == 'all':
            chosen = [m for m in range(nmodes) if active[m]]
        elif opt == 'auto':
            ranked = np.argsort(mod_total)[::-1]
            chosen = [int(m) for m in ranked if active[m]][:3]
        else:
            raise ValueError("heatmap_modes string must be 'auto' or 'all'")
    else:
        chosen = [int(m) for m in heatmap_modes]
    if chosen:
        heatmap_dir = os.path.join(out_dir, 'heatmaps')
        os.makedirs(heatmap_dir, exist_ok=True)
        for m in chosen:
            filename = os.path.join(heatmap_dir, 'pair_intensity_mode_%03d.dat' % m)
            with open(filename, 'w') as f:
                f.write('# Pair intensity heatmap, mode %d '
                        '(omega = %.2f cm-1, %.4f eV)\n'
                        % (m, ph_cm[m], ph_eV[m]))
                f.write('# laser_energy = %.6f eV\n' % wL_eV)
                f.write('# gnuplot:  splot "%s" u 1:2:5 with image\n'
                        % os.path.basename(filename))
                f.write('# 1:lambda1  2:lambda2  3:|c_res|^2  4:|c_ares|^2  5:|c_total|^2\n')
                for l1 in range(nexc):
                    for l2 in range(nexc):
                        f.write('%6d  %6d  %14.6e  %14.6e  %14.6e\n'
                                % (l1, l2,
                                   pair_res[m, l1, l2],
                                   pair_ares[m, l1, l2],
                                   pair_total[m, l1, l2]))
                    f.write('\n')

    # 6. README and plot.gp
    with open(os.path.join(out_dir, 'README.txt'), 'w') as f:
        f.write(
            'Excitonic Raman diagnostic dump\n'
            'Generated by save_raman_components(...).\n\n'
            'Laser energy : %.6f eV\n'
            'Broadening   : %.4f eV (electronic Gamma, full width)\n'
            'Nmodes       : %d   Nexc : %d\n\n'
            'Files\n-----\n'
            'components.npz\n'
            '    Compressed numpy archive with every array of the diagnostic\n'
            '    dict (use np.load(...) in Python).\n\n'
            'summary_modes.dat\n'
            '    One row per phonon mode. Columns:\n'
            '    1 mode_idx | 2 ph_freq (cm-1) | 3 ph_freq (eV)\n'
            '    4 |R_res|^2_sum_ab | 5 |R_ares|^2_sum_ab | 6 |R_total|^2_sum_ab\n'
            '    7 radiation_factor sqrt(|wL-wph|/wL) | 8 active(1)/skipped(0)\n\n'
            'summary_excitons.dat\n'
            '    One row per exciton state. Columns:\n'
            '    1 exc_idx | 2 E (eV) | 3 |Dx| | 4 |Dy| | 5 |Dz| | 6 |D|^2\n'
            '    7 |1/(wL - E + iG/2)| | 8 |1/(wL + E - iG/2)|\n\n'
            'pairs_per_mode/pairs_mode_NNN.dat\n'
            '    Top-N most contributing (lambda1, lambda2) pairs per active\n'
            '    mode. Ranked by |c_res|^2 + |c_ares|^2.\n'
            '    Columns: rank, lambda1, lambda2, E1, E2, |g|^2, |c_res|^2,\n'
            '             |c_ares|^2, |c_total|^2\n\n'
            'heatmaps/pair_intensity_mode_NNN.dat\n'
            '    Full (nexc x nexc) pair-intensity table for selected modes,\n'
            '    in gnuplot matrix-with-coordinates format (blank line between\n'
            '    rows). Columns: lambda1, lambda2, |c_res|^2, |c_ares|^2,\n'
            '    |c_total|^2\n\n'
            'Conventions\n-----------\n'
            '  lambda1 = first-vertex exciton (couples to incoming photon)\n'
            '  lambda2 = second-vertex exciton (couples to outgoing photon)\n'
            '  exc_ph[m, init, fin] = yambopy raw absorption-form matrix\n'
            '  element. The Raman formula uses exc_ph[m, lambda2, lambda1].\n\n'
            'See plot.gp in this folder for example gnuplot commands.\n'
            % (wL_eV, float(components['broad_eV']), nmodes, nexc))

    with open(os.path.join(out_dir, 'plot.gp'), 'w') as f:
        f.write(
            '# Quick gnuplot recipe for the Raman diagnostic dump.\n'
            '# Run inside this folder:  gnuplot plot.gp\n\n'
            'set terminal pdfcairo size 9in,7in enhanced font "Helvetica,11"\n'
            'set output "diagnostic.pdf"\n'
            'set grid\n\n'
            'set multiplot layout 2,2 title '
            '"Excitonic Raman diagnostic at omega_L = %.3f eV"\n\n'
            '# (1) Mode-resolved spectrum\n'
            'set title "Per-mode |R|^2 (sum over polarisations)"\n'
            'set xlabel "Raman shift (cm^{-1})"\n'
            'set ylabel "|R|^2 (arb.u.)"\n'
            'plot "summary_modes.dat" u 2:6 w impulses lw 2 t "total", \\\n'
            '     ""                   u 2:4 w impulses lw 1 lt 3 t "resonant only"\n\n'
            '# (2) Exciton oscillator strengths\n'
            'set title "Exciton |D|^2 vs E"\n'
            'set xlabel "Exciton energy (eV)"\n'
            'set ylabel "|D|^2"\n'
            'plot "summary_excitons.dat" u 2:6 w impulses lw 1.5 notitle\n\n'
            '# (3) Resonance pattern\n'
            'set title "Resonance: |1/(omega_L - E +/- iGamma/2)|"\n'
            'set xlabel "Exciton energy (eV)"\n'
            'set ylabel "|1/denom|"\n'
            'plot "summary_excitons.dat" u 2:7 w l lw 1.5 t "resonant", \\\n'
            '     ""                     u 2:8 w l lw 1   t "anti-resonant"\n\n'
            '# (4) Top-pair contributions for the brightest mode\n'
            'set title "Top 20 (lambda1, lambda2) contributions"\n'
            'set xlabel "Rank"\n'
            'set ylabel "|c_total|^2"\n'
            'set logscale y\n'
            'plot for [f in system("ls pairs_per_mode/pairs_mode_*.dat | head -1")] \\\n'
            '     f u 1:9 w impulses lw 2 t f\n\n'
            'unset multiplot\n'
            'unset output\n\n'
            '# Pair-intensity heatmap (if heatmaps/ exists), uncomment:\n'
            '#   set terminal pdfcairo size 6in,5in\n'
            '#   set output "heatmap.pdf"\n'
            '#   set title "|c_total|^2(lambda1, lambda2)"\n'
            '#   set xlabel "lambda1"; set ylabel "lambda2"\n'
            '#   set palette defined (0 "white", 0.5 "orange", 1 "red")\n'
            '#   set logscale cb\n'
            '#   splot "heatmaps/pair_intensity_mode_012.dat" u 1:2:5 w image\n'
            % wL_eV)


def save_resonance_heatmap(components, mode_idx, out_path=None):
    """
    Write a (nexc x nexc) heatmap of the bare resonance denominators for a
    single phonon mode — no dipoles, no phonon matrix elements.

    Each cell (l1, l2) contains:

        res_map[l1, l2]  = |1/(wL - E_l1 + iG/2)| * |1/(wL - E_l2 - wph + iG/2)|
        ares_map[l1, l2] = |1/(wL + E_l1 - iG/2)| * |1/(wL + E_l2 - wph - iG/2)|

    This shows the pure double-resonance landscape: which (l1, l2) pairs are
    simultaneously resonant at both photon vertices, irrespective of whether
    they are optically bright or phonon-connected.

    Parameters
    ----------
    components : dict
        As returned by `exc_raman_components_oneph`.
    mode_idx : int
        Index of the phonon mode to plot.
    out_path : str or None
        File path for the output .dat file.  If None, defaults to
        'resonance_heatmap_mode_<mode_idx>.dat' in the current directory.

    Returns
    -------
    res_map  : (nexc, nexc) float ndarray
    ares_map : (nexc, nexc) float ndarray
    """
    import os

    inv_v1_res  = np.abs(np.asarray(components['inv_denom_res_v1']))   # (nexc,)
    inv_v1_ares = np.abs(np.asarray(components['inv_denom_ares_v1']))
    inv_v2_res  = np.abs(np.asarray(components['inv_denom_res_v2']))   # (nmodes, nexc)
    inv_v2_ares = np.abs(np.asarray(components['inv_denom_ares_v2']))

    m = int(mode_idx)
    # outer product: axis-0 = l1 (first vertex), axis-1 = l2 (second vertex)
    res_map  = np.outer(inv_v1_res,  inv_v2_res[m])   # (nexc, nexc)
    ares_map = np.outer(inv_v1_ares, inv_v2_ares[m])

    ph_eV  = np.asarray(components['ph_energies_eV'])
    wL_eV  = float(components['laser_energy_eV'])
    exc_eV = np.asarray(components['exc_energies_eV'])
    nexc   = exc_eV.shape[0]

    if out_path is None:
        out_path = 'resonance_heatmap_mode_%03d.dat' % m

    with open(out_path, 'w') as f:
        f.write('# Bare resonance heatmap — mode %d '
                '(omega = %.2f cm-1, %.4f eV)\n'
                % (m, ph_eV[m] * _EV_TO_CM1, ph_eV[m]))
        f.write('# laser_energy = %.6f eV\n' % wL_eV)
        f.write('# No dipoles, no phonon matrix elements.\n')
        f.write('# res_map[l1,l2]  = |1/(wL-E_l1+iG/2)| * |1/(wL-E_l2-wph+iG/2)|\n')
        f.write('# ares_map[l1,l2] = |1/(wL+E_l1-iG/2)| * |1/(wL+E_l2-wph-iG/2)|\n')
        f.write('# gnuplot:  splot "%s" u 1:2:3 with image   (res)\n'
                % os.path.basename(out_path))
        f.write('# gnuplot:  splot "%s" u 1:2:4 with image   (ares)\n'
                % os.path.basename(out_path))
        f.write('# 1:l1  2:l2  3:res_map  4:ares_map  5:E_l1_eV  6:E_l2_eV\n')
        for l1 in range(nexc):
            for l2 in range(nexc):
                f.write('%6d  %6d  %14.6e  %14.6e  %14.6f  %14.6f\n'
                        % (l1, l2,
                           res_map[l1, l2], ares_map[l1, l2],
                           exc_eV[l1], exc_eV[l2]))
            f.write('\n')

    return res_map, ares_map


def save_dipole_g_dipole_heatmap(components, mode_idx, out_path=None):
    """
    Write a (nexc x nexc) heatmap of the bare dipole-phonon-dipole product for
    a single phonon mode — no energy denominators.

    Each cell (l1, l2) contains:

        dip_g_dip[l1, l2] = |D|²[l1]  *  |g[m, l2, l1]|²  *  |D|²[l2]

    where |D|²[l] = Σ_a |D_emi[a, l]|² is the total oscillator strength of
    exciton l (summed over Cartesian components), and |g[m, l2, l1]|² is the
    phonon matrix element squared for mode m connecting l1 -> l2.

    This shows which pairs are simultaneously optically bright and
    phonon-connected, regardless of their resonance with the laser.
    Comparing this with the full pair-intensity heatmap reveals how much of
    the structure comes from the resonance denominators vs. the intrinsic
    optical / phonon coupling.

    Parameters
    ----------
    components : dict
        As returned by `exc_raman_components_oneph`.
    mode_idx : int
        Index of the phonon mode to plot.
    out_path : str or None
        File path for the output .dat file.  If None, defaults to
        'dip_g_dip_heatmap_mode_<mode_idx>.dat' in the current directory.

    Returns
    -------
    dip_g_dip : (nexc, nexc) float ndarray
    """
    import os

    D_emi    = np.asarray(components['exc_dip_emi'])    # (3, nexc)
    exc_ph   = np.asarray(components['exc_ph_abs'])     # (nmodes, nexc, nexc) [mode, init, fin]
    exc_eV   = np.asarray(components['exc_energies_eV'])
    ph_eV    = np.asarray(components['ph_energies_eV'])
    wL_eV    = float(components['laser_energy_eV'])

    m   = int(mode_idx)
    nexc = exc_eV.shape[0]

    D_sq      = np.sum(np.abs(D_emi)**2, axis=0)        # (nexc,)  |D|² per exciton
    g_sq      = np.abs(exc_ph[m])**2                    # (nexc, nexc) [init, fin]
    # g_sq[l1, l2] = |g[m, init=l1, fin=l2]|²
    # We want dip_g_dip[l1, l2] = D_sq[l1] * g_sq[l2, l1] * D_sq[l2]
    # (l1 is first-vertex, l2 is second-vertex, phonon connects l1->l2
    #  in the yambopy raw convention [mode, init, fin] where init=l2, fin=l1
    #  matches exc_ph_emi[m, l2, l1] used in term_res)
    dip_g_dip = D_sq[:, None] * g_sq.T * D_sq[None, :]  # (nexc, nexc)

    if out_path is None:
        out_path = 'dip_g_dip_heatmap_mode_%03d.dat' % m

    with open(out_path, 'w') as f:
        f.write('# Dipole-phonon-dipole heatmap — mode %d '
                '(omega = %.2f cm-1, %.4f eV)\n'
                % (m, ph_eV[m] * _EV_TO_CM1, ph_eV[m]))
        f.write('# laser_energy = %.6f eV  (no denominators)\n' % wL_eV)
        f.write('# dip_g_dip[l1,l2] = |D|^2[l1] * |g[m,l2,l1]|^2 * |D|^2[l2]\n')
        f.write('# gnuplot:  splot "%s" u 1:2:3 with image\n'
                % os.path.basename(out_path))
        f.write('# 1:l1  2:l2  3:dip_g_dip  4:|D|^2_l1  5:|g|^2  6:|D|^2_l2  '
                '7:E_l1_eV  8:E_l2_eV\n')
        for l1 in range(nexc):
            for l2 in range(nexc):
                f.write('%6d  %6d  %14.6e  %14.6e  %14.6e  %14.6e  %14.6f  %14.6f\n'
                        % (l1, l2,
                           dip_g_dip[l1, l2],
                           D_sq[l1],
                           g_sq[l2, l1],
                           D_sq[l2],
                           exc_eV[l1], exc_eV[l2]))
            f.write('\n')

    return dip_g_dip


# =============================================================================
# Interference between exciton pairs
#
# Every (lambda1, lambda2) pair adds a complex, polarisation-resolved
# amplitude c[a, b, l1, l2] to the Raman tensor R[a, b] = sum_pairs c.
# The measured intensity is |sum c|^2 (coherent), whereas the pair-intensity
# heatmaps show sum |c|^2 (incoherent). The difference is interference.
#
# c has rank-1 structure in (a, b) for each pair,
#     c = p * ( conj(g) * dipS_res[:, l1] (x) dipSp_res[:, l2]
#             +      g  * dipS_ares[:, l1] (x) dipSp_ares[:, l2] ),
# so every quantity below reduces to (nexc x nexc) matrix algebra and is
# evaluated in row blocks. The (nmodes, 3, 3, nexc, nexc) per-pair tensor is
# never built, which keeps this usable at thousands of exciton states.
# =============================================================================

def _select_interference_modes(modes, active, I_coh):
    """Resolve the `modes` argument of `exc_raman_interference_oneph`."""
    import warnings

    act = np.flatnonzero(active)
    if modes is None or (isinstance(modes, str) and modes.strip().lower() == 'auto'):
        ranked = act[np.argsort(I_coh[act])[::-1]]
        return [int(m) for m in ranked[:3]]
    if isinstance(modes, str):
        if modes.strip().lower() == 'all':
            return [int(m) for m in act]
        raise ValueError("modes must be 'auto', 'all', None or a list of "
                         "mode indices (got %r)" % modes)

    sel = [int(m) for m in np.atleast_1d(modes)]
    bad = [m for m in sel if m < 0 or m >= len(active)]
    if bad:
        raise ValueError("mode indices %s out of range 0..%d"
                         % (bad, len(active) - 1))
    skipped = [m for m in sel if not active[m]]
    if skipped:
        warnings.warn("modes %s are below the acoustic threshold and carry "
                      "no Raman tensor; skipped." % skipped, stacklevel=3)
    return [m for m in sel if active[m]]


def _select_mode_groups(mode_groups, active, ph_eV, tol_cm1):
    """Resolve the `mode_groups` argument of `exc_raman_interference_oneph`."""
    import warnings

    if mode_groups is None or isinstance(mode_groups, str):
        if isinstance(mode_groups, str):
            raise ValueError("mode_groups must be a list of lists of mode "
                             "indices, e.g. [[7, 8]] (got %r)" % mode_groups)
        return []
    groups = list(mode_groups)
    if not groups:
        return []
    if all(np.ndim(g) == 0 for g in groups):          # flat list = one group
        groups = [groups]

    out = []
    for g in groups:
        g = [int(m) for m in np.atleast_1d(g)]
        if not g:
            raise ValueError("mode_groups contains an empty group")
        if len(set(g)) != len(g):
            raise ValueError("mode group %s contains a mode twice" % g)
        bad = [m for m in g if m < 0 or m >= len(active)]
        if bad:
            raise ValueError("mode group %s: indices %s out of range 0..%d"
                             % (g, bad, len(active) - 1))
        dead = [m for m in g if not active[m]]
        if dead:
            raise ValueError("mode group %s: modes %s are below the acoustic "
                             "threshold and carry no Raman tensor" % (g, dead))
        spread = (ph_eV[g].max() - ph_eV[g].min()) * _EV_TO_CM1
        if spread > tol_cm1:
            warnings.warn("mode group %s spans %.2f cm-1 (> %.2f cm-1); are "
                          "these really degenerate partners?"
                          % (g, spread, tol_cm1), stacklevel=3)
        out.append(tuple(sorted(g)))
    return out


def exc_raman_interference_oneph(components, modes='auto', group_edges=None,
                                 n_groups=8, curve_points=500,
                                 curve_max_pairs=None, block_elems=2_000_000,
                                 mode_groups=None, group_freq_tol_cm1=1.0):
    """
    Quantify interference between exciton pairs in the excitonic Raman tensor.

    Works on the dict returned by `exc_raman_components_oneph` (it does NOT
    need ``save_per_pair=True``) and reports, per phonon mode, how the
    complex amplitudes of the (lambda1, lambda2) pairs combine.

    Definitions (mode m, polarisations a, b)
    -----------------------------------------
    c[a,b,l1,l2]      amplitude of one pair;  R[a,b] = sum_pairs c  (the tensor)
    I_coh             sum_ab |R|^2                 measured intensity
    I_incoh           sum_pairs sum_ab |c|^2       intensity without interference
    coherence_ratio   I_coh / I_incoh              1: independent pairs,
                                                   <1: destructive, >1: constructive

    Resonant and anti-resonant parts of the SAME pair are added coherently
    inside c; "incoherent" refers to adding different pairs.

    Parameters
    ----------
    components : dict
        As returned by `exc_raman_components_oneph` (or re-loaded from its
        components.npz). Needs the dressed dipoles, `exc_ph_abs`,
        `radiation_factor`, `ram_fac` and `active_modes`.
    modes : 'auto', 'all', None or list of int, optional
        Modes for which the maps, curve and blocks are built. 'auto'/None =
        the 3 brightest active modes. The scalar ratios are always computed
        for every active mode. Indices are 0-based.
    group_edges : 1D array of float, optional
        Exciton-energy edges (eV) of the windows used for the block
        decomposition; windows are half-open [e_i, e_i+1). States outside all
        windows are reported through ``residual_share``.
    n_groups : int or None, optional
        If `group_edges` is not given, use this many uniform windows spanning
        all exciton energies. None/0 skips the block decomposition. Default 8.
    curve_points : int, optional
        Number of (log-spaced) samples stored for the cumulative curve.
    curve_max_pairs : int or None, optional
        Use only the strongest this-many pairs for the cumulative curve. The
        full ordering needs an argsort over nexc^2 pairs (~16 nexc^2 bytes);
        capping it saves memory. The curve then does not reach I_coh
        (``curve_complete`` is False).
    block_elems : int, optional
        Approximate number of matrix elements processed per row block. Lower
        it to reduce peak memory; results do not depend on it.
    mode_groups : list of lists of int, optional
        Degenerate phonon partners to analyse together, e.g. [[7, 8]] for an
        E-type mode (a flat list such as [7, 8] is one group). Different
        phonon modes are different final states, so a group adds the modes
        INCOHERENTLY: I = sum_m I_m, and every map uses
        sum_m Re(sum_ab conj(R_m) c_m) and sum_m sum_ab |c_m|^2. The result is
        independent of the eigenvector basis chosen inside the degenerate
        subspace, unlike the single-mode maps. Group members need not be
        listed in `modes`.
    group_freq_tol_cm1 : float, optional
        Warn if the phonon frequencies inside a group differ by more than
        this (cm-1). Default 1.0.

    Returns
    -------
    dict with keys
        I_coh, I_incoh, coherence_ratio : (nmodes,) float  (NaN if inactive)
        analysed_modes                  : (n,) int
        group_edges_eV                  : edges used, or None
        per_mode : {mode: dict} with, for each analysed mode,
            phase_weight      (nexc, nexc)  Re(sum_ab conj(R) c) / I_coh.
                              Signed share of the intensity carried by each
                              pair; sums to exactly 1. Negative = cancels.
            incoherent_share  (nexc, nexc)  sum_ab|c|^2 / I_incoh; sums to 1.
                              Share the pair would have without interference.
            interference      (nexc, nexc)  phase_weight - sum_ab|c|^2 / I_coh.
                              Half of the pair's cross terms with all other
                              pairs, over I_coh: >0 constructive, <0
                              destructive. Sums to 1 - I_incoh/I_coh. Same
                              convention as block_interference.
            curve_n_pairs, curve_I_coh, curve_I_incoh : running intensities
                              when adding pairs from strongest to weakest.
            curve_complete    whether the curve covers every pair.
            block_intensity   (ng, ng) sum_ab |A_ij|^2, block (i, j) alone.
            block_share       (ng, ng) Re(sum_ab conj(R) A_ij) / I_coh.
            block_interference (ng, ng) block_share - block_intensity/I_coh:
                              interference of block (i, j) with all others.
            block_coherence   I_coh / sum_ij block_intensity.
            residual_share    share carried by states outside the windows.
        mode_groups                     : list of tuples, groups analysed
        per_group : {label: dict}, label like '007_008', with the same keys
            as per_mode (sums over the group's modes as described under
            `mode_groups`) plus modes, I_coh, I_incoh, coherence_ratio.
        plus laser_energy_eV, broad_eV, exc_energies_eV, ph_energies_eV,
        active_modes copied from `components`.

    Notes
    -----
    Index convention matches the heatmaps: l1 = first-vertex exciton (rows),
    l2 = second-vertex exciton (columns), g = exc_ph_abs[m, l2, l1].

    Memory per analysed mode is ~16 nexc^2 bytes for the two maps plus ~16
    nexc^2 transiently for the curve ordering (0.4 + 0.8 GB at 7200 states).

    Interference means subtracting nearly equal amplitudes; run
    `exc_raman_components_oneph` with ``precision='d'`` for quantitative
    results. A warning is issued for single-precision input.
    """
    import warnings

    dSr = np.asarray(components['dipS_res'])
    if dSr.dtype == np.complex64:
        warnings.warn(
            "components were produced with precision='s'. Interference means "
            "subtracting nearly equal amplitudes, so single precision can "
            "exaggerate weak cancellation; rerun exc_raman_components_oneph "
            "with precision='d' for quantitative results.", stacklevel=2)

    dSr    = dSr.astype(np.complex128)                              # (3, nexc)
    dSa    = np.asarray(components['dipS_ares'], dtype=np.complex128)
    dSpr   = np.asarray(components['dipSp_res'])                    # (nmodes, 3, nexc)
    dSpa   = np.asarray(components['dipSp_ares'])
    G      = np.asarray(components['exc_ph_abs'])                   # (nmodes, nexc, nexc) [m, l2, l1]
    active = np.asarray(components['active_modes']).astype(bool)
    pref   = (np.asarray(components['radiation_factor'], dtype=float)
              * float(components['ram_fac']))
    exc_eV = np.asarray(components['exc_energies_eV'], dtype=float)
    ph_eV  = np.asarray(components['ph_energies_eV'], dtype=float)
    nmodes, nexc = G.shape[0], G.shape[1]

    rb = max(1, min(nexc, int(block_elems) // max(nexc, 1)))        # rows per block

    # Mode-independent first-vertex factors.
    Ar = np.sum(np.abs(dSr)**2, axis=0)
    Aa = np.sum(np.abs(dSa)**2, axis=0)
    Sx = np.sum(np.conj(dSr) * dSa, axis=0)                          # res-ares overlap

    def _gblock(m, s, e):
        """Rows s:e of gT[l1, l2] = G[m, l2, l1], as a double-precision copy."""
        return np.ascontiguousarray(G[m][:, s:e].T, dtype=np.complex128)

    def _second_vertex(m):
        return (dSpr[m].astype(np.complex128), dSpa[m].astype(np.complex128))

    def _tensor(m):
        """R[a,b] = sum_pairs c, accumulated block-wise in double precision."""
        dpr, dpa = _second_vertex(m)
        R = np.zeros((3, 3), dtype=np.complex128)
        for s in range(0, nexc, rb):
            e  = min(s + rb, nexc)
            gb = _gblock(m, s, e)
            R += dSr[:, s:e] @ (np.conj(gb) @ dpr.T) + dSa[:, s:e] @ (gb @ dpa.T)
        return pref[m] * R

    def _pair_pass(m, R, want_maps, P):
        """Incoherent total; optionally the pair maps and block accumulators."""
        p = pref[m]
        dpr, dpa = _second_vertex(m)
        Br = np.sum(np.abs(dpr)**2, axis=0)
        Ba = np.sum(np.abs(dpa)**2, axis=0)
        Sy = np.sum(np.conj(dpr) * dpa, axis=0)
        Rc = np.conj(R)

        inc = w = Zr = Za = None
        if want_maps:
            inc = np.empty((nexc, nexc))
            w   = np.empty((nexc, nexc))
        if P is not None:
            Zr = np.zeros((P.shape[0], 3, nexc), dtype=np.complex128)
            Za = np.zeros_like(Zr)

        I_inc = 0.0
        for s in range(0, nexc, rb):
            e  = min(s + rb, nexc)
            gb = _gblock(m, s, e)
            # sum_ab |c|^2 = p^2 [ |g|^2 (Ar Br + Aa Ba) + 2 Re(g^2 Sx Sy) ]
            incb  = (gb.real**2 + gb.imag**2) * (np.outer(Ar[s:e], Br)
                                                 + np.outer(Aa[s:e], Ba))
            incb += 2.0 * np.real(gb * gb * np.outer(Sx[s:e], Sy))
            incb *= p * p
            I_inc += incb.sum()
            if want_maps:
                inc[s:e] = incb
                # sum_ab conj(R) c = p [ conj(g) u_res + g u_ares ]
                ur = dSr[:, s:e].T @ Rc @ dpr
                ua = dSa[:, s:e].T @ Rc @ dpa
                w[s:e] = p * np.real(np.conj(gb) * ur + gb * ua)
            if P is not None:
                Zr += (P[:, None, s:e] * dSr[None, :, s:e]) @ np.conj(gb)
                Za += (P[:, None, s:e] * dSa[None, :, s:e]) @ gb
        return I_inc, inc, w, Zr, Za

    def _curve(ms, inc):
        """
        Running coherent / incoherent intensity, strongest pairs first.
        `ms` are modes added incoherently: pairs are ordered by `inc` (their
        own intensity summed over ms) and the running coherent intensity is
        sum_m sum_ab |sum_first-n c_m|^2.
        """
        vert   = [(pref[m],) + _second_vertex(m) + (G[m],) for m in ms]
        flat   = inc.ravel()
        npairs = flat.size
        if curve_max_pairs is not None and int(curve_max_pairs) < npairs:
            k     = max(1, int(curve_max_pairs))
            part  = np.argpartition(flat, npairs - k)[npairs - k:]
            order = part[np.argsort(flat[part])[::-1]]
        else:
            order = np.argsort(flat)[::-1]
        nuse = order.size

        npts = max(2, min(int(curve_points), nuse))
        ks   = np.unique(np.append(
            np.geomspace(1, nuse, num=npts).astype(np.int64), nuse))
        I_coh_k = np.empty(ks.size)
        I_inc_k = np.empty(ks.size)

        chunk = max(1, int(block_elems) // (9 * len(ms)))
        S = [np.zeros((3, 3), dtype=np.complex128) for _ in ms]
        run_inc, ptr = 0.0, 0
        for s in range(0, nuse, chunk):
            idx = order[s:s + chunk]
            l1, l2 = idx // nexc, idx % nexc
            cs_all = []
            for k, (p, dpr, dpa, Gm) in enumerate(vert):
                g = Gm[l2, l1].astype(np.complex128)
                c = (np.conj(g)[:, None, None] * dSr[:, l1].T[:, :, None] * dpr[:, l2].T[:, None, :]
                     + g[:, None, None]        * dSa[:, l1].T[:, :, None] * dpa[:, l2].T[:, None, :])
                c *= p
                cs = np.cumsum(c, axis=0)
                cs += S[k]
                cs_all.append(cs)
            run = run_inc + np.cumsum(flat[idx])
            j0 = ptr
            while ptr < ks.size and ks[ptr] <= s + idx.size:
                ptr += 1
            if ptr > j0:
                pos = ks[j0:ptr] - (s + 1)
                I_coh_k[j0:ptr] = sum(np.sum(np.abs(cs[pos])**2, axis=(1, 2))
                                      for cs in cs_all)
                I_inc_k[j0:ptr] = run[pos]
            S, run_inc = [cs[-1] for cs in cs_all], run[-1]
        return ks, I_coh_k, I_inc_k, nuse == npairs

    # ---- energy windows for the block decomposition ----
    P = edges = None
    if group_edges is not None:
        edges = np.asarray(group_edges, dtype=float)
        if edges.ndim != 1 or edges.size < 2 or np.any(np.diff(edges) <= 0):
            raise ValueError("group_edges must be a strictly increasing 1D "
                             "array with at least 2 values (eV)")
    elif n_groups:
        ng = int(n_groups)
        if ng < 1:
            raise ValueError("n_groups must be a positive integer or None")
        lo, hi = float(exc_eV.min()), float(exc_eV.max())
        if hi <= lo:
            hi = lo + 1e-6
        edges = np.linspace(lo, hi, ng + 1)
        edges[-1] = np.nextafter(hi, np.inf)        # keep the top state inside
    if edges is not None:
        ng  = edges.size - 1
        gid = np.digitize(exc_eV, edges) - 1
        ok  = (gid >= 0) & (gid < ng)
        P = np.zeros((ng, nexc))
        P[gid[ok], np.flatnonzero(ok)] = 1.0

    groups = _select_mode_groups(mode_groups, active, ph_eV, group_freq_tol_cm1)

    # ---- tensor and coherent intensity for every active mode ----
    I_coh = np.zeros(nmodes)
    tensors = {}
    for m in np.flatnonzero(active):
        tensors[m] = _tensor(m)
        I_coh[m] = np.sum(np.abs(tensors[m])**2)

    sel = _select_interference_modes(modes, active, I_coh)
    I_inc = np.zeros(nmodes)

    def _analyse(ms):
        """
        Maps, cumulative curve and blocks for the modes `ms` added
        incoherently (a single mode is the one-element case). Also fills
        I_inc for these modes. Returns (entry, I_coh, I_incoh) of the set.
        """
        Ic = float(sum(I_coh[m] for m in ms))
        w = inc = b_int = b_raw = None
        for m in ms:
            R = tensors[m]
            I_inc[m], inc_m, w_m, Zr, Za = _pair_pass(m, R, True, P)
            if w is None:
                w, inc = w_m, inc_m
            else:
                w += w_m
                inc += inc_m
            del w_m, inc_m
            if P is not None:
                dpr, dpa = _second_vertex(m)
                A = pref[m] * (np.einsum('iaL,bL,jL->ijab', Zr, dpr, P)
                               + np.einsum('iaL,bL,jL->ijab', Za, dpa, P))
                a_int = np.sum(np.abs(A)**2, axis=(2, 3))
                a_raw = np.real(np.einsum('ab,ijab->ij', np.conj(R), A))
                b_int = a_int if b_int is None else b_int + a_int
                b_raw = a_raw if b_raw is None else b_raw + a_raw
        Ii = float(sum(I_inc[m] for m in ms))

        entry = {}
        ks, curve_coh, curve_inc, complete = _curve(ms, inc)
        entry.update(curve_n_pairs=ks, curve_I_coh=curve_coh,
                     curve_I_incoh=curve_inc, curve_complete=complete)

        if Ic > 0:
            w /= Ic
            intf = w - inc / Ic                     # half the cross terms / I_coh
        else:
            w[:] = np.nan
            intf = np.full_like(w, np.nan)
        if Ii > 0:
            inc /= Ii
        else:
            inc[:] = np.nan
        entry['phase_weight']     = w
        entry['incoherent_share'] = inc
        entry['interference']     = intf

        if P is not None:
            if Ic > 0:
                b_share = b_raw / Ic
                b_intf  = b_share - b_int / Ic
                resid   = 1.0 - b_share.sum()
            else:
                b_share = b_intf = np.full_like(b_int, np.nan)
                resid   = np.nan
            tot = b_int.sum()
            entry.update(block_intensity=b_int, block_share=b_share,
                         block_interference=b_intf,
                         block_coherence=(Ic / tot) if tot > 0 else np.nan,
                         residual_share=resid)
        return entry, Ic, Ii

    in_groups = {m for g in groups for m in g}
    per_mode = {}
    for m in np.flatnonzero(active):
        if m in sel:
            per_mode[int(m)] = _analyse([m])[0]
        elif m not in in_groups:
            I_inc[m] = _pair_pass(m, tensors[m], False, None)[0]

    per_group = {}
    for g in groups:
        entry, Ic, Ii = _analyse(list(g))
        entry.update(modes=np.array(g, dtype=np.int64), I_coh=Ic, I_incoh=Ii,
                     coherence_ratio=(Ic / Ii) if Ii > 0 else np.nan)
        per_group['_'.join('%03d' % m for m in g)] = entry

    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np.where(I_inc > 0, I_coh / I_inc, np.nan)

    return {
        'laser_energy_eV' : float(components['laser_energy_eV']),
        'broad_eV'        : float(components['broad_eV']),
        'exc_energies_eV' : exc_eV,
        'ph_energies_eV'  : ph_eV,
        'active_modes'    : active.astype(np.int64),
        'I_coh'           : I_coh,
        'I_incoh'         : I_inc,
        'coherence_ratio' : ratio,
        'analysed_modes'  : np.array(sel, dtype=np.int64),
        'group_edges_eV'  : edges,
        'per_mode'        : per_mode,
        'mode_groups'     : groups,
        'per_group'       : per_group,
    }


def save_raman_interference(interference, out_dir='raman_interference'):
    """
    Write the result of `exc_raman_interference_oneph` to disk.

    Files
    -----
    interference.npz                      everything (per-mode keys suffixed _mode_NNN)
    summary_interference.dat              one row per mode: I_coh, I_incoh, ratio
    phase_weight_mode_NNN.dat             (nexc x nexc) maps, gnuplot image format
    cumulative_mode_NNN.dat               running intensities vs number of pairs
    blocks_mode_NNN.dat                   energy-window block decomposition
    summary_groups.dat                    one row per degenerate-mode group
    phase_weight_group_NNN_MMM.dat        } same files for each group in
    cumulative_group_NNN_MMM.dat          } `mode_groups` (modes added
    blocks_group_NNN_MMM.dat              } incoherently)
    README.txt
    """
    import os

    os.makedirs(out_dir, exist_ok=True)
    wL     = float(interference['laser_energy_eV'])
    exc_eV = np.asarray(interference['exc_energies_eV'])
    ph_eV  = np.asarray(interference['ph_energies_eV'])
    active = np.asarray(interference['active_modes']).astype(bool)
    edges  = interference.get('group_edges_eV')
    per    = interference['per_mode']
    per_g  = interference.get('per_group', {})
    nexc   = exc_eV.size
    analysed = set(int(m) for m in interference['analysed_modes'])

    # ---- npz ----
    flat = {k: np.asarray(v) for k, v in interference.items()
            if k not in ('per_mode', 'per_group', 'mode_groups', 'group_edges_eV')}
    if edges is not None:
        flat['group_edges_eV'] = np.asarray(edges)
    for m, d in per.items():
        for k, v in d.items():
            flat['%s_mode_%03d' % (k, m)] = np.asarray(v)
    for label, d in per_g.items():
        for k, v in d.items():
            flat['%s_group_%s' % (k, label)] = np.asarray(v)
    np.savez_compressed(os.path.join(out_dir, 'interference.npz'), **flat)

    # ---- summary ----
    with open(os.path.join(out_dir, 'summary_interference.dat'), 'w') as f:
        f.write('# Exciton-pair interference at laser_energy = %.6f eV\n' % wL)
        f.write('# I_coh = sum_ab |R|^2,  I_incoh = sum_pairs sum_ab |c|^2\n')
        f.write('# ratio < 1: destructive, > 1: constructive, nan: inactive mode\n')
        f.write('# 1:mode_idx  2:ph_freq_cm-1  3:ph_freq_eV  4:active  5:I_coh  '
                '6:I_incoh  7:I_coh/I_incoh  8:analysed\n')
        for m in range(ph_eV.size):
            f.write('%6d  %12.4f  %14.6e  %d  %14.6e  %14.6e  %12.5f  %d\n'
                    % (m, ph_eV[m] * _EV_TO_CM1, ph_eV[m], int(active[m]),
                       interference['I_coh'][m], interference['I_incoh'][m],
                       interference['coherence_ratio'][m], int(m in analysed)))

    # ---- summary of degenerate-mode groups ----
    if per_g:
        with open(os.path.join(out_dir, 'summary_groups.dat'), 'w') as f:
            f.write('# Degenerate-mode groups at laser_energy = %.6f eV\n' % wL)
            f.write('# modes of a group are different final states and add '
                    'incoherently: I = sum_m I_m\n')
            f.write('# ratio < 1: destructive, > 1: constructive\n')
            f.write('# 1:group_idx  2:ph_freq_cm-1_mean  3:ph_freq_eV_mean  4:n_modes  '
                    '5:I_coh  6:I_incoh  7:I_coh/I_incoh  8:modes\n')
            for n, (label, d) in enumerate(per_g.items()):
                ms = [int(m) for m in d['modes']]
                f.write('%6d  %12.4f  %14.6e  %4d  %14.6e  %14.6e  %12.5f  %s\n'
                        % (n, ph_eV[ms].mean() * _EV_TO_CM1, ph_eV[ms].mean(), len(ms),
                           d['I_coh'], d['I_incoh'], d['coherence_ratio'],
                           ','.join(str(m) for m in ms)))

    lam = np.arange(nexc)
    jobs = [('mode_%03d' % m, d,
             'mode %d (omega = %.2f cm-1, %.4f eV), laser_energy = %.6f eV'
             % (m, ph_eV[m] * _EV_TO_CM1, ph_eV[m], wL),
             interference['I_coh'][m], interference['I_incoh'][m])
            for m, d in per.items()]
    for label, d in per_g.items():
        ms = [int(m) for m in d['modes']]
        jobs.append(('group_%s' % label, d,
                     'modes %s summed incoherently (omega = %s cm-1), laser_energy = %.6f eV'
                     % ('+'.join(str(m) for m in ms),
                        ', '.join('%.2f' % (ph_eV[m] * _EV_TO_CM1) for m in ms), wL),
                     d['I_coh'], d['I_incoh']))

    for tag, d, title, Ic, Ii in jobs:
        # ---- phase-weight maps (vectorised: one savetxt per l1 row) ----
        w, sh, it = d['phase_weight'], d['incoherent_share'], d['interference']
        with open(os.path.join(out_dir, 'phase_weight_%s.dat' % tag), 'w') as f:
            f.write('# Pair phase-weight heatmap, %s\n' % title)
            f.write('# phase_weight     = Re(sum_ab conj(R) c) / I_coh   (sums to 1; <0 cancels)\n')
            f.write('# incoherent_share = sum_ab |c|^2 / I_incoh        (sums to 1)\n')
            f.write('# interference     = phase_weight - sum_ab |c|^2 / I_coh\n')
            f.write('#                    (half the cross terms / I_coh; >0 constructive, '
                    '<0 destructive; sums to 1 - I_incoh/I_coh)\n')
            f.write('# 1:l1  2:l2  3:phase_weight  4:incoherent_share  5:interference  '
                    '6:E_l1_eV  7:E_l2_eV\n')
            for l1 in range(nexc):
                blk = np.column_stack([np.full(nexc, l1), lam, w[l1], sh[l1],
                                       it[l1], np.full(nexc, exc_eV[l1]), exc_eV])
                np.savetxt(f, blk, fmt='%6d  %6d  %14.6e  %14.6e  %14.6e  %14.6f  %14.6f')
                f.write('\n')

        # ---- cumulative curve ----
        with open(os.path.join(out_dir, 'cumulative_%s.dat' % tag), 'w') as f:
            f.write('# Cumulative intensity adding pairs strongest-first, %s\n' % title)
            f.write('# complete = %s  (False: curve_max_pairs truncated it)\n'
                    % bool(d['curve_complete']))
            f.write('# final I_coh = %.6e   final I_incoh = %.6e\n' % (Ic, Ii))
            f.write('# 1:n_pairs  2:I_coh_running  3:I_incoh_running  4:I_coh_running/I_coh\n')
            frac = d['curve_I_coh'] / Ic if Ic > 0 else np.full(d['curve_I_coh'].size, np.nan)
            np.savetxt(f, np.column_stack([d['curve_n_pairs'], d['curve_I_coh'],
                                           d['curve_I_incoh'], frac]),
                       fmt='%12d  %14.6e  %14.6e  %12.6f')

        # ---- blocks ----
        if 'block_intensity' in d and edges is not None:
            ng = len(edges) - 1
            with open(os.path.join(out_dir, 'blocks_%s.dat' % tag), 'w') as f:
                f.write('# Energy-window block decomposition, %s\n' % title)
                f.write('# windows are half-open [E_lo, E_hi); i = first-vertex window, '
                        'j = second-vertex window\n')
                f.write('# block_coherence = I_coh / sum_ij intensity_alone = %.6f\n'
                        % d['block_coherence'])
                f.write('# residual_share (states outside all windows) = %.6e\n'
                        % d['residual_share'])
                f.write('# 1:i  2:j  3:E_lo_i  4:E_hi_i  5:E_lo_j  6:E_hi_j  '
                        '7:intensity_alone  8:share  9:interference\n')
                for i in range(ng):
                    for j in range(ng):
                        f.write('%4d  %4d  %10.5f  %10.5f  %10.5f  %10.5f  '
                                '%14.6e  %14.6e  %14.6e\n'
                                % (i, j, edges[i], edges[i + 1], edges[j], edges[j + 1],
                                   d['block_intensity'][i, j], d['block_share'][i, j],
                                   d['block_interference'][i, j]))
                    f.write('\n')

    with open(os.path.join(out_dir, 'README.txt'), 'w') as f:
        f.write(
            'Exciton-pair interference dump (exc_raman_interference_oneph)\n\n'
            'Laser energy : %.6f eV   Nexc : %d\n\n'
            'Each exciton pair (l1, l2) adds a complex amplitude c to the Raman\n'
            'tensor R = sum c. Intensity is |sum c|^2 (coherent); the pair heatmaps\n'
            'show sum |c|^2 (incoherent). Their ratio measures interference.\n\n'
            'summary_interference.dat  I_coh, I_incoh and I_coh/I_incoh per mode\n'
            'phase_weight_mode_NNN.dat signed share of the intensity per pair\n'
            '                          (col 3, sums to 1), no-interference share\n'
            '                          (col 4, sums to 1) and the interference\n'
            '                          (col 5) = col 3 - |c|^2/I_coh: half the\n'
            '                          cross terms, >0 constructive, <0 destructive\n'
            'cumulative_mode_NNN.dat   running I_coh / I_incoh as pairs are added\n'
            '                          strongest first; where col 2 and col 3\n'
            '                          separate, interference sets in\n'
            'blocks_mode_NNN.dat       same decomposition coarse-grained into\n'
            '                          exciton-energy windows\n'
            'summary_groups.dat        I_coh, I_incoh and ratio per group of\n'
            '                          degenerate phonon modes (mode_groups)\n'
            '*_group_NNN_MMM.dat       same files for a group: its modes are\n'
            '                          different final states and are added\n'
            '                          incoherently (sum_m I_m); independent of\n'
            '                          the eigenvector basis of the degenerate\n'
            '                          subspace, unlike the single-mode files\n'
            'interference.npz          all arrays\n\n'
            'l1 = first-vertex exciton, l2 = second-vertex exciton (0-based).\n'
            'gnuplot: splot "phase_weight_mode_NNN.dat" u 1:2:3 w image\n'
            '         (use a diverging palette centred on 0)\n'
            % (wL, nexc))
