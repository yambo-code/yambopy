##
# Authors: PMI, MN
##

import warnings

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


def _per_mode_intensity(R, pol_in=None, pol_out=None):
    """Per-mode Raman intensity from a Raman tensor.

    Parameters
    ----------
    R : (..., nmodes, 3, 3) complex ndarray
        Raman tensor, optionally carrying leading axes (e.g. a
        laser-energy probe window of shape (n_probe, nmodes, 3, 3)).
    pol_in, pol_out : array-like of shape (3,), optional
        Incoming / outgoing photon polarisation 3-vectors (complex
        allowed, normalised internally). If either is None the
        unpolarised Eq. (D1) sum ``Σ_{μν} |R[..,λ,μ,ν]|²`` is used,
        else ``|Σ_{μν} e_out[μ] R[..,λ,μ,ν] e_in[ν]|²``.

    Returns
    -------
    (..., nmodes) real ndarray
        Intensity per mode, preserving any leading axes of ``R``.
    """
    R = np.asarray(R)
    if (pol_in is None) or (pol_out is None):
        return np.sum(np.abs(R)**2, axis=(-1, -2))
    e_in  = np.asarray(pol_in,  dtype=complex)
    e_out = np.asarray(pol_out, dtype=complex)
    e_in  = e_in  / np.linalg.norm(e_in)
    e_out = e_out / np.linalg.norm(e_out)
    amp = np.einsum('m, ...lmn, n -> ...l', e_out, R, e_in)
    return np.abs(amp)**2


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
        Number of laser-energy probe points used for *windowed probing*
        of the resonance profile. Default 1 (single laser energy ω_L,
        identical to the previous behaviour).

        For ``n_probe > 1`` the per-mode intensity is integrated over a
        symmetric window of laser energies centred on ω_L:

            ω_L + k·h ,  k = -(n-1)/2 … +(n-1)/2 ,  h = γ/4 ,  γ = `broad`

        following the windowed-probe recipe (step h tied to the known
        broadening γ; integral via Simpson, which requires an odd number
        of points). This catches a resonance sitting near ω_L even when
        ω_L itself lands on a tail. The per-mode intensity fed to the
        spectrum is then ``∫ I_λ(ω_L) dω_L`` over the window (Simpson),
        a conserved-area collapse.

        Requires ``raman_tensor is None`` (the tensor must be recomputed
        at every probe energy) and ``n_probe`` odd. A warning is issued
        if the implied window reach (``(n-1)/8`` in units of γ) is much
        smaller than the ~10γ recommended for a Lorentzian resonance.

    Returns
    -------
    energy_grid : (nE,) float ndarray
        Energy axis in `energy_units`.
    intensity   : (nE,) float ndarray
        Raman intensity at each energy.
    """
    # ---- Compute or accept the Raman tensor, then per-mode intensity ------
    nmodes  = int(ph_energies.shape[0])
    n_probe = int(n_probe)
    if n_probe < 1:
        raise ValueError("n_probe must be >= 1 (got %d)" % n_probe)

    if raman_tensor is not None:
        if n_probe > 1:
            raise ValueError(
                "windowed probing (n_probe > 1) recomputes the Raman tensor "
                "at several laser energies and is incompatible with a "
                "precomputed `raman_tensor`")
        R = np.asarray(raman_tensor)
        if R.ndim == 4 and R.shape[0] == 1:
            R = R[0]
        if R.shape != (nmodes, 3, 3):
            raise ValueError("raman_tensor shape %s != (%d, 3, 3)" %
                             (tuple(R.shape), nmodes))
        I_mode = _per_mode_intensity(R, pol_in, pol_out)       # (nmodes,)
    else:
        for name, arr in (('el_energies', el_energies),
                          ('elec_dipoles', elec_dipoles),
                          ('eph_g', eph_g),
                          ('cell_vol', cell_vol)):
            if arr is None:
                raise ValueError("`%s` must be provided when `raman_tensor` is None" % name)

        if n_probe == 1:
            laser_grid = np.array([float(laser_energy)])
        else:
            if n_probe % 2 == 0:
                raise ValueError("n_probe must be odd (Simpson integration "
                                 "needs an even number of intervals); got %d"
                                 % n_probe)
            if broad <= 0.:
                raise ValueError("windowed probing needs broad (gamma) > 0")
            h      = broad / 4.0                               # step = gamma/4
            n_half = (n_probe - 1) // 2
            laser_grid = float(laser_energy) + np.arange(-n_half, n_half + 1) * h
            reach_g = n_half * h / broad                       # = n_half/4, in units of gamma
            if reach_g < 10.0:
                warnings.warn(
                    "windowed-probe reach is %.2f*gamma (n_probe=%d, h=gamma/4); "
                    "the Lorentzian resonance tail recommends ~10-30*gamma "
                    "(n_probe >= 81). The window only samples the immediate "
                    "neighbourhood of the laser energy." % (reach_g, n_probe),
                    stacklevel=2)

        R_win = ip_resonant_raman_tensor_oneph(
                    laser_grid, ph_energies, el_energies, elec_dipoles, eph_g,
                    cell_vol, broad=broad, ph_freq_threshold=ph_freq_threshold)
        I_win = _per_mode_intensity(R_win, pol_in, pol_out)    # (n_probe, nmodes)

        if n_probe == 1:
            I_mode = I_win[0]                                  # (nmodes,)
        else:
            from scipy.integrate import simpson
            I_mode = simpson(I_win, dx=h, axis=0)              # (nmodes,) conserved area

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
                               save_per_pair=False):
    """
    Diagnostic version of `exc_resonant_raman_tensor_oneph`. Computes the
    Raman tensor at a SINGLE laser energy in pure numpy and returns every
    intermediate quantity that goes into it.

    The returned dict is intended to be fed to `save_raman_components`,
    which writes it both as a single .npz archive and as a folder of
    column-format .dat files suitable for gnuplot.

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
        Memory ~ nmodes * 9 * nexc^2 * 16 bytes. Off by default.

    Returns
    -------
    dict
        See README produced by `save_raman_components` for the
        complete list of keys.
    """
    cm1_to_Ha = 0.12398e-3 / ha2ev

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

    dipS_res   = D_abs * inv_denom_res_v1[None, :]                            # (3, nexc)
    dipS_ares  = D_emi * inv_denom_ares_v1[None, :]
    dipSp_res  = D_emi[None, :, :] * inv_denom_res_v2[:, None, :]             # (nmodes, 3, nexc)
    dipSp_ares = D_abs[None, :, :] * inv_denom_ares_v2[:, None, :]

    exc_ph_abs = np.asarray(exc_ph_mat_el, dtype=complex)   # yambopy raw (mode, init, fin)
    exc_ph_emi = np.conj(exc_ph_abs)                         # Stokes / emission

    radiation_factor = np.sqrt(np.abs(wL - ph_Ha) / wL)
    active_modes     = np.abs(ph_Ha) > ph_thresh_Ha
    prefactor        = radiation_factor * ram_fac            # (nmodes,) real

    term_res  = np.einsum('al, mLl, mbL -> mab',
                          dipS_res,  exc_ph_emi, dipSp_res,  optimize=True)
    term_ares = np.einsum('al, mLl, mbL -> mab',
                          dipS_ares, exc_ph_abs, dipSp_ares, optimize=True)

    term_res  = term_res  * prefactor[:, None, None]
    term_ares = term_ares * prefactor[:, None, None]
    term_res[~active_modes]  = 0.0
    term_ares[~active_modes] = 0.0
    raman_tensor = term_res + term_ares

    # Per-pair factorized intensity (sum over polarizations); cheap & always returned.
    A_res     = np.sum(np.abs(dipS_res)**2,  axis=0)
    A_ares    = np.sum(np.abs(dipS_ares)**2, axis=0)
    B_res     = np.sum(np.abs(dipSp_res)**2,  axis=1)
    B_ares    = np.sum(np.abs(dipSp_ares)**2, axis=1)
    g_sq_pair = (np.abs(exc_ph_abs)**2).transpose(0, 2, 1)
    pref_sq   = (prefactor**2)[:, None, None]

    pair_intensity_res  = (A_res[None, :, None]  * g_sq_pair * B_res[:, None, :]) * pref_sq
    pair_intensity_ares = (A_ares[None, :, None] * g_sq_pair * B_ares[:, None, :]) * pref_sq
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
        'exc_ph_emi'              : exc_ph_emi,
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
        per_pair_res  = np.einsum('al, mLl, mbL -> mablL',
                                  dipS_res,  exc_ph_emi, dipSp_res,  optimize=True)
        per_pair_ares = np.einsum('al, mLl, mbL -> mablL',
                                  dipS_ares, exc_ph_abs, dipSp_ares, optimize=True)
        per_pair_res  = per_pair_res  * prefactor[:, None, None, None, None]
        per_pair_ares = per_pair_ares * prefactor[:, None, None, None, None]
        per_pair_res[~active_modes]  = 0.0
        per_pair_ares[~active_modes] = 0.0
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
