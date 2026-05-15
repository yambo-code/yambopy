##
# Authors: PMI, MN
##

import numpy as np
from tqdm import tqdm
from yambopy.units import ha2ev


def ip_resonant_raman_oneph(laser_energies, ph_energies, el_energies,
                            elec_dipoles, eph_g, n_val, cell_vol,
                            broad=0.1, ph_freq_threshold=5.0,
                            eph_bands=None):
    """
    Independent-particle one-phonon resonant Raman tensor at q=0 (Stokes).

    Computes the Stokes branch (phonon emission) of the third-order
    light-matter / electron-phonon perturbation expression at the
    independent-particle level. Resonant and anti-resonant time orderings
    are both included; the anti-Stokes channel (phonon absorption) is not.
    Both the electron-scattering and hole-scattering vertices are summed
    with the standard relative minus sign.

    The electron-phonon matrix elements (``eph_g``) may cover a band
    window that is *smaller than* (or otherwise different from) the one
    spanned by ``el_energies`` / ``elec_dipoles``. Use the ``eph_bands``
    argument to specify, for each band axis of ``eph_g``, which band of
    ``el_energies`` it corresponds to. Bands in ``el_energies`` that are
    not listed in ``eph_bands`` contribute to the optical (dipole) sums
    but not to the phonon-mediated scattering.

    Mirrors `compute_Raman_oneph_ip` from
    https://github.com/muralidhar-nalabothula/PhdScripts/blob/main/exph/raman.py
    (which is the special case ``eph_bands = np.arange(nb)``).

    Returns
    -------
    raman_tensor : (nfreqs, nmodes, 3, 3) complex ndarray
        Raman tensor R^nu_{alpha,beta}(omega_L). For incoming polarization
        e_in and outgoing polarization e_out the (un-normalized) intensity
        is |sum_{ab} e_out[a] R[w, m, a, b] e_in[b]|^2.

    Parameters
    ----------
    laser_energies : (nfreqs,) float ndarray
        Incoming laser energies in eV.
    ph_energies : (nmodes,) float ndarray
        Phonon energies at q=0 in eV.
    el_energies : (nk, nb) float ndarray
        Single-particle (KS or QP) energies in eV. The first ``n_val``
        bands are valence, the remaining ``nc = nb - n_val`` are
        conduction (yambo convention).
    elec_dipoles : (3, nk, nc, nv) complex ndarray
        Velocity-gauge electronic dipoles <c k | v | v k> in atomic
        units, as stored in yambo's ``ndb.dipoles``.
    eph_g : (nmodes, nk, nb_eph, nb_eph) complex ndarray
        Electron-phonon matrix elements at q=0 in Hartree, with index
        order (mode, k, final_band, initial_band) and the standard
        ``1/sqrt(2*omega_nu)`` normalization. The matrix is square in
        its last two axes; ``nb_eph`` can be less than or equal to
        ``nb``.
    n_val : int
        Number of valence bands inside the ``el_energies`` window.
    cell_vol : float
        Unit-cell volume in bohr^3.
    broad : float, optional
        Lorentzian broadening Gamma in eV. Default: 0.1.
    ph_freq_threshold : float, optional
        Acoustic-mode cutoff in cm^-1; modes with |omega_nu| below this
        contribute zero Raman tensor. Default: 5.
    eph_bands : array-like of int, optional
        Length ``nb_eph``. ``eph_bands[i]`` is the index (within the
        band axis of ``el_energies``, in [0, nb)) of the band that the
        i-th band axis of ``eph_g`` corresponds to. Indices in
        ``[0, n_val)`` are interpreted as valence, indices in
        ``[n_val, nb)`` as conduction. If ``None`` (default), ``eph_g``
        must cover all ``nb`` bands of ``el_energies`` in their natural
        order.
    """
    # cm^-1 -> Ha (Murali's constant: 1 cm^-1 = 0.12398e-3 eV)
    cm1_to_eV = 0.12398e-3
    cm1_to_Ha = cm1_to_eV / ha2ev

    # ---- Input validation (raise so it survives `python -O`) ---------------
    if el_energies.ndim != 2:
        raise ValueError("el_energies must be (nk, nb)")
    nk, nb = el_energies.shape
    nv = int(n_val)
    nc = nb - nv
    if not (0 < nv < nb):
        raise ValueError("n_val must satisfy 0 < n_val < nb; got n_val=%d, nb=%d"
                         % (nv, nb))
    if tuple(elec_dipoles.shape) != (3, nk, nc, nv):
        raise ValueError("elec_dipoles shape %s != (3, %d, %d, %d)"
                         % (tuple(elec_dipoles.shape), nk, nc, nv))
    nmodes = ph_energies.shape[0]
    if eph_g.ndim != 4 or eph_g.shape[0] != nmodes or eph_g.shape[1] != nk:
        raise ValueError("eph_g shape %s incompatible with (nmodes=%d, nk=%d, nb_eph, nb_eph)"
                         % (tuple(eph_g.shape), nmodes, nk))
    if eph_g.shape[2] != eph_g.shape[3]:
        raise ValueError("eph_g must be square in its last two axes; got %s"
                         % (tuple(eph_g.shape),))
    nb_eph = eph_g.shape[2]

    if eph_bands is None:
        if nb_eph != nb:
            raise ValueError(
                "eph_g spans nb_eph=%d bands but el_energies spans nb=%d; "
                "pass `eph_bands` (length nb_eph) to specify the mapping."
                % (nb_eph, nb))
        eph_bands = np.arange(nb)
    else:
        eph_bands = np.asarray(eph_bands, dtype=int)
        if eph_bands.shape != (nb_eph,):
            raise ValueError("eph_bands must have shape (%d,); got %s"
                             % (nb_eph, tuple(eph_bands.shape)))
        if (eph_bands < 0).any() or (eph_bands >= nb).any():
            raise ValueError("eph_bands must contain band indices in [0, %d); got min=%d, max=%d"
                             % (nb, int(eph_bands.min()), int(eph_bands.max())))
        if len(np.unique(eph_bands)) != nb_eph:
            raise ValueError("eph_bands must not contain duplicate band indices")

    if laser_energies.ndim != 1:
        raise ValueError("laser_energies must be 1D")
    if cell_vol <= 0.:
        raise ValueError("cell_vol must be positive")

    # ---- Unit conversion ---------------------------------------------------
    laser_Ha     = np.asarray(laser_energies) / ha2ev
    ph_Ha        = np.asarray(ph_energies)    / ha2ev
    El_Ha        = np.asarray(el_energies)    / ha2ev
    # Murali's convention: half of input broadening enters the denominator
    broad_Ha     = (broad / ha2ev) / 2.0
    ph_thresh_Ha = ph_freq_threshold * cm1_to_Ha

    # Vertical c-v gaps minus i*Gamma : delta[k, c, v] = eps_c(k) - eps_v(k) - i*Gamma
    delta = (El_Ha[:, nv:, None] - El_Ha[:, None, :nv]) - 1j * broad_Ha   # (nk, nc, nv)

    # d_abs = d* (incoming-photon vertex), d_emi = d (outgoing-photon vertex)
    d_abs = elec_dipoles.conj()
    d_emi = elec_dipoles

    # ---- Identify elph valence / conduction subsets ------------------------
    v_mask      = (eph_bands < nv)                        # (nb_eph,) bool
    c_mask      = ~v_mask
    eph_v_local = eph_bands[v_mask]                       # indices in [0, nv)
    eph_c_local = eph_bands[c_mask] - nv                  # indices in [0, nc)
    n_eph_v     = int(v_mask.sum())
    n_eph_c     = int(c_mask.sum())

    # cc and vv blocks of eph_g, with the same ordering as eph_c_local /
    # eph_v_local (so the contractions below stay consistent).
    gcc = eph_g[:, :, c_mask, :][:, :, :, c_mask]         # (nmodes, nk, n_eph_c, n_eph_c)
    gvv = eph_g[:, :, v_mask, :][:, :, :, v_mask]         # (nmodes, nk, n_eph_v, n_eph_v)

    # ---- Sliced dipoles / gaps per channel ---------------------------------
    # Electron channel: c restricted to elph_c, v spans full dipole-v.
    if n_eph_c > 0:
        delta_e = delta[:, eph_c_local, :]                # (nk, n_eph_c, nv)
        d_abs_e = d_abs[:, :, eph_c_local, :]             # (3, nk, n_eph_c, nv)
        d_emi_e = d_emi[:, :, eph_c_local, :]
    # Hole channel: c spans full dipole-c, v restricted to elph_v.
    if n_eph_v > 0:
        delta_h = delta[:, :, eph_v_local]                # (nk, nc, n_eph_v)
        d_abs_h = d_abs[:, :, :, eph_v_local]             # (3, nk, nc, n_eph_v)
        d_emi_h = d_emi[:, :, :, eph_v_local]

    # ---- Output ------------------------------------------------------------
    nfreqs  = laser_Ha.size
    Ram     = np.zeros((nfreqs, nmodes, 3, 3), dtype=np.complex128)
    ram_fac = 1.0 / nk / np.sqrt(cell_vol)

    for w_idx in tqdm(range(nfreqs), desc="IP one-phonon Raman"):
        wL = laser_Ha[w_idx]

        # First-vertex dressed dipoles (no phonon shift in the denominator)
        if n_eph_c > 0:
            dipS_res_e  = d_abs_e / (wL - delta_e)[None, ...]   # (3, nk, n_eph_c, nv)
            dipS_ares_e = d_emi_e / (wL + delta_e)[None, ...]
        if n_eph_v > 0:
            dipS_res_h  = d_abs_h / (wL - delta_h)[None, ...]   # (3, nk, nc, n_eph_v)
            dipS_ares_h = d_emi_h / (wL + delta_h)[None, ...]

        for m in range(nmodes):
            wph = ph_Ha[m]
            if abs(wph) <= ph_thresh_Ha:
                continue

            R_e = 0.0
            R_h = 0.0

            # ---- Electron channel (g_cc) -------------------------------
            if n_eph_c > 0:
                dipSp_res_e  = d_emi_e / (wL - delta_e - wph)[None, ...]
                dipSp_ares_e = d_abs_e / (wL + delta_e - wph)[None, ...]
                gcc_m = gcc[m]                                   # (nk, n_eph_c, n_eph_c)
                # tmp[a, k, c_final, v] = sum_{c_init} g[k, c_final, c_init].conj
                #                         * dipSp[a, k, c_init, v]
                tmp_e_res  = np.einsum('kij,lkjv->lkiv', gcc_m.conj(), dipSp_res_e)
                tmp_e_ares = np.einsum('kij,lkjv->lkiv', gcc_m,        dipSp_ares_e)
                # Outer contraction over (k, c_final, v).
                R_e = (np.einsum('akcv,bkcv->ab', dipS_res_e,  tmp_e_res)
                     + np.einsum('akcv,bkcv->ab', dipS_ares_e, tmp_e_ares))

            # ---- Hole channel (g_vv) -----------------------------------
            if n_eph_v > 0:
                dipSp_res_h  = d_emi_h / (wL - delta_h - wph)[None, ...]
                dipSp_ares_h = d_abs_h / (wL + delta_h - wph)[None, ...]
                gvv_m = gvv[m]                                   # (nk, n_eph_v, n_eph_v)
                # tmp[a, k, c, v_init] = sum_{v_final} dipSp[a, k, c, v_final]
                #                         * g[k, v_final, v_init].conj
                tmp_h_res  = np.einsum('lkcv,kvw->lkcw', dipSp_res_h,  gvv_m.conj())
                tmp_h_ares = np.einsum('lkcv,kvw->lkcw', dipSp_ares_h, gvv_m)
                R_h = (np.einsum('akcv,bkcv->ab', dipS_res_h,  tmp_h_res)
                     + np.einsum('akcv,bkcv->ab', dipS_ares_h, tmp_h_ares))

            # Electron minus hole, with outgoing-photon radiation factor
            # and BZ / volume normalisation.
            Ram[w_idx, m] = (R_e - R_h) * np.sqrt(abs(wL - wph) / wL) * ram_fac

    return Ram
