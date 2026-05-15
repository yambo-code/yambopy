##
# Authors: PMI, MN
##

import numpy as np
from tqdm import tqdm
from yambopy.units import ha2ev


def ip_resonant_raman_oneph(laser_energies, ph_energies, el_energies,
                            elec_dipoles, eph_g, n_val, cell_vol,
                            broad=0.1, ph_freq_threshold=5.0):
    """
    Independent-particle one-phonon resonant Raman tensor at q=0 (Stokes).

    Computes the Stokes branch (phonon emission) of the third-order
    light-matter / electron-phonon perturbation expression at the
    independent-particle level. Resonant and anti-resonant time orderings
    are both included; the anti-Stokes channel (phonon absorption) is not.
    Both the electron-scattering and hole-scattering vertices are summed
    with the standard relative minus sign.

    Mirrors `compute_Raman_oneph_ip` from
    https://github.com/muralidhar-nalabothula/PhdScripts/blob/main/exph/raman.py

    The function takes all physical ingredients as input arrays and only
    evaluates the formula; it does not load yambo databases.

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
        Single-particle (KS or QP) energies in eV. The first ``n_val`` bands
        are valence, the remaining ``nc = nb - n_val`` are conduction
        (yambo convention).
    elec_dipoles : (3, nk, nc, nv) complex ndarray
        Velocity-gauge electronic dipoles <c k | v | v k> in atomic units,
        as stored in yambo's ``ndb.dipoles``. The conduction index runs over
        the ``nc`` states above ``n_val``; the valence index over the
        ``n_val`` states below.
    eph_g : (nmodes, nk, nb, nb) complex ndarray
        Electron-phonon matrix elements at q=0 in Hartree, with index order
        (mode, k, final_band, initial_band) and the standard
        ``1/sqrt(2*omega_nu)`` normalization (yambo / LetzElPhC convention).
    n_val : int
        Number of valence bands.
    cell_vol : float
        Unit-cell volume in bohr^3.
    broad : float, optional
        Lorentzian broadening Gamma in eV. Default: 0.1.
    ph_freq_threshold : float, optional
        Acoustic-mode cutoff in cm^-1; modes with |omega_nu| below this
        contribute zero Raman tensor. Default: 5.
    """
    # cm^-1 -> Ha (Murali's constant: 1 cm^-1 = 0.12398e-3 eV)
    cm1_to_eV = 0.12398e-3
    cm1_to_Ha = cm1_to_eV / ha2ev

    # Shape checks
    #assert el_energies.ndim == 2, "el_energies must be (nk, nb)"
    #nk, nb = el_energies.shape
    #nv = int(n_val)
    #nc = nb - nv
    #assert 0 < nv < nb, "n_val must satisfy 0 < n_val < nb"
    #assert elec_dipoles.shape == (3, nk, nc, nv), \
    #    f"elec_dipoles shape {elec_dipoles.shape} != (3, {nk}, {nc}, {nv})"
    #nmodes = ph_energies.shape[0]
    #assert eph_g.shape == (nmodes, nk, nb, nb), \
    #    f"eph_g shape {eph_g.shape} != ({nmodes}, {nk}, {nb}, {nb})"
    #assert laser_energies.ndim == 1
    #assert cell_vol > 0.

    # Unit conversion to Hartree
    #laser_Ha = np.asarray(laser_energies) / ha2ev
    #ph_Ha    = np.asarray(ph_energies)    / ha2ev
    El_Ha    = np.asarray(el_energies)    / ha2ev
    # Murali's convention: half of input broadening enters the denominator
    broad_Ha = (broad / ha2ev) / 2.0
    ph_thresh_Ha = ph_freq_threshold * cm1_to_Ha

    #TOCHECK sign on broadening 
    # Vertical c-v gaps minus i*Gamma : delta[k, c, v] = eps_c(k) - eps_v(k) - i*Gamma
    delta = (El_Ha[:, nv:, None] - El_Ha[:, None, :nv]) - 1j * broad_Ha   # (nk, nc, nv)

    # d_abs = d* (incoming-photon vertex), d_emi = d (outgoing-photon vertex)
    d_abs = elec_dipoles.conj()
    d_emi = elec_dipoles

    # c-c and v-v blocks of the e-ph matrix elements: index order [final, initial]
    gcc = eph_g[:, :, nv:, nv:]   # (nmodes, nk, nc, nc)
    gvv = eph_g[:, :, :nv, :nv]   # (nmodes, nk, nv, nv)

    # Output
    nfreqs = laser_Ha.size
    Ram = np.zeros((nfreqs, nmodes, 3, 3), dtype=np.complex128)

    # Murali's overall prefactor: 1 / N_k / sqrt(V)
    ram_fac = 1.0 / nk / np.sqrt(cell_vol)

    for w_idx in tqdm(range(nfreqs), desc="IP one-phonon Raman"):
        wL = laser_Ha[w_idx]

        # First-vertex dressed dipoles (no phonon shift in denominator)
        dipS_res  = d_abs / (wL - delta)[None, ...]   # (3, nk, nc, nv)
        dipS_ares = d_emi / (wL + delta)[None, ...]

        for m in range(nmodes):
            wph = ph_Ha[m]
            if abs(wph) <= ph_thresh_Ha:
                continue

            # Second-vertex dressed dipoles (Stokes => -wph in the denominator)
            dipSp_res  = d_emi / (wL - delta - wph)[None, ...]
            dipSp_ares = d_abs / (wL + delta - wph)[None, ...]

            gcc_m = gcc[m][None, ...]   # (1, nk, nc, nc), broadcasts over polarization
            gvv_m = gvv[m][None, ...]   # (1, nk, nv, nv)

            # Electron-channel minus hole-channel (matmul on the last two axes)
            tmp_res  = gcc_m.conj() @ dipSp_res  - dipSp_res  @ gvv_m.conj()
            tmp_ares = gcc_m        @ dipSp_ares - dipSp_ares @ gvv_m

            # Sum over (k, c, v) -> Raman tensor element (alpha, beta)
            Ram[w_idx, m] = np.einsum('akcv,bkcv->ab', dipS_res,  tmp_res ) \
                          + np.einsum('akcv,bkcv->ab', dipS_ares, tmp_ares)

            # Outgoing-photon radiation factor and BZ/volume normalization
            Ram[w_idx, m] *= np.sqrt(abs(wL - wph) / wL) * ram_fac

    return Ram
