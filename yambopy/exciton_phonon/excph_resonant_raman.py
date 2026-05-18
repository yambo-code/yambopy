##
# Authors: PMI, MN
##

import numpy as np
from tqdm import tqdm
from yambopy.units import ha2ev


#@citation("PHYSICAL REVIEW B 113, 085201 (2026)")
def ip_resonant_raman_oneph(laser_energies, ph_energies, el_energies,
                            elec_dipoles, eph_g, cell_vol,
                            broad=0.1, ph_freq_threshold=5.0):
    """
    ! 1-phonon Raman tensor at independent-particle (IP) level
    ! Phys. Rev. B 113, 085201 – Published 2 February, 2026
    ! DOI: 10.1103/ty8m-mgml

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
