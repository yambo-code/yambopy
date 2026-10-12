##
# Authors: (adapted from excph_luminescence.py by MN FP)
##

import numpy as np
from yambopy.units import ha2ev
from yambopy.tools.funcs import bose
from tqdm import tqdm
from yambopy.tools.citations import citation

@citation("P. Lechifflart et al. Phys. Rev. Materials 7, 024006 (2023)")
def exc_ph_absorption(ph_temp,ph_energies,exc_energies,exc_dipoles,exc_ph_mat_el,cell_volume,\
                      exc_energies_in=None,ph_channels='b',nexc_out='all',nexc_in='all',\
                      emin=0,emax=10,estep=0.01,broad=0.025,broad_0=None,\
                      spin_degen=2,pol_dir=None,off_diagonal=True,include_direct=True,\
                      mode_resolved=False,renorm_direct=True,renorm_satellites=True,\
                      return_renorm=False,max_chunk_elements=20_000_000):
    """
    This function calculates the phonon-assisted optical absorption (phonon satellites
    in the macroscopic dielectric function) of indirect materials, i.e. the first-order
    correction to the BSE dielectric function due to exciton-phonon coupling.

    It is the absorption counterpart of `exc_ph_luminescence` and uses the same input data.

    Implemented equations (atomic units, Omega = cell volume, N_q = number of q-points):

      eps_BSE(w) = 1 - 4*pi*s/Omega * sum_l |T_l|^2 / (w - E_{0,l} + i*eta)                    (Eq. 1)

      D eps(w)   = - 4*pi*s/(Omega*N_q) * sum_{q,mu,b,+-}  | sum_l T_l^* D^{+-}_{bl,q mu} |^2
                                       * F^{+-}_{q mu} / ( w - [E_{q,b} +- w_{q mu}] + i*eta )   (Eq. 2)

      D^{+-}_{bl,q mu} = G^{q mu}_{bl} / ( E_{q,b} - E_{0,l} +- w_{q mu} + i*delta )           (Eq. 3)
      F^{+-}_{q mu}    = 1/2 +- 1/2 + N_{q mu}(T)

    '+' = phonon emission (satellite above E_{q,b}), '-' = phonon absorption (below).

    Renormalisation (spectral weight transferred from the exciton peaks to the satellites):

      R_l = 1/N_q * sum_{q,mu,b,+-} F^{+-}_{q mu} |D^{+-}_{bl,q mu}|^2 ,      Z_l = 1 - R_l

      * renorm_direct=True     : direct term uses Z_l |T_l|^2 / (w - E_{0,l} + i*eta)
                                 (as (1-Renorm-Renorm_PH_abs) in Yambo EXCPH_X_phass)
      * renorm_satellites=True : the intermediate exciton l carries weight Z_l also in the
                                 satellites, i.e. T_l -> sqrt(Z_l) T_l inside |sum_l T_l^* D_bl|^2
      R_l is diagonal in l and summed over the phonon channels selected with ph_channels.
      If R_l > 1 (perturbation theory breaking down) a warning is printed and Z_l is set to 0.

    Notes on the prefactor of Eq. 2:
      * Eq. 2 must carry the same -4*pi*s/Omega factor as Eq. 1 (otherwise Im[D eps] < 0
        and the two terms are in different units) plus 1/N_q from the q-average.
      * exc_dipoles from yambopy.bse.excitondipoles.exc_dipoles_pol already contain the
        1/sqrt(N_k) normalisation, so no extra 1/N_k appears in Eq. 1.
      * exc_dipoles are *emission* dipoles T_l = <0|r|l>; the absorption vertex is <l|r|0> = T_l^*.
        exc_ph_mat_el[q,mu,l,b] = <b,q|dV_{q mu}|l,0> (as computed by exciton_phonon_matelem),
        so the absorption amplitude is sum_l G_{bl} T_l^* / (...).

    Returns
    ----------
    * Energy axis in specified range (in eV)
    * Total dielectric function eps(w) (complex): 1 + direct BSE + phonon satellites
      (if include_direct=False: 1 + phonon satellites only)
    * Phonon-assisted correction D eps(w) (complex)
      - shape [nfreqs] or, if mode_resolved=True, [nfreqs,nmodes]
    * if return_renorm=True: R_l [nchannels,nexc_in] (row order: emission, absorption
      as selected in ph_channels)

    Parameters
    ----------
    ph_temp : float
        Lattice temperature in kelvin
    ph_energies : float ndarray
        Phonon energies in eV [nqpts,nmodes]
    exc_energies : float ndarray
        Exciton energies in eV [nqpts,nexc_out] (final states at momentum q)
    exc_dipoles : cmplx ndarray
        Exciton (emission) dipoles in a.u. (bohr), as from exc_dipoles_pol
        * shape [3,nexc_in]: cartesian dipoles (polarisation average or projection, see pol_dir)
        * shape [nexc_in]: dipoles already projected along the field
    exc_ph_mat_el : cmplx ndarray
        Exciton-phonon matrix elements in hartree [nqpts,nmodes,nexc_in,nexc_out]
    cell_volume : float
        Unit cell volume in bohr^3 (e.g. YamboLatticeDB.lat_vol)
    exc_energies_in : float ndarray, optional
        Q=0 exciton energies in eV [nexc_in]. Default is exc_energies[0].
    ph_channels : str
        'e' -> only phonon emission, 'a' -> only phonon absorption, 'b' -> both (default)
    nexc_out, nexc_in : int
        Number of finite-q / optical excitons in the sums. Default: all available states.
    emin, emax, estep : float
        Photon energy range in eV.
    broad : float
        eta in Eqs. 1-2 (half width at half maximum of the Lorentzians) in eV.
    broad_0 : float
        delta in Eq. 3, in eV. Default is broad. (Paper value: 0.001 eV)
    spin_degen : int
        Spin degeneracy factor s (2 for non-spin-polarised singlets, as in YamboExcitonDB.get_chi)
    pol_dir : array_like(3), optional
        Light polarisation direction (cartesian). If None, average over x,y,z.
    off_diagonal : bool
        True  -> coherent sum |sum_l T_l^* D_bl|^2 (gauge-consistent, includes off-diagonal terms)
        False -> incoherent sum sum_l |T_l|^2 |D_bl|^2 (diagonal approximation, for comparison)
    include_direct : bool
        Add the direct BSE term of Eq. 1 to the returned total eps.
    mode_resolved : bool
        If True, D eps is returned resolved in phonon modes [nfreqs,nmodes].
    renorm_direct : bool
        Multiply the direct BSE term by Z_l = 1 - R_l. Default True.
    renorm_satellites : bool
        Weight the intermediate exciton l in the satellites by Z_l. Default True.
    return_renorm : bool
        If True, also return the renormalisation factors R_l per phonon channel.
    max_chunk_elements : int
        Memory control for the frequency loop.
    """
    # Checks
    assert exc_energies.shape[0]==ph_energies.shape[0], "Q-point mismatch between excitons and phonons"
    Nqpts, nmodes = ph_energies.shape
    if exc_energies_in is None: exc_energies_in = exc_energies[0]
    exc_energies    = np.real(exc_energies)
    exc_energies_in = np.real(exc_energies_in)
    assert ph_channels in ['b', 'e', 'a'], "Allowed values for phonon channels are 'b', 'e', 'a'"
    nexc_out_avail = min(exc_energies.shape[1],exc_ph_mat_el.shape[3])
    nexc_in_avail  = min(len(exc_energies_in),exc_dipoles.shape[-1],exc_ph_mat_el.shape[2])
    if nexc_out =='all': nexc_out = nexc_out_avail
    if nexc_in  =='all': nexc_in  = nexc_in_avail
    assert nexc_out <= nexc_out_avail, "less exciton states than requested (nexc_out)"
    assert nexc_in  <= nexc_in_avail,  "less exciton states than requested (nexc_in)"
    exc_energies    = exc_energies[:,:nexc_out]
    exc_energies_in = exc_energies_in[:nexc_in]
    exc_dipoles     = np.asarray(exc_dipoles)[...,:nexc_in]
    G = exc_ph_mat_el[...,:nexc_in,:nexc_out]   # <b,q|dV|l,0>, NOT conjugated for absorption
    assert G.shape[1]==nmodes, "number of modes mismatch between phonon energies and matrix elements"
    if broad_0 is None: broad_0 = broad
    # Acoustic phonons at q=0 are excluded as exc_ph_mat_el[0,0:3,...]=0 and bose(0)=0

    # Dipoles: absorption vertex <l|r|0> = conj(T_l). Shape -> [npol,nexc_in]
    if exc_dipoles.ndim==1:
        dip_abs = np.conj(exc_dipoles)[None,:]
    elif pol_dir is not None:
        e = np.asarray(pol_dir,dtype=float); e = e/np.linalg.norm(e)
        dip_abs = np.einsum('r,ri->i',e,np.conj(exc_dipoles))[None,:]
    else:
        dip_abs = np.conj(exc_dipoles)              # average over 3 directions
    npol = dip_abs.shape[0]

    # Prefactor (atomic units)
    prefactor = 4.*np.pi*spin_degen/cell_volume

    ph_occ   = bose(ph_energies,ph_temp)            # [nq,nmodes]
    broad_Ha = broad/ha2ev

    channels = []
    if ph_channels in ['e','b']: channels.append(+1)  # phonon emission:   F = N+1, pole E + w_ph
    if ph_channels in ['a','b']: channels.append(-1)  # phonon absorption: F = N,   pole E - w_ph

    def satellite_strengths(s):
        """ D^{s}_{bl,q mu} [nq,nmodes,nexc_in,nexc_out] (dimensionless) and F^{s} [nq,nmodes] """
        sat_den = (exc_energies[:,None,None,:]-exc_energies_in[None,None,:,None] \
                   + s*ph_energies[:,:,None,None] + 1j*broad_0)/ha2ev
        return G/sat_den, ph_occ + (1.+s)/2.        # F = 1/2 +- 1/2 + N

    # First pass: renormalisation factors R_l (one row per channel)
    renorm = np.zeros((len(channels),nexc_in))
    for ic,s in enumerate(channels):
        D, F = satellite_strengths(s)
        renorm[ic] = np.einsum('qm,qmio->i',F,np.abs(D)**2,optimize=True)/Nqpts
    Z = 1. - np.sum(renorm,axis=0)
    if np.any(Z<0.):
        print(f"[WARNING] Renormalisation factor R > 1 for exciton(s) {np.where(Z<0.)[0]}: "
              "perturbation theory breaks down, setting Z=0 for them.")
        Z = np.clip(Z,0.,None)
    print("Renormalisation factors R (rows: "+", ".join('emission' if c>0 else 'absorption' for c in channels)+"):")
    for i in range(nexc_in): print(f"  exc {i+1:4d}: "+"  ".join(f"{r:.5e}" for r in renorm[:,i])+f"   Z = {Z[i]:.5f}")

    # Dipoles entering the satellites
    dip_sat = dip_abs*np.sqrt(Z)[None,:] if renorm_satellites else dip_abs

    # Second pass: satellite weights and poles (frequency independent)
    weights, poles = [], []
    for s in channels:
        D, F = satellite_strengths(s)
        if off_diagonal:
            A = np.einsum('pi,qmio->pqmo',dip_sat,D,optimize=True)
            S = np.sum(np.abs(A)**2,axis=0)/npol
        else:
            S = np.einsum('pi,qmio->qmo',np.abs(dip_sat)**2,np.abs(D)**2,optimize=True)/npol
        S = S*F[:,:,None]                           # [nq,nmodes,nexc_out] in bohr^2
        weights.append(S)
        poles.append((exc_energies[:,None,:]+s*ph_energies[:,:,None])/ha2ev)
    weights = np.stack(weights)                     # [nch,nq,nmodes,nexc_out]
    poles   = np.stack(poles)
    # Flatten keeping modes as last axis -> [N,nmodes] one-hot weights
    weights = np.moveaxis(weights,2,-1).reshape(-1,nmodes)
    poles   = np.moveaxis(poles,2,-1).reshape(-1,nmodes)
    N = weights.shape[0]

    # Photon energies
    light_energies = np.arange(emin,emax,estep)
    nfreqs = len(light_energies)
    W = light_energies/ha2ev

    # D eps(w) = - pref/Nq * sum S/(w - pole + i eta)
    deps = np.zeros((nfreqs,nmodes),dtype=complex)
    chunk = max(1,int(max_chunk_elements//max(1,N*nmodes)))
    for i0 in tqdm(range(0,nfreqs,chunk),desc="Ph-assisted absorption"):
        w = W[i0:i0+chunk,None,None]
        deps[i0:i0+chunk] = np.sum(weights[None]/(w-poles[None]+1j*broad_Ha),axis=1)
    deps *= -prefactor/Nqpts

    # Direct BSE term (Eq. 1)
    eps = np.ones(nfreqs,dtype=complex)
    if include_direct:
        T2 = np.sum(np.abs(dip_abs)**2,axis=0)/npol                      # [nexc_in]
        if renorm_direct: T2 = Z*T2
        E0 = exc_energies_in/ha2ev
        eps -= prefactor*np.sum(T2[None,:]/(W[:,None]-E0[None,:]+1j*broad_Ha),axis=1)
    eps += np.sum(deps,axis=1)

    if not mode_resolved: deps = np.sum(deps,axis=1)
    if return_renorm: return light_energies, eps, deps, renorm
    return light_energies, eps, deps
