##
# Authors: MN FP
##

import numpy as np
from yambopy.units import ha2ev
from yambopy.tools.funcs import bose,boltzman_f
from tqdm import tqdm
from numba import njit,prange

@njit(cache=True, nogil=True, parallel=True)
def exc_ph_luminescence(ph_temp,ph_energies,exc_energies,exc_dipoles,exc_ph_mat_el,\
                        exc_energies_in=None,exc_temp=None,ph_channels='b',\
                        PL_energy_prefactor='PT',nexc_out='all',nexc_in='all',\
                        emin=0,emax=10,estep=0.01,broad=0.1,broad_0=None):
    """
    This function calculates the phonon-assisted satellites (phonon replica) in the 
    luminescence of indirect materials (i.e., phonon-mediated exciton recombination process).

    It uses a perturbation theory derivation that includes the off-diagonal elements
    of the exciton-phonon matrix elements.

    Exciton populations are approximated with a Boltzmann distribution function.

    See:
    - Phys. Rev. Materials 7, 024006 (2023)
    - Phys. Rev. Letters 131, 206902 (2023) Equations from (1) to (4)

    Returns:
    ----------
    * Energy axis in specified range (in eV)
    * Luminescence intensity

    Parameters
    ----------
    ph_temp : float
        Lattice temperature in kelvin
    ph_energies : float ndarray
        Phonon energies in eV [nqpts,nmodes]
    exc_energies : float ndarray
        Exciton energies in eV [nqpts,nexc_out] 
        Typically, YamboExcitonDB.eigenvalues are used
    exc_dipoles : cmplx ndarray
        Exciton dipole matrix elements in a.u. (bohr) [nexc_in]
        Typically, YamboExcitonDB.l_residuals are used.
    exc_ph_mat_el : cmplx ndarray
        Exciton-phonon coupling matrix elements in a.u. (hartree) [nqpts,nmodes,nexc_in,nexc_out]
        Typically, computed via exciton_phonon_matelem
    exc_energies_in : float ndarray, optional
        If given, zero-momentum states will be taken from this array. Default is exc_energies.
    exc_temp : float
        If given, effective temperature controlling excitonic population. Default is ph_temp.
    ph_channels : str
        'e'  -> only phonon emission channel
        'a' -> only phonon absorption channel
        'b' -> both channels (default)
    PL_energy_prefactor : str
        'PT' -> prefactor from time-dependent perturbation theory (default)
        'RS' -> prefactor from van Roosbroeck-Shockley relation (not including refraction index).
        None or any other str -> no prefactor is used
    nexc_out : int
        Number of finite-q excitons in the sum. Default: all available states.
    nexc_in : int
        Number of optical excitons in the sum. Default: all available states.
    emin, emax, esteps : float
        Laser energy range parameters in eV.
    broad : float
        Broadening parameter in eV.
    broad_0 : float
        Broadening parameter used inside satellite oscillator strengths in eV. Default is broad.
    """
    def get_PL_satellite(W,ph_sign=-1):
        """ Actual calculation of PL satellites at each frequency.
            Uses the same variables as main function.

            W --> laser energy value in hartree
            ph_sign -->
                * -1 : phonon emission
                * +1 : phonon absorption
        """
        # Satellite oscillator strenghts
        ## Dimensions: [nqpts, nmodes, nexc_in, nexc_out]
        satellite_energy = exc_energies[:,None,None,:]-exc_energies_in[None,None,:,None] \
                            + ph_sign*ph_energies[:,:,None,None] + 1j*broad_0
        satellite_energy_Ha = satellite_energy/ha2ev
        satellite_weight = exc_ph_mat_el/satellite_energy_Ha
    
        # Phonon emission/absorption factor: if em --> ph_occ+1, if abs --> ph_occ
        ph_fac = ph_occ - (ph_sign-1)/2.
    
        # Sum over exc_in and then square (this includes off-diagonal excph SE terms)
        T = np.einsum('i,qmio->qmo',exc_dipoles,satellite_weight,optimize=True)
        # Include occupation factors
        T = np.abs(T)**2 * exc_occ[:,None,:]*ph_fac[:,:,None]

        # Energy differences [nqpts,nmodes,nexc_out]
        pole_energy = exc_energies[:,None,:]+ph_sign*ph_energies[:,:,None]
        pole_energy_Ha = pole_energy[:,:,:]/ha2ev
        
        # PL energy prefactors
        # Note: with 'PT+RS' as argument both will be considered
        light_fac = 1.
        if 'PT' in PL_energy_prefactor: # [nqpts,nmodes,nexc_out]
            light_fac = light_fac * 1./pole_energy_Ha
        if 'RS' in PL_energy_prefactor: # [nqpts,nmodes,nexc_out]
            light_fac = light_fac * W * (W-ph_sign*ph_energies[:,:,None])**2.

        # Energy conservation
        delta_funct = 1./((W-pole_energy_Ha)**2.+ broad_Ha**2.)
        E = light_fac*delta_funct
        # Dimensions
        E = E * broad_Ha/np.pi/Nqpts * ha2ev**2

        # Putting everything together
        return np.einsum('qmo,qmo->', E, T, optimize=True) 

    # Checks
    assert exc_energies.shape[0]==ph_energies.shape[0], "Q-point mismatch between excitons and phonons"
    Nqpts = ph_energies.shape[0]
    if exc_energies_in is None: exc_energies_in = exc_energies[0]
    if np.iscomplexobj(exc_energies):    exc_energies    = exc_energies.real
    if np.iscomplexobj(exc_energies_in): exc_energies_in = exc_energies_in.real
    assert ph_channels in ['b', 'e', 'a'], "Allowed values for phonon channels are 'b', 'e', 'a'"
    nexc_out_avail = min(exc_energies.shape[1],exc_ph_mat_el.shape[3])
    nexc_in_avail  = min(len(exc_energies_in),len(exc_dipoles),exc_ph_mat_el.shape[2])
    assert nexc_out <= nexc_out_avail, "less exciton states than requested (nexc_out)"
    assert nexc_in  <= nexc_in_avail,  "less exciton states than requested (nexc_in)"
    if nexc_out =='all': nexc_out = nexc_out_avail
    if nexc_in  =='all': nexc_in  = nexc_in_avail
    exc_energies    = exc_energies[:,:nexc_out]
    exc_energies_in = exc_energies_in[:nexc_in]
    # We conj because we need it for photon emission
    exc_ph_mat_el = np.conj( exc_ph_mat_el[...,:nexc_in,:nexc_out] ) 
    assert ph_energies.shape[1]==exc_ph_mat_el.shape[1], "number of modes mismatch between phonon energies and matrix elements"
    broad = broad/2
    if broad_0 is None: broad_0 = broad
    if exc_temp is None: exc_temp = ph_temp
    
    # Occupation functions
    exc_min_energy = np.min(exc_energies)
    exc_occ = boltzman_f(exc_energies-exc_min_energy,exc_temp)
    ph_occ  = bose(ph_energies,ph_temp)

    # Laser energies
    light_energies = np.arange(emin,emax,estep,dtype=np.float32)
    nfreqs = len(light_energies)
    light_energies_Ha = light_energies/ha2ev#light_energies[None,None,None,:]/ha2ev
    broad_Ha = broad/ha2ev

    # Calculation
    PL_satellites = np.zeros(nfreqs)
    for w in tqdm(prange(nfreqs)):
        # Accumulate phonon emission satellites
        if ph_channels=='e' or ph_channels=='b': 
            PL_satellites[w] += get_PL_satellite(light_energies_Ha[w],ph_sign=-1)
        # Accumulate phonon absorption satellites
        if ph_channels=='a' or ph_channels=='b': 
            PL_satellites[w] += get_PL_satellite(light_energies_Ha[w],ph_sign=+1)

    return light_energies, PL_satellites