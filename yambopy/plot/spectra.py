#
# License-Identifier: GPL
#
# Copyright (C) 2024 The Yambo Team
#
# Authors: FP
#
# This file is part of the yambopy project
#
import numpy as np

def get_spectra(energies,constant=1.,weights=None,residuals=None,broadening=0.01,emin=0,emax=10,estep=0.01):
    """
    DOS and spectral functions.

    This class computes the generic quantity: 

            D(w)= C sum\_i^N_i wa\_i |D\_i|^2 delta_b(w - dE\_i)
    
    Inputs:
    * dE_i:    pole energies in eV, must be an array of shape (N_momenta, N_states) or (N_states)
    * C:       dimensional constant [default: 1.]
    * wa_i:    weights of the i-th contribution (e.g. because of sums over k or q points) [default: 1./N_momenta]
    * |D_i|^2: matrix elements squared [default: 1.], must be an array of shape (N_momenta, N_states) or (N_states)
    * b:       Lorentzian broadening of the delta function [default: 0.01 eV, note that b is the HWHM, i.e. 0.5*FWHM]
    * emin, emax, estep:  energy range parameters in eV [default: (0,10, 0.01)]
    
    Outputs:
    * w:       energy range in eV
    * D(w):    DOS or spectral function

    Possible applications:
    * DOS:                 D(w)    =   sum\_vk  w_\k              delta_b(w - e\_vk)
    * JDOS:                D(w)    =   sum\_cvk w_\k              delta_b(w - (e_\ck-e\_vk))
    * IP abs.:             D(w)    = C sum\_cvk w_\k |D_\cvk|^2   delta_b(w - (e_\ck-e\_vk))
    * BSE abs.:            D(w)    = C sum\_l        |R_\l|^2     delta_b(w - E_\l)
    * ph. DOS:             D(w)    =   sum\_mq  w_\q              delta_b(w - w\_mq)
    * ph. Eliashberg fct.: D_ik(w) =   sum\_mq  w_\q T^\ik_\mq(T) delta_b(w - w\_mq)  [ T^\ik_\mq(T)=(\sum_j|g^\ijk^\mq|^2)(2N_\mq(T)+1), only FM term]
    """
    # Energy range
    w = np.arange(emin,emax,estep)

    # Check if we are integrating over momenta or not
    try:
        nmomenta = energies.shape[0]
        nstates  = energies.shape[1]
    except IndexError:
        nmomenta = None
        nstates = energies.shape[0]

    # This is always the case for single-particle quantities
    if nmomenta is not None: 

        if weights is None:   weights   = np.ones(nmomenta)/nmomenta
        if residuals is None: residuals = np.ones((nmomenta,nstates))

        # This is a (nmomenta, nstates, len(w)) array
        energy_poles = w[None, None, :] - energies[:, :, None]

        # Compute the DOS / spectral function
        D = constant*np.einsum('k,kn,knw->w', weights, residuals, broadening / (energy_poles**2 + broadening**2), optimize=True)  

    # Case of the q=0 BSE absorption spectrum 
    else:

        if weights is None:   weights   = np.ones(nstates)
        if residuals is None: residuals = np.ones(nstates)

        # This is a (nstates, len(w)) array
        energy_poles = w[None, :] - energies[:, None]

        # Compute the DOS / spectral function
        D = constant*np.einsum('n,n,nw->w', weights, residuals, broadening / (energy_poles**2 + broadening**2),optimize=True)

    return w,D
