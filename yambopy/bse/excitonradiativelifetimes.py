#
# License-Identifier: GPL
#
# Copyright (C) 2024 The Yambo Team
#
# Authors: RR, FP
#
# This file is part of the yambopy project
#
from yambopy.units import *

def ExcRadLifetimes(yexcdb,statelist=None,degen_step=0.001,gauge='length',verbosity=0,no_cutoff=False):
    """
    This function calculates the excitonic radiative lifetimes according to Eq. (2)
    of Chen et al. (2018): https://pubs.acs.org/doi/10.1021/acs.nanolett.8b01114

    Limitations of the current implementation:
    1. Calculation at Q->0 
    2. Calculation for T=0
    3. Direction of emitted photon parallel to exciton dipole direction
    4. [FP: Only 2D systems with layer plane being xy?]

    :: Input
    - yexcdb : YamboExcitonDB object
    - statelist : list of states for which to compute the lifetime
    - degen_step : merge states as degenerate below this threshold
    - gauge : gauge of the exciton dipole
    - verbosity : print information on merged states
    - nocutoff : if True, set q0_norm=1e-5 even if Coulomb cutoff was used in the BSE run

    TODO: 
        - [FP] Careful check about units + gauge + point 4.
        - Enable T/=0 (solve 2.) 
        - recompute internally exciton residual as e_dir.A.dip (solve 3.)
    """

    if (statelist is None): statelist = np.arange(1,yexcdb.nexcitons+1)
    Omega = yexcdb.lattice.lat_vol
    a1 =    yexcdb.lattice.lat[0]
    a2 =    yexcdb.lattice.lat[1]
    Area = np.linalg.norm(np.cross(a1,a2))
    muS2 = 0
    excE = sorted( np.ma.asarray(yexcdb.eigenvalues.real) )
    excI = np.ma.asarray(yexcdb.l_residual * yexcdb.r_residual)

    tau0_tot = np.zeros(len(statelist))
    merged_states = np.empty(len(statelist), dtype='object')
    
    # q0 norm factor 
    # [BEWARE: if you are using the Coulomb cutoff, you need to have read ndb.cutoff
    #          in yexcdb ]
    if yexcdb.q_cutoff is None or no_cutoff: q0_norm = 1e-5
    else:                                    q0_norm = yexcdb.q_cutoff

    for l, state in enumerate(statelist):
        state_internal = state

        # find states within degen_step window from state
        mask = np.logical_and(excE>=excE[state_internal]-degen_step,excE<=excE[state_internal]+degen_step)
        states_idx = np.where(mask==True)[0]

        #compute gamma i.e. the radiative decay rate
        gamma0 = 0
        for i, st in enumerate(states_idx):
            ES = excE[st]/ha2ev
            # [FP]: is this gauge treatment correct?
            if (gauge=='length'):
                # [RR] if you inspect the Yambo code you might expect 
                # another 1/((2*np.pi)**3) but I think that 
                # d3k_factor/((2np.pi)**3) is actually 1/Omega
                muS2 = excI[st]/(q0_norm**2) 
            elif (gauge == 'velocity'):
                muS2 = excI[st]/(ES**2)
            gg = 4.*np.pi*ES*(muS2/yexcdb.nkpoints)/(Area*speed_of_light)
            gamma0 += gg
        #compute tau, i.e. the radiative lifetime in seconds
        tau0_tot[l] = autime2s/(gamma0.real)
        merged_states[l] = '{}<->{}'.format(min(states_idx)+1,max(states_idx)+1)

    if (verbosity):
        print('=== Exciton radiative lifetimes ===')
        print('#')
        print(f"Global Gauge = {gauge}")
        print(f'Energy degeneration step = {degen_step} eV')
        print(f'q0 norm = {q0_norm}')
        print('#')
        print(f"# {'State':>5} {'Energy (eV)':>15} {'$tau_0$ (s)':>20} {'Merged states':>20}")
        for l,state in enumerate(statelist):
            print(f"{state:>5} {excE[state]:>15} {tau0_tot[l]:>20} {merged_states[l]:>20}")
    return tau0_tot
