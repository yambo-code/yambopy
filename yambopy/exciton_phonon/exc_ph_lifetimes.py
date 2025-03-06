#
# License-Identifier: GPL
#
# Copyright (C) 2024 The Yambo Team
#
# Authors: FP
#
# This file is part of the yambopy project
#
"""
This class calculates the exciton-phonon lifetimes using the following expression:

    [[MATH]]

Input:

    :: G2[Nexc_in,Nq,Nexc_sum,Nph] --> array of |G_exc-ph|^2 [UNITS]
    :: Eexc_in[Nexc_in]            --> Energies of Q=0/initial exc. states [UNITS]
    :: Eexc_sum[Nq,Nexc_sum]       --> Energies of scattered/internal exc. states [UNITS]
    :: Eph[Nq,Nmodes]              --> Phonon energies [UNITS]
    :: T                           --> Lattice temperature [UNITS] 

Implementation steps:
    - Load exc. energies IN and SUM
    - Possible stretching factor E_iq = E_iq + alpha*|iq|
    - deg. finder Eexc_sum at q=0, take average of deg. states
    - Load phonon energies
    - Load exc-ph matrix elements
    - Info: minimum exc energy at which q-point, Fan threshold, broadening
    - Lifetime formula evaluation in full BZ: loops q, exc_in, exc_sum, ph; expression is |G|^2 * POLE; additional output q-dependent lifetime and poles
    - deg. finder Eexc_in at q=0, take average of lifetimes for that state [Hidden option]
    - Pole: 
        - if ph_E<FAN_thresh pole is zero
        - if E_in(q=0) degenerate with E_sum(q=0) pole is zero [Hidden option]
        - if E_in is exactly the same state a E_out, pole is zero [Hidden option]
        - pole includes occupation functions [Hidden option to remove excitonic one]
Notes:
    - At q=0, degenerate states are set at exactly equal energy values (average)
    - By default the sum over exc. states does not include self-scattering 
    - By default self-scattering is excluded over full degenerate subspaces
    - By default lifetime values over degenerate states are averaged
"""
