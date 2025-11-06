#!/usr/bin/python3
#
from qepy.supercell import Supercell
import numpy as np
import sys

"""
Map phonon in a supercell
no_invar_ph = remove phonon modes invariant under inversion symmetry
                  modulo a reciprocal lattice vector
"""

def Map_Phonons(qe_input, qe_dyn, R, no_invar_ph=None, sc_fname=None, dyn_fname=None, debug=None):           # R is the diagonal supercell
    
    print(" \n\n\n * * * Map phonons in a supercell * * *\n")
    print(" This code works only without symmetries!!! \n")

    if debug:
        print(" Supercell : ",str(R))

    #Check and map phonons
    if qe_dyn.nqpoints != np.prod(R):
        print("Error: number of q-points not compatible with supercell ") 
        print("      ",str(qe_dyn.nqpoints),"  vs ",str(np.prod(R)))
        sys.exit(0)
        # Better check can be implemented

    sc=Supercell(qe_input)
    sc.d_sup(R)

    #write supercell to file
    sc.qe_d.write(sc.qe_d.filename+'_sc')

    #map_phonons
    qe_mapped=qe_dyn.expand_in_supercell(sc.qe_d)
    qe_mapped.write_modes(filename="matdyn_sc.modes")
