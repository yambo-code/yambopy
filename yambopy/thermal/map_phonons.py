#!/usr/bin/python3
#
from qepy.supercell import Supercell
import numpy as np

def Map_Phonons(qe_input, qe_dyn, R, sc_fname=None, dyn_fname=None, debug=None):           # R is the supercell
    
    print(" \n\n\n * * * Map phonons in a supercell * * *\n")
    print(" This code works only without symmetries!!! \n")

    #Check and map phonons
    if qe_dyn.nqpoints == np.prod(R):
        print("Error: number of q-points not compatible with supercell ") 
    


    sc=Supercell(qe_input)
    sc.d_sup(R)

    #write supercell to file
    sc.qe_d.write(sc.qe_d.filename+'_sc')

