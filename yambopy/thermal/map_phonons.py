#!/usr/bin/python3
#
from qepy.supercell import Supercell

def Map_Phonons(self,qe_input, qe_dyn, R, supercell_name=None, dynmat_name=None):           # R is the supercell
    
    print(" \n\n\n * * * Map phonons in a supercell * * *\n")
    print(" This code works only without symmetries!!! \n")

    sc=Supercell(qe_input)
    sc.d_sup(R)

    #write supercell to file
    sc.qe_d.write(qe_d.filename+'_sc')

    #Check and map phonons
