#!/usr/bin/python3
#
# Copyright (c) 2017-2018, E. Cannuccia and C. Attaccalite
# All rights reserved.
#
# Calculate thermal lines on single phonon modes or 
# on linear combination of the different phonons
#
# 27/02
# Added a threshold in phonon frequencies to discard the negative and zero frequencies
# (replaces the excluded_freq argument in the previous version)

import numpy as np

from yambopy.zeros import default_freq_thr

def generate_ZG_conf(qe_input, qe_dyn, T=0.0, folder="ZG", freq_thr = default_freq_thr, modes=None, debug=None,skip_ortho_check=True, minus_sign=False):

    atoms      = qe_input.get_atoms("bohr")
    new_atoms  = np.empty((qe_dyn.natoms,3),dtype=float)
    masses     = qe_input.get_masses()

    # Check ortogonaly of the phonon eigenvectors
    # 
    if not skip_ortho_check:
        if not qe_dyn.check_orthogonality():
            print("Error phonon eigenvectors not orthogonal!!! ")
            exit(0)
    else:  
            print("Orthogonality is not checked!!! ")
    
    #
    # Folder name with temperature
    #
    folder=folder+str(T)+"K"
    #
    if modes == None:
        modes=range(0, qe_dyn.nmodes)

    new_filename = qe_input.filename  # default file name

    qe_new=qe_input.copy()

    masses=qe_input.get_masses()
#    qe_dyn.normalize_with_masses(masses,skip_ortho_check)

