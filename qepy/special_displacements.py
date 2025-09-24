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
import math

from yambopy.zeros import default_freq_thr
from yambopy.tools.funcs import bose
from yambopy.units import au2kelvin,amu2au

def generate_ZG_conf(qe_input, qe_dyn, T=0.0, folder="ZG", freq_thr = default_freq_thr, modes=None, debug=None,skip_ortho_check=True, minus_sign=False):

    print("\n\n * * * Special Displacement Generation * * * \n\n")

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

    qe_new =qe_input.copy()

    masses=qe_input.get_masses()
    # qe_dyn.normalize_with_masses(masses) 
    # Not need anymore
    # already included in qe.displace
    #
    # Fix Gauge sign of eigenvalues   
    # see page 075125-3 of PRB 94, 075125 (2006)
    #
    for im in modes:
        if  qe_dyn.eiv[0,im,0] < 0.0:
            qe_dyn.eiv[0,im,:] *= -1.0
    #
    cart_mode=np.zeros([qe_dyn.natoms,3],float)
    for im in modes:
       #
       # Gaussian width
       #
       # see Eq. 12 of arXiv:1512.06377v1
       #
       w_au = qe_dyn.get_phonon_freq(0,im+1,unit='Ha')
       if w_au > freq_thr:
           q_0  = 1.0/math.sqrt(2.0*w_au)
           q_T  = q_0*math.sqrt(1.0+2.0*bose(w_au,T/au2kelvin))
       else:
           continue
       
       if debug is not None:
           print("Mode: "+str(im))
           print("W and T atomic units : %14.10f, %14.10f " % (w_au,T/au2kelvin))
           print("W in cm-1 %14.10f " % qe_dyn.get_phonon_freq(0,im+1,unit='cm-1'))
           print("Amplitude at T=0     : %14.10f " % q_0)
           print("Amplitude at finite T: %14.10f " % q_T)

       #
       i_sign=1.0
       if minus_sign:
           i_sign=-1.0
       if (im % 2 ) ==0:
           delta =  i_sign*q_T
       else:
           delta = -i_sign*q_T
    
       for a in range(qe_dyn.natoms):
           e = qe_dyn.eiv[0,im,a*3:(a+1)*3]
           cart_mode[a][:]=cart_mode[a][:]+e.real*delta/math.sqrt(amu2au) #/np.sqrt(masses[a]*amu2au)

    qe_new.displace(cart_mode,1.0)

    suffix="_ZG" 
    if minus_sign:
        suffix=suffix+"m"

    qe_new.control['prefix']=qe_input.control['prefix'].strip("'")+suffix

    if not debug:
        qe_new.write(str(new_filename)+suffix)
    else:
        print("ZG line: ")
        print(qe_new.get_atoms())
