# License-Identifier: GPL
#
# Copyright (C) 2024 The Yambo Team
#
# Authors: FP
#
# This file is part of the yambopy project
#
import numpy as np
from yambopy import YamboLatticeDB,YamboExcitonDB,YamboDipolesDB
from yambopy.units import ha2ev

def exciton_dipoles(blongdir,lattice_path,dipoles_path=None,bse_path=None,kplot=False,check=False):
    """
    This function computes the dipoles D in the excitonic basis (i.e., the residuals)
    at Q=0 starting from the transition-space expression:

    D_{a} = \sum_{cvk} A^{a}_{cvk} Efield \cdot r_cvk 

    :: A      --> BSE eigenvectors
    :: Efield --> electric field
    :: r      --> dipole matrix elements

    This is useful to change Efield direction without rerunning Yambo

    NB: in Yambo, q=0 is actually q=q0_def_norm=1e-5, so if you are checking for 
        consistency with the residuals in BS_diago db you have to multiply 
        by q0_def_norm

    Input:
    * blongdir:   electric field polarization direction
    * exc_states: list of (degenerate) states, e.g. [0,1], [2,3], or [4]
                  degenerate states are summed over squares
    * lattice_path, dipoles_path, bse_path : paths to required databases
    * kplot [default=False]: get k-resolved values (e.g. for plotting)
    * check [default=False]: compare with BSE residuals

    Output:
    * Values of |D|^2 in bohr^{-2}
    * [if kplot=True] Value of D(k) in bohr^{-1} (cmplx)
    """
    if bse_path is None:     bse_path = lattice_path
    if dipoles_path is None: dipoles_path = lattice_path
    if not isinstance(blongdir,np.ndarray): blongdir = np.array(blongdir)

    # Load k-space info
    ylat = YamboLatticeDB.from_db_file(filename=lattice_path+'/ns.db1')
    # Load full BSE database
    yexc = YamboExcitonDB.from_db_file(ylat,filename=bse_path+'/ndb.BS_diago_Q1')
    # Turn off default [1,1,1] dipole projection upon expansion
    ydip = YamboDipolesDB(ylat,save=dipoles_path,filename='ndb.dipoles',project=False)

    # Dipoles are dimensioned as (k,c,v) not (k,v,c) so we switch the table
    table_kcv = yexc.table
    table_kcv[:,[1,2]] = yexc.table[:,[2,1]]

    # Field direction
    field_dir = blongdir/np.linalg.norm(blongdir)

    # Dipole projection along field direction
    dipoles = np.einsum('x,kxcv->kcv',field_dir,ydip.dipoles)

    # Rotate to exciton basis (no-loop fast sum)
    dip_exc = np.sum(yexc.eigenvectors*dipoles[tuple(table_kcv[:,:3].T-1)],axis=1)
    dip_exc = np.abs(dip_exc)**2.
    
    # If k resolution required, we have to do it manually
    if kplot:
        dip_exc_k = np.zeros( (yexc.nexcitons,ylat.nkpoints), dtype=complex)
        for it in range(yexc.nexcitons):
            ik,ic,iv = table_kcv[it,:3]-1
            dip_exc_k[:,ik]+=yexc.eigenvectors[:,it]*dipoles[ik,ic,iv]

    if check==True:
        q0_def_norm=1e-5 # Optical-limit field in Yambo residuals
        dip_exc_reference = np.abs(yexc.l_residual*yexc.r_residual)
        err = np.abs(dip_exc-dip_exc_reference/q0_def_norm**2.)
        av_err  = np.mean(err)
        max_err = np.max(err)
        print(f"Selected field_dir: {blongdir}")
        print(f"Total Average | Max errors: {av_err} | {max_err}")

    if not kplot: return dip_exc
    else:         return dip_exc, dip_exc_k
