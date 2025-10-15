# License-Identifier: GPL
#
# Copyright (C) 2024 The Yambo Team
#
# Authors: FP RR
#
# This file is part of the yambopy project
#
import numpy as np
from yambopy import YamboLatticeDB,YamboExcitonDB,YamboDipolesDB
from yambopy.units import ha2ev
import os

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
    * lattice_path, dipoles_path, bse_path : paths to required databases
    * kplot [default=False]: get k-resolved values (e.g. for plotting)
    * check [default=False]: compare with BSE residuals

    Output:
    * Values of |D|^2 in bohr^2
    * [if kplot=True] Value of D(k) in bohr (cmplx)
    """
    if bse_path is None:     bse_path = lattice_path
    if dipoles_path is None: dipoles_path = lattice_path
    if not isinstance(blongdir,np.ndarray): blongdir = np.array(blongdir)

    # Load k-space info
    ylat = YamboLatticeDB.from_db_file(filename=lattice_path+'/ns.db1')
    # Load full BSE database
    yexc = YamboExcitonDB.from_db_file(ylat,filename=bse_path+'/ndb.BS_diago_Q1')
    # Turn off default [1,1,1] dipole projection upon expansion
    ydip = YamboDipolesDB.from_db_file(ylat,filename=f'{dipoles_path}/ndb.dipoles',project=False)

    # Dipoles are dimensioned as (k,c,v) not (k,v,c) so we switch the table
    table_kcv = yexc.table
    table_kcv[:,[1,2]] = yexc.table[:,[2,1]]

    # Field direction
    field_dir = blongdir/np.linalg.norm(blongdir)

    # Dipole projection along field direction
    dipoles = np.einsum('x,kxcv->kcv',field_dir,ydip.dipoles)

    # Rotate to exciton basis (no-loop fast sum)
    dip_exc = np.sum(yexc.eigenvectors*dipoles[tuple(table_kcv[:,:3].T-1)],axis=1)
    dip_exc_squared = np.abs(dip_exc)**2.
    
    # If k resolution required, we have to do it manually
    if kplot:
        dip_exc_k = np.zeros( (yexc.nexcitons,ylat.nkpoints), dtype=complex)
        for it in range(yexc.nexcitons):
            ik,ic,iv = table_kcv[it,:3]-1
            dip_exc_k[:,ik]+=yexc.eigenvectors[:,it]*dipoles[ik,ic,iv]

    if check==True:
        q0_def_norm=1e-5 # Optical-limit field in Yambo residuals
        dip_exc_reference = np.abs(yexc.l_residual*yexc.r_residual)
        err = np.abs(dip_exc_squared-dip_exc_reference/q0_def_norm**2.)
        av_err  = np.mean(err)
        max_err = np.max(err)
        print(f"Selected field_dir: {blongdir}")
        print(f"Total Average | Max errors: {av_err} | {max_err}")

    if not kplot: return dip_exc_squared
    else:         return dip_exc_squared, dip_exc_k

def exc_dipoles_pol(lattice_path,dipoles_path=None,bse_path=None,save_files=True,dip_file="exc_dipoles.npy",overwrite=False):
    """
    This function computes the dipoles D in the excitonic basis like the one
    above,  but the output is different.

    Output:
    * Complex D_{a,r=x,y,z} for photon emission not projected along Efield 
    direction -- this is used, e.g., for polarization averages in PL spectra.
 
    Input:
    * lattice_path, dipoles_path, bse_path : paths to required databases
    * save_files : whether to save files in `dip_file` .npy database
    * overwrite : if False and `dip_file` is found, load from file
    """

    # Check if we just need to load
    if os.path.exists(dip_file) and overwrite==False:
        print(f'Loading EXCDIP matrices from {dip_file}...')
        dip_exc_loaded = np.load(dip_file)
        return dip_exc_loaded

    if bse_path is None:     bse_path = lattice_path
    if dipoles_path is None: dipoles_path = lattice_path

    # Load k-space info
    ylat = YamboLatticeDB.from_db_file(filename=lattice_path+'/ns.db1')
    # Load full BSE database at Q=0
    yexc = YamboExcitonDB.from_db_file(ylat,filename=bse_path+'/ndb.BS_diago_Q1')
    # Read dipoles in bands range | don't project | don't expand
    # bands range is fixed by BSE calculation
    ydip = YamboDipolesDB.from_db_file(ylat,filename=f'{dipoles_path}/ndb.dipoles',\
                                         bands_range=yexc.bs_bands,project=False,expand=False)
 
    # Expand dipoles
    rot_mats = ylat.sym_car[ylat.kmap[:,1], ...]
    dip_expanded = np.einsum('kij,kjcv->kicv',rot_mats,ydip.dipoles[ylat.kmap[:,0],...],
                             optimize=True)
    time_rev_s = (ylat.kmap[:, 1] >= ylat.sym_car.shape[0]/(int(ylat.time_rev)+1))
    dip_expanded[time_rev_s] = dip_expanded[time_rev_s].conj()

    # Rotate dipoles in exc. basis [n,nblks,nspin,k,c,v] -> [n,k,c,v]
    BS_wfc = np.squeeze( yexc.get_Akcv() ) # Works in TDA and no spin pol
    dip_exc = np.einsum('nkcv,kicv->in',BS_wfc,dip_expanded,
                        optimize=True).astype(dtype=ydip.dipoles.dtype)

    if save_files:
        if dip_file[-4:]!='.npy': dip_file = dip_file+'.npy'
        print(f'Exciton dipoles file saved to {dip_file}')
        np.save(dip_file,dip_exc)

    return dip_exc
