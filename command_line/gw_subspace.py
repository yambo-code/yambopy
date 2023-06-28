import os
import numpy as np
from netCDF4 import *
from glob import glob
import argparse

"""
Script to calculate off-diago corrections of yambo ndb.QP databases in order to plot band structure.
For reference, check Appendix A of PHYSICAL REVIEW RESEARCH 2, 043105 (2020)

Inputs:
 1. --fld_diag='path/to/folder/with/diago/dbs' [e.g., the diag/]
 2. --fld_offdiag='path/to/folder/with/offdiagoold/dbs' [e.g., the offdiag/ndb.QP]

This script will prompt the user to go through with updating the dbs.
"""

def get_db(fldr):
    """
    Read NetCDF4 database
    """
    # FP: Not that this function could be replaced by YamboQPDB calls as:
    # $> from yambopy import *
    # $> ds = YamboQPDB.from_db(folder='fldr')
    try: ds=Dataset(f'{fldr}/ndb.QP','r+', format = 'NETCDF4')
    except: raise FileNotFoundError(f'{fldr} Database not found')
    table = np.asarray(ds['QP_table']).T  # [n m k]
    kpts = np.asarray(ds['QP_kpts']).T    # [kx, ky, kz]
    E = np.asarray(ds['QP_E']) # [Re(<n|Sigma|m>) Im(<n|Sigma|m>)]
    Eo = np.asarray(ds['QP_Eo']) # [E_DFT]
    Z = np.asarray(ds['QP_Z'])
    return ds,table, kpts, E, Eo,Z

def create_newdb(flddiag,fldoffdiag):
    """
    Print info for the user and ask if they want to go through with the change
    """
    # Size of the subspace to diagonalize
    dsdiag, table_large, kpts_large, E_large, Eo_large, Z_large = get_db(flddiag)
    print("diago database found:")
    print("=====================")

    dsoffdiag, table_small, kpts_small, E_small, Eo_small, Z_small = get_db(fldoffdiag)
    print("offdiago database found:")
    print("=====================")
    
    N = int(max(table_small[:,0]) - min(table_small[:,0]) + 1)
    # First state
    E_small_diag = np.zeros_like(E_small)
    E_large_new = E_large;
    for ik in range(1,len(kpts_small)+1): 
        cond = np.where(np.equal(table_small[:,2],ik))
        E_mat = np.reshape(E_small[:,0][cond],(N,N),order='F') + 1j*np.reshape(E_small[:,1][cond],(N,N),order='F')
        Eo_mat = np.reshape(Eo_small[cond],(N,N),order='F')
        #   Remove Z
        Z_mat = np.reshape(Z_small[:,0][cond],(N,N),order='F')
        E_mat = Eo_mat + (E_mat - Eo_mat)/Z_mat    
        E_eig_diag,dumb = np.linalg.eig(E_mat)
        E_eig = np.diag(E_eig_diag)
        Emat_order = np.argsort(np.real(np.diag(E_mat)))
        Eeig_order = np.argsort(np.real(np.diag(E_eig)))
        E_mat_diag_Re=np.real(E_eig_diag[Eeig_order])
        E_mat_diag_Im=np.imag(E_eig_diag[Eeig_order])
        #   Recompose the diagonal matrix, preserving the order of the original
        #   eigenstates, we assume that the corrections do not induce a
        #   crossover.
        dumb = np.copy(Emat_order)
        Emat_order2 = np.argsort(dumb)
        E_eig = np.diag(E_mat_diag_Re[Emat_order2])+1j*np.diag(E_mat_diag_Im[Emat_order2])
        #   % Put the Z correction back
        E_eig = Eo_mat + Z_mat*(E_eig - Eo_mat)
        #print(np.concatenate(np.reshape(np.real(E_eig),N*N)))
        #print(np.reshape(np.imag(E_eig),(N*N,1)).flatten())
        E_small_diag[:,0][cond] = np.reshape(np.real(E_eig),(N*N,1),order='F').flatten()
        E_small_diag[:,1][cond] = np.reshape(np.imag(E_eig),(N*N,1),order='F').flatten()
        #E_small_diag=E_small_diag.T
        for n in range(int(min(table_small[:,0])),int(max(table_small[:,0])+1)):
            cond_large2 = np.where(np.equal(table_large[:,2],ik),1,0)
            cond_large0 = np.where(np.equal(table_large[:,0],n),1,0)
            cond_small0 = np.where(np.equal(table_small[:,0],n),1,0)
            cond_small1 = np.where(np.equal(table_small[:,1],n),1,0)
            cond_small2 = np.where(np.equal(table_small[:,2],ik),1,0)
            index_large_new = np.where(cond_large2&cond_large0)[0].tolist()
            index_small_diag = np.where(cond_small2 & cond_small1 & cond_small0)[0].tolist()
            E_large_new[index_large_new] = E_small_diag[index_small_diag]

    filename_out = 'ndb.QP'
    os.system(f'cp {flddiag}/ndb.QP {filename_out}')
    dsout=Dataset(filename_out,'a', format = 'NETCDF4')
    print (f'Writing in {filename_out} \n')
    print (E_large_new)
    var = dsout['QP_E']                     # [Re(<n|Sigma|m>) Im(<n|Sigma|m>)]
    var =E_large_new 
    dsout.close()
    dsoffdiag.close()
    dsdiag.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a new diago dbs after a GW subspace calculation numbers in yambo databases')
    parser.add_argument('-d','--fld_diag', type=str, default="./diago", help='<Required> Path to folder with the ndb.QP diago-database   (Default is ./diago)',required=True)
    parser.add_argument('-o','--fld_offdiag', type=str,default ="./offdiago",help='<Required> Path to folder with the ndb.QP offdiagodiago-database   (Default is ./diago)',required=True)
    args = parser.parse_args()

    fld_diag = args.fld_diag
    fld_offdiag = args.fld_offiag

    create_newdb(fld_diag,fld_offdiag)
