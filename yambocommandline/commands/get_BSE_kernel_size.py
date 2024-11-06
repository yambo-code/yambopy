import numpy as np
import argparse

"""
Script to calculate the expected BSE kernel size, in GB, to be loaded by each MPI task

Size_kernel_GB = (nk*nv*nc*ncpl)**2*np*4*2)/(1024.)**3.

Usage:

    :: -nk,--nkpts_bz      = Number of kpoints in the full BZ
    :: -nv,--nvalence,     = Number of valence bands included in the kernel
    :: -nc,--nconduction   = Number of conduction bands included in the kernel
    :: -ndp,--ndoubleprec  = Double precision [Default: False]
    :: -ncpl,--ncoupling   = BSE kernel in coupling mode (i.e., no TDA) [Default: False ]

"""

def get_BSE_kernel_size(nk,nv,nc,np,ncpl):

    bytes_to_GB = (1024.)**3.
    float_bytes = 4.
    ncmplx = 2.

    BSE_kernel_dimension = nk*nv*nc
    BSE_kernel_size = (BSE_kernel_dimension*ncpl)**2.*np*float_bytes*ncmplx/bytes_to_GB

    print('============= BSE kernel size ============')
    print(f' Kernel dimension (resonant part): {BSE_kernel_dimension}')
    print( ' Expected full kernel size [GB]  : %.3f'%BSE_kernel_size)
    print('==========================================')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate expected size of BSE kernel to be loaded by each MPI task (in GB)')
    parser.add_argument('-nk',   '--nkpts_bz', type=int, help='Number of kpoints in the full BZ',required=True)
    parser.add_argument('-nv',   '--nvalence', type=int, help='Number of valence bands included in the kernel',required=True)
    parser.add_argument('-nc','--nconduction', type=int, help='Number of conduction bands included in the kernel',required=True)
    parser.add_argument('-ndp', '--ndoubleprec', type=bool, help='double precision',action='store_false')
    parser.add_argument('-ncpl','--ncoupling', type=bool, help='BSE kernel in coupling mode (i.e., no TDA)',action='store_false')

    args = parser.parse_args()

    nk   = args.nkpts_bz
    nv   = args.nvalence
    nc   = args.nconduction
    np   = 1 if args.ndoubleprec else 2
    ncpl = 1 if args.ncoupling else 2

    get_BSE_kernel_size(nk,nv,nc,np,ncpl)
