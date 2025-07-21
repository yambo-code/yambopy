import warnings
from numba import njit, prange
import os
from netCDF4 import Dataset
import torch as pytorch
from yambopy import YamboLatticeDB
from yambopy.dbs.wfdb import YamboWFDB
from yambopy.dbs.excitondb import YamboExcitonDB
from yambopy.dbs.dipolesdb import YamboDipolesDB
from yambopy.units import *
try:
    from pykdtree.kdtree import KDTree 
    ## pykdtree is much faster and is recommanded
    ## pip install pykdtree
    ## useful in Dmat computation
except ImportError as e:
    from scipy.spatial import KDTree
from yambopy.kpoints import build_ktree, find_kpt

from tqdm import tqdm
warnings.filterwarnings('ignore')

class ExcitonDipole(object):
    def __init__(self, path=None, save='SAVE', latdb=None, wfdb=None, \
                ydipdb=None, bands_range=[], BSE_dir='bse', \
                 DIP_dir='gw',save_files=True):
        self.SAVE_dir  = os.path.join(path, save)
        self.BSE_dir   = os.path.join(path,BSE_dir)
        self.DIP_dir   = os.path.join(path,DIP_dir) # usually dip_dir is in gw run
        self.latdb = latdb
        self.wfdb = wfdb
        self.ydipdb = ydipdb
        
        self.read(latdb=latdb, wfdb=wfdb, \
                  ydipdb=ydipdb, bands_range=bands_range)
        
        self.save_files =save_files # whether the user wants to save files in .npy database

    def read(self, latdb=None, wfdb=None,\
             ydipdb = None, bands_range = []):
        # Open the ns.db1 database to get essential data
        SAVE_dir = self.SAVE_dir
        # readlatdb        
        try:
            ns_db1_fname = os.path.join(SAVE_dir, 'ns.db1')
            if latdb :
                if not hasattr(latdb,'ibz_kpoints'): latdb.expand_kpoints()
                self.ydb = latdb
            else :
                self.ydb = YamboLatticeDB.from_db_file(ns_db1_fname, Expand=True)        
        except Exception as e:
            raise IOError(f'Cannot read ns.db1 file: {e}')

        self.lat_vecs = self.ydb.lat
        self.nibz = self.ydb.ibz_nkpoints
        self.symm_mats = self.ydb.sym_car
        self.ele_time_rev = self.ydb.time_rev
        self.blat_vecs = self.ydb.rlat.T

        #readwfdb
        try:
            ns_wfdb_fname = os.path.join(SAVE_dir, 'ns.wf')
            if wfdb :
                if not hasattr(latdb,'save_Dmat'): wfdb.Dmat()
                self.wfdb = wfdb
            else :
                self.wfdb = YamboWFDB(filename = ns_wfdb_fname, Expand=True, latdb=self.ydb, bands_range=bands_range)  
        except Exception as e:
            raise IOError(f'Cannot read ns.wf file: {e}')
        #Read dimensions
        self.nkpoints = self.wfdb.nkpoints # Number of k-points in iBZ
        self.nspin    = self.wfdb.nspin       # Number of spin components
        self.nspinor  = self.wfdb.nspinor   # Number of spinor components
        self.nbands   = self.wfdb.nbands     # Number of bands
       
        nbnds = max(bands_range)-min(bands_range)
        start_bnd_idx = 0
        end_bnd = start_bnd_idx + nbnds
        self.bands_range = bands_range

        # set kmap
        kmap = np.zeros((self.wfdb.nkBZ,2), dtype=int)
        kmap[:,0]=self.ydb.kpoints_indexes
        kmap[:,1]=self.ydb.symmetry_indexes
        self.kmap=kmap
        # read exciton database
        self.bs_bands, self.BS_eigs, self.BS_wfcs, self.excQpt = self.read_excdb(self.BSE_dir)

        #read YamboDipolesDb
        try:
            ndb_dipoles_fname = os.path.join(self.DIP_dir, 'ndb.dipoles')
            if ydipdb :
                self.ydipdb = ydipdb
            else :
                self.ydipdb = YamboDipolesDB(self.ydb, save='',filename = ndb_dipoles_fname, dip_type='iR',\
                                           field_dir=[1,1,1],project=False, expand=False,bands_range=bands_range,\
                                            )        
        except Exception as e:
            raise IOError(f'Cannot read ndb.dipoles file: {e}')
        if(self.ydipdb.spin == 2):
            self.ele_dips = self.ydipdb.dipoles.conjugate().transpose(1,2,3,4,0)
        if(self.ydipdb.spin == 1):
            self.ele_dips = self.ydipdb.dipoles.conjugate().transpose(0,2,3,1)

        ### build a kdtree for kpoints
        print('Building kD-tree for kpoints')
        self.kpts = self.ydb.red_kpoints
        kpt_tree = build_ktree(self.kpts)
        self.kpt_tree = kpt_tree
        ### fomd tje omdoces pf qpoints in kpts
        self.qidx_in_kpts = find_kpt(self.kpt_tree, self.kpts)
        # Remember b_lat @ red_kpoints = car_kpoints -> red_kpoins = inv(b_lat) @ car_kpoints
        # (inv_blat) = (ydb.lat.T)
        temp = np.matmul(self.symm_mats, self.blat_vecs)  # shape (n, j, l)
        # temp: (n, j, l)
        # lat_vecs: (i, j)
        # reshape lat_vecs for batched matmul: (1, i, j)
        # use matmul: (1, i, j) @ (n, j, l) → (n, i, l)
        sym_red = np.matmul(self.lat_vecs[None, :, :], temp)  # result (n, i, l)
        self.sym_red = np.rint(sym_red).astype(int)

    def read_excdb(self, BSE_dir):
        """Read yambo exciton database for each Q-point"""
        bs_bands = [] # bands involved in BSE
        BS_eigs  = [] #eigenenergies BSE
        BS_wfcs = [] # exciton wavefunctions
        excQpt  = [] #Q-point of BSE -> The q o A^{\lambda Q}_{cvk}
        for iq in tqdm(range(self.nibz), desc="Loading Ex-wfcs "):
            try:
                bse_db_iq = YamboExcitonDB.from_db_file(self.ydb, folder=BSE_dir,filename=f'ndb.BS_diago_Q{iq+1}')
            except Exception as e:
                raise IOError(f'Cannot read ndb.BS_diago_Q{iq} file: {e}')
            bs_bands=bse_db_iq.nbands
            tmp_eigs = bse_db_iq.eigenvalues
            tmp_wfcs = bse_db_iq.get_Akcv()
            tmp_qpt = self.ydb.lat @ bse_db_iq.car_qpoint
            BS_eigs.append(tmp_eigs)
            BS_wfcs.append(tmp_wfcs)
            excQpt.append(tmp_qpt)
        return bs_bands, (np.array(BS_eigs)/ha2ev).astype(self.wfdb.wf.dtype), np.array(BS_wfcs).astype(self.wfdb.wf.dtype), excQpt
    
    def compute_Exdipole(self):
        from time import time
        """Compute exciton-photon coupling matrix elements"""
        time_exdip = 0
        time_eldip_io = 0
        time_ex_rot = 0
        time_exdip_io = 0
        ### build a kdtree for kpoints
        print('Building kD-tree for kpoints')
        kpt_tree = build_ktree(self.kpts)
        self.kpt_tree = kpt_tree
        ### fomd tje omdoces pf qpoints in kpts
        self.qidx_in_kpts = find_kpt(self.kpt_tree, self.kpts)

        print('Computing Exciton-photon matrix elements')
        self.ex_dip = self.exe_dipoles(self.ele_dips, self.BS_wfcs[0],
                            self.kmap,self.symm_mats,self.ele_time_rev)
        np.save('Ex-dipole', self.ex_dip)
        
        print('** Timings **')
        print(f'ELPH IO   : {time_eldip_io:.4f} s')
        print(f'Ex-ph IO  : {time_exdip_io:.4f} s')
        print(f'Ex rot    : {time_ex_rot:.4f} s')
        print(f'Ex-ph     : {time_exdip:.4f} s')
        print('*' * 30, ' Program ended ', '*' * 30)

    def exe_dipoles(self,ele_dipoles, exe_wfc_gamma, kmap, symm_mats, time_rev):
        ## exe_dipoles are for photon emission
        ## ele_dipoles are in iBZ (k,c,v,pol)
        nv = exe_wfc_gamma.shape[-1]
        rot_mats = symm_mats[kmap[:, 1], ...]
        dip_expanded = np.einsum('kij,kcvj->kcvi', rot_mats, ele_dipoles[kmap[:, 0],
                                                                        ...])
        time_rev_s = (kmap[:, 1] >= symm_mats.shape[0] / (int(time_rev) + 1))
        dip_expanded[time_rev_s] = dip_expanded[time_rev_s].conj()
        return np.einsum('nkcv,kcvi->in',
                        exe_wfc_gamma.conj(),
                        dip_expanded,
                        optimize=True).conj().astype(dtype=self.ydipdb.dipoles.dtype)  ##(pol,nexe)