import warnings
from numba import njit, prange
import os
from netCDF4 import Dataset
from yambopy.letzelphc_interface.lelphcdb import LetzElphElectronPhononDB
from yambopy.dbs.latticedb import YamboLatticeDB
from yambopy.dbs.wfdb import YamboWFDB
from yambopy.dbs.excitondb import YamboExcitonDB
from yambopy.dbs.dipolesdb import YamboDipolesDB
from yambopy.bse.exciton_matrix_elements import exciton_X_matelem
from yambopy.units import *
from yambopy.bse.rotate_excitonwf import rotate_exc_wf
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

class ExcitonPhonon(object):
    def __init__(self, path=None, save='SAVE', lelph_db=None, latdb=None, wfdb=None, \
                 ydipdb=None, bands_range=[], BSE_dir='bse', LELPH_dir='lelph', \
                 DIP_dir='gw',save_files=True):
        if path is None:
            path = os.getcwd()        
        self.path = path
        self.SAVE_dir  = os.path.join(path, save)
        self.BSE_dir   = os.path.join(path,BSE_dir)
        self.LELPH_dir = os.path.join(path,LELPH_dir)
        self.DIP_dir   = os.path.join(path,DIP_dir) # usually dip_dir is in gw run
        self.latdb = latdb
        self.lelph_db = lelph_db
        self.wfdb = wfdb
        self.ydipdb = ydipdb
        
        self.read(lelph_db=lelph_db, latdb=latdb, wfdb=wfdb, \
                  ydipdb=ydipdb, bands_range=bands_range)
        
        self.save_files =save_files # whether the user wants to save files in .npy database

    def read(self, lelph_db=None, latdb=None, wfdb=None,\
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
        self.Dmats = self.wfdb.Dmat()[:,:,0,:,:]
        #self.nbands = max(bands_range) - self.min_bnd
        self.bands_range = bands_range
        

        # set kmap
        kmap = np.zeros((self.wfdb.nkBZ,2), dtype=int)
        kmap[:,0]=self.ydb.kpoints_indexes
        kmap[:,1]=self.ydb.symmetry_indexes
        self.kmap=kmap
        # read exciton database
        self.bs_bands, self.BS_eigs, self.BS_wfcs, self.excQpt = self.read_excdb(self.BSE_dir)

        #read LetzElPhC
        try:
            ndb_lelph_fname = os.path.join(self.LELPH_dir, 'ndb.elph')
            if lelph_db :
                self.lelph_db = lelph_db
            else :
                self.lelph_db = LetzElphElectronPhononDB(filename = ndb_lelph_fname)        
        except Exception as e:
            raise IOError(f'Cannot read ndb.elph file: {e}')        
        
        self.kpts = self.lelph_db.kpoints
        self.qpts = self.lelph_db.qpoints
        self.elph_bnds_range = self.lelph_db.bands
        self.ph_freq = self.lelph_db.ph_energies/ha2ev # yambopy gives energies in Ry, I work in Hartree
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
    
    def compute_Exph(self, gamma_only = True):
        from time import time
        """Compute exciton-phonon coupling matrix elements"""
        time_exph = 0
        time_elph_io = 0
        time_ex_rot = 0
        time_exph_io = 0
        ### compute ex-ph matrix elements:
        ex_ph = np.zeros((len(self.excQpt), self.qidx_in_kpts.shape[0], self.lelph_db.nm, \
                          self.BS_eigs.shape[1],self.BS_eigs.shape[1]), dtype = self.BS_eigs.dtype)
        kpts_ibz = self.wfdb.kpts_iBZ
        ### build a kdtree for kpoints
        print('Building kD-tree for kpoints')
        kpt_tree = build_ktree(self.kpts)
        self.kpt_tree = kpt_tree
        ### fomd tje omdoces pf qpoints in kpts
        self.qidx_in_kpts = find_kpt(self.kpt_tree, self.kpts)
        
        self.qplusQ_in_kpts = np.zeros((len(self.excQpt), self.kpts.shape[0]), dtype = int)
        for iQ in range(0,len(self.excQpt)):
            self.qplusQ_in_kpts[iQ] = find_kpt(self.kpt_tree, self.kpts + self.excQpt[iQ])
        
        print('Computing Exciton-phonon matrix elements for phonon absorption ...')
        
        if gamma_only:
            print("Only at the Gamma point")
            iQ = 0  # Only run for Gamma
            for i in tqdm(range(self.qidx_in_kpts.shape[0]), desc="Exciton-phonon Gamma"):
                # < S|dv|0>
                # Basic check
                kq_diff = self.qpts[i] - self.kpts[self.qidx_in_kpts[i]]
                kq_diff = kq_diff - np.rint(kq_diff)
                assert (np.linalg.norm(kq_diff) < 1e-5)

                tik = time()
                # Read elph_matrix elements
                _, eph_mat_iq = self.lelph_db.read_iq(i, bands_range=self.bands_range, convention='standard')
                eph_mat_iq = eph_mat_iq[:, :, 0, :, :].transpose(1, 0, 3, 2)
                time_elph_io += time() - tik

                # Get rotated exciton wavefunctions
                tik = time()
                ik_ibz, isym = self.kmap[self.qidx_in_kpts[i]]
                ik_ibz_qplusQ, isym_qplusQ = self.kmap[self.qplusQ_in_kpts[iQ,i]]

                is_sym_time_rev = isym >= self.symm_mats.shape[0] // (int(self.ele_time_rev) + 1)

                Ak = rotate_exc_wf(
                    self.BS_wfcs[ik_ibz], self.sym_red[isym], self.kpts.data,
                    kpts_ibz[ik_ibz], self.Dmats[isym], is_sym_time_rev, ktree=self.kpt_tree
                )
                Akq = rotate_exc_wf(
                    self.BS_wfcs[ik_ibz_qplusQ], self.sym_red[isym_qplusQ], self.kpts.data,
                    kpts_ibz[ik_ibz_qplusQ], self.Dmats[isym_qplusQ], is_sym_time_rev, ktree=self.kpt_tree
                )
                time_ex_rot += time() - tik

                # Compute exciton-phonon matrix element
                tik = time()
                ex_ph_tmp = exciton_X_matelem(
                    self.excQpt[ik_ibz_qplusQ],
                    self.kpts[self.qidx_in_kpts[i]],
                    Akq,
                    Ak,
                    eph_mat_iq,
                    self.kpts,
                    contribution='b'
                )
                ex_ph[iQ, i, ...] = ex_ph_tmp
            time_exph = time_exph + time() - tik
        else:
            for iQ in tqdm(range(len(self.excQpt)), desc="Exciton-phonon all Q"):
                for i in tqdm(range(self.qidx_in_kpts.shape[0]), desc=f"Exciton-phonon Q[{iQ}]"):
                        # < S|dv|0>
                        ## first do a basic check
                        kq_diff = self.qpts[i] - self.kpts[self.qidx_in_kpts[i]]
                        kq_diff = kq_diff - np.rint(kq_diff)
                        assert (np.linalg.norm(kq_diff) < 10**-5)
                        tik = time()
                        ## read elph_matrix elements
                        _,eph_mat_iq = self.lelph_db.read_iq(i, bands_range = self.bands_range,convention='standard') # this has to be standard for exciton_X_matelem
                        eph_mat_iq=eph_mat_iq[:,:,0,:,:].transpose(1,0,3,2) # take spin 0 and I don't know why MN swap initial and final bands
                        time_elph_io = time_elph_io + time() - tik
                        ## get rotated ex-wfc
                        tik = time()
                        ik_ibz, isym = self.kmap[self.qidx_in_kpts[i]]  ## get the ibZ kpt and symmetry matrix for this q/k point
                        ik_ibz_qplusQ, isym_qplusQ = self.kmap[self.qplusQ_in_kpts[iQ,i]]

                        is_sym_time_rev = False
                        if (isym >= self.symm_mats.shape[0] / (int(self.ele_time_rev) + 1)):
                            is_sym_time_rev = True
                        Ak = rotate_exc_wf(self.BS_wfcs[ik_ibz], self.sym_red[isym], self.kpts.data, \
                            kpts_ibz[ik_ibz], self.Dmats[isym], is_sym_time_rev, ktree=self.kpt_tree
                            )#kpts_ibz is the momentum of exciton #
                        Akq = rotate_exc_wf(self.BS_wfcs[ik_ibz_qplusQ], self.sym_red[isym_qplusQ], self.kpts.data,\
                                            kpts_ibz[ik_ibz_qplusQ], self.Dmats[isym_qplusQ], is_sym_time_rev, ktree = self.kpt_tree)
                        ## compute ex-ph matrix
                        time_ex_rot = time_ex_rot + time() - tik
                        tik = time()
                        ex_ph_tmp = exciton_X_matelem(
                                        self.excQpt[ik_ibz_qplusQ], # Q+q (for iQ = 0 -> q)
                                        self.kpts[self.qidx_in_kpts[i]], # q
                                        Akq, # <S', Q+q|
                                        Ak, #|S,Q>
                                        eph_mat_iq, # <k+q| dvscf| k>
                                        self.kpts, 
                                        contribution='b'
                        )
                        ex_ph[iQ, i,...] = ex_ph_tmp
            time_exph = time_exph + time() - tik

        ## SAving exciton phonon matrix elements
        if(self.save_files):
            print('Saving exciton phonon matrix elements')
            tik_exph = time()
            ex_ph = np.array(ex_ph).astype(self.BS_wfcs.dtype)
            #     # (iq, modes, initial state (i), final-states (f)) i.e <f|dv_Q|i> for phonon absoption
            np.save('Ex-ph', ex_ph)
            time_exph_io = time_exph_io + time() - tik_exph
        else:
            print('Loading exciton phonon matrix elements')
            tik_exph = time()
            ex_ph = np.load('Ex-ph.npy')
            time_exph_io = time_exph_io + time() - tik_exph
        self.ex_ph = ex_ph

        print('** Timings **')
        print(f'ELPH IO   : {time_elph_io:.4f} s')
        print(f'Ex-ph IO  : {time_exph_io:.4f} s')
        print(f'Ex rot    : {time_ex_rot:.4f} s')
        print(f'Ex-ph     : {time_exph:.4f} s')
        print('*' * 30, ' Program ended ', '*' * 30)