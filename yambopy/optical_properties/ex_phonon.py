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
from time import time
warnings.filterwarnings('ignore')

class ExcitonPhonon(object):
    """
    This class contains the methods to compute the exciton-phonon matrix elements.
    
    Parameters
    ----------
    path : str, optional
        The path where the Yambo calculation is performed. Default is the current
        working directory.
    save : str, optional
        The name of the folder which contains the Yambo save folder. Default is 'SAVE'.
    lelph_db : LetzElphElectronPhononDB, optional
        The LetzElphElectronPhononDB object which contains the electron-phonon matrix
        elements. If not provided, it will be read from the lelph database.
    latdb : YamboLatticeDB, optional
        The YamboLatticeDB object which contains the lattice information. If not
        provided, it will be read from the ns.db1 file.
    wfdb : YamboWFDB, optional
        The YamboWFDB object which contains the wavefunction information. If not
        provided, it will be read from the ns.wf file.
    ydipdb : YamboDipolesDB, optional
        The YamboDipolesDB object which contains the dipole information. If not
        provided, it will be read from the dipoles database.
    bands_range : list or tuple, optional
        The range of bands for which the exciton-phonon matrix elements will be
        computed. Default is all bands.
    BSE_dir : str, optional
        The name of the folder which contains the BSE calculation. Default is 'bse'.
    LELPH_dir : str, optional
        The name of the folder which contains the electron-phonon matrix elements.
        Default is 'lelph'.
    DIP_dir : str, optional
        The name of the folder which contains the dipole information. Default is 'gw'.
    save_files : bool, optional
        If True, the matrix elements will be saved in .npy files. Default is True.
    
    Attributes
    ----------
    SAVE_dir : str
        The path of the SAVE folder.
    BSE_dir : str
        The path of the BSE folder.
    LELPH_dir : str
        The path of the folder which contains the electron-phonon matrix elements.
    DIP_dir : str
        The path of the folder which contains the dipole information.
    latdb : YamboLatticeDB
        The YamboLatticeDB object which contains the lattice information.
    lelph_db : LetzElphElectronPhononDB
        The LetzElphElectronPhononDB object which contains the electron-phonon matrix
        elements.
    wfdb : YamboWFDB
        The YamboWFDB object which contains the wavefunction information.
    ydipdb : YamboDipolesDB
        The YamboDipolesDB object which contains the dipole information.
    save_files : bool
        If True, the matrix elements will be saved in .npy files.
    """
    def __init__(self, path=None, save='SAVE', lelph_db=None, latdb=None, wfdb=None, \
                 ydipdb=None, bands_range=[], BSE_dir='bse', LELPH_dir='lelph', \
                 DIP_dir='gw',save_files=True,neig=-1):
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
        self.neig = neig
        
        self.read(lelph_db=lelph_db, latdb=latdb, wfdb=wfdb, \
                  ydipdb=ydipdb, bands_range=bands_range)
        
        self.save_files =save_files # whether the user wants to save files in .npy database

        # Tricks for caching
        if not hasattr(self, '_Ak_cache'):
            self._Ak_cache = {}
        if not hasattr(self, '_eph_mat_cache'):
            self._eph_mat_cache = {}

    def read(self, lelph_db=None, latdb=None, wfdb=None,\
             ydipdb = None, bands_range = []):
        # Open the ns.db1 database to get essential data
        """
        Read in the YamboLatticeDB, YamboWFDB, LetzElphElectronPhononDB, YamboDipolesDB
        objects and set the kmap, bs_bands, BS_eigs, BS_wfcs, excQpt, kpts, qpts, 
        elph_bnds_range, and ph_freq attributes.

        Parameters
        ----------
        lelph_db : LetzElphElectronPhononDB
            The LetzElphElectronPhononDB object which contains the electron-phonon matrix
            elements.
        latdb : YamboLatticeDB
            The YamboLatticeDB object which contains the lattice information.
        wfdb : YamboWFDB
            The YamboWFDB object which contains the wavefunction information.
        ydipdb : YamboDipolesDB
            The YamboDipolesDB object which contains the dipole information.
        bands_range : list
            A list of two integers which define the range of bands to be read.
        """
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
                if not hasattr(wfdb,'save_Dmat'): wfdb.Dmat()
                self.wfdb = wfdb
            else :
                self.wfdb = YamboWFDB(filename = ns_wfdb_fname, latdb=self.ydb, bands_range=bands_range)  
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
        self.bs_bands, self.BS_eigs, self.BS_wfcs, self.excQpt = self.read_excdb(self.BSE_dir, self.neig)

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

    def read_excdb(self, BSE_dir, neig):
        """
        Read yambo exciton database for each Q-point from the specified BSE directory.

        Parameters
        ----------
        BSE_dir : str
            The directory containing the BSE calculation data.

        Returns
        -------
        bs_bands : list
            The list of bands involved in the BSE calculation.
        BS_eigs : numpy.ndarray
            The eigenenergies of the BSE, converted to the appropriate units.
        BS_wfcs : numpy.ndarray
            The exciton wavefunctions.
        excQpt : list
            The Q-points associated with the BSE, calculated as the product of the
            lattice vectors and the Cartesian Q-points.
        
        Raises
        ------
        IOError
            If a required BSE database file cannot be read.
        """
        bs_bands = [] # bands involved in BSE
        BS_eigs  = [] #eigenenergies BSE
        BS_wfcs = [] # exciton wavefunctions
        excQpt  = [] #Q-point of BSE -> The q o A^{\lambda Q}_{cvk}
        for iq in tqdm(range(self.nibz), desc="Loading Ex-wfcs "):
            try:
                bse_db_iq = YamboExcitonDB.from_db_file(self.ydb, folder=BSE_dir,filename=f'ndb.BS_diago_Q{iq+1}',neigs=neig)
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
        """
        Compute exciton-phonon matrix elements.

        Parameters
        ----------
        gamma_only : bool
            If True, only compute the matrix elements for the gamma point.
            Otherwise, compute matrix elements for all q-points.

        Notes
        -----
        This function is the main entry point for the computation of exciton-phonon
        matrix elements. It first performs some precomputations, then performs the
        main computation, and finally saves the results to disk.

        The precomputations involve preparing a k-D tree of the k-points, finding
        the indices of the q-points in the k-point list, and finding the indices
        of q+Q in the k-point list.

        The main computation involves computing the matrix elements of the
        exciton-phonon interaction. This is done by first computing the matrix
        elements for the gamma point (if gamma_only is True), and then computing
        the matrix elements for all q-points. The matrix elements are computed
        using the _compute_gamma_only_exph and _compute_full_exph functions.

        The results are then saved to disk using the _save_or_load_exph function.

        Finally, the timings are printed.

        """
        from time import time

        # Timing key/values
        self.timings = {
        'elph_io': 0,
        'ex_rot': 0,
        'exph': 0,
        'exph_io': 0,            
        }

        # Precomputations
        self._prepare_kdtree()
        self._find_qidx_in_kpts()
        self._find_qplusQ_indices()

        # Main computation
        if gamma_only:
            ex_ph = self._compute_gamma_only_exph()
        else:
            ex_ph = self._compute_full_exph()

        # I/O
        self.ex_ph = self._save_or_load_exph(ex_ph)
        
        # Profiling
        print('** Timings **')
        for key, val in self.timings.items():
            print(f'{key.upper():<10}: {val:.4f} s')
        print('*' * 30, ' Program ended ', '*' * 30)

    def _prepare_kdtree(self):
        """
        Build a k-D tree of the k-points.

        This function builds a k-D tree of the k-points, which is used to
        efficiently find the indices of the q-points and q+Q in the k-point
        list.

        The timings are stored in the 'elph_io' key of the timings dictionary.

        """
        print('Build KD-tree for k-points...')
        self.kpt_tree = build_ktree(self.kpts)
    
    def _find_qidx_in_kpts(self):
        """
        Find the indices of the q-points in the k-point list.

        This function finds the indices of the q-points in the k-point list
        using the k-D tree built in _prepare_kdtree.

        The result is stored in the `qidx_in_kpts` attribute.

        """
        self.qidx_in_kpts = find_kpt(self.kpt_tree,self.kpts)

    def _find_qplusQ_indices(self):
        """
        Find the indices of the q+Q points in the k-point list.

        This function finds the indices of the q+Q points in the k-point list
        using the k-D tree built in _prepare_kdtree.

        The result is stored in the `qplusQ_in_kpts` attribute.

        """
        nq = self.kpts.shape[0]
        nQ = len(self.excQpt)
        self.qplusQ_in_kpts = np.zeros((nQ, nq), dtype=int)
        for iQ in range(nQ):
            self.qplusQ_in_kpts[iQ] = find_kpt(self.kpt_tree, self.kpts + self.excQpt[iQ])

    def _compute_gamma_only_exph(self):
        'ex-ph matrix elements only at the \u0393 point'
        nq = self.qidx_in_kpts.shape[0]
        nb = self.BS_eigs.shape[1]
        nm = self.lelph_db.nm
        ex_ph = np.zeros((1, nq, nm, nb, nb), dtype=self.BS_eigs.dtype)

        print("Computing exciton-phonon matrix elements at \u0393 point...")

        for iq in tqdm(range(nq), desc="Exciton-phonon \u0393"):
            # Validate q -> k mapping
            kq_diff = self.qpts[iq] - self.kpts[self.qidx_in_kpts[iq]]
            kq_diff -= np.rint(kq_diff)
            assert np.linalg.norm(kq_diff) < 1e-5, f"Inconsistent q→k mapping at q index {iq}"

            # Load elph matrix
            t0 = time()
            eph_mat_iq = self._load_elph_matrix(iq)
            self.timings['elph_io'] += time() - t0

            # Rotate exciton wavefunctions
            t0 = time()
            _, Akq = self._rotate_exciton_pair(iq, iQ=0)
            self.timings['ex_rot'] += time() - t0
            Ak0,_ = self._rotate_exciton_pair(0,0)
            # Compute ex-ph element
            t0 = time()
            ex_ph[0,iq] = exciton_X_matelem(
                                    self.excQpt[0],         # Q
                                    self.lelph_db.qpoints[iq],
                                    Akq, Ak0 ,eph_mat_iq, self.wfdb.kBZ,
                                    contribution='b', diagonal_only=False, ktree=self.wfdb.ktree,
                                    )
            self.timings['exph'] += time() - t0

        return ex_ph

    def _compute_full_exph(self):
        """
        Compute exciton-phonon matrix elements for all Q points.

        Parameters
        ----------
        None

        Returns
        -------
        ex_ph : ndarray, shape (nQ, nq, nm, nb, nb)
            Exciton-phonon matrix elements, where nQ is the number of Q points,
            nq is the number of q points, nm is the number of phonon modes, and
            nb is the number of bands.
        """
        nq = self.qidx_in_kpts.shape[0]
        nQ = len(self.excQpt)
        nb = self.BS_eigs.shape[1]
        nm = self.lelph_db.nm
        ex_ph = np.zeros((nQ, nq, nm, nb, nb), dtype=self.BS_eigs.dtype)
        print("Computing exciton-phonon matrix elements for all Q points...")

        for iQ in tqdm(range(nQ), desc="Exciton-phonon Q"):
            for iq in range(nq):
                # Check q consistency
                kq_diff = self.qpts[iq] - self.kpts[self.qidx_in_kpts[iq]]
                kq_diff -= np.rint(kq_diff)
                assert np.linalg.norm(kq_diff) < 1e-5, f"Inconsistent q→k at q={iq}, Q={iQ}"

                # Load elph matrix
                t0 = time()
                eph_mat_iq = self._load_elph_matrix(iq)
                self.timings['elph_io'] += time() - t0

                # Rotate exciton wavefunctions
                t0 = time()
                Ak, Akq = self._rotate_exciton_pair(iq, iQ)
                Ak0, Akq0 = self._rotate_exciton_pair(iq=0, iQ = iQ)

                self.timings['ex_rot'] += time() - t0

                # Compute ex-ph element
                t0 = time()
                #Evaluate the exciton-phonon matrix element ⟨Akq|g(q)|Ak⟩.
                ex_ph[iQ, iq] = exciton_X_matelem(
                                    self.excQpt[iQ],         # Q
                                    self.lelph_db.qpoints[iq],  # q
                                    Akq, Ak0, eph_mat_iq, self.wfdb.kBZ,
                                    contribution='b', diagonal_only=False, ktree=self.wfdb.ktree
                                    )
                self.timings['exph'] += time() - t0

        return ex_ph

    def _load_elph_matrix(self, iq):
        """
        Load the electron-phonon matrix for a given q-point index.

        This function retrieves the electron-phonon matrix from the cache if available;
        otherwise, it reads the matrix from the lelph database and caches it. The matrix
        is transformed to maintain compatibility by selecting spin 0 and transposing the axes.

        Parameters
        ----------
        iq : int
            The index of the q-point for which the electron-phonon matrix is to be loaded.

        Returns
        -------
        ndarray
            The transformed electron-phonon matrix for the specified q-point.
        """
        if iq in self._eph_mat_cache:
            return self._eph_mat_cache[iq]

        _, eph_mat_iq = self.lelph_db.read_iq(iq, convention='standard')
        self._eph_mat_cache[iq] = eph_mat_iq[:, :, 0, :, :].transpose(1, 0, 3, 2)
        # Select spin 0 and transpose axes: (spin, mode, m, n) → (mode, m, n)
        # Original shape: (nb1, nb2, 2 spins, nm, nk) → we keep only spin 0 and swap bands for compatibility
        return eph_mat_iq[:, :, 0, :, :].transpose(1, 0, 3, 2)
    
    def _rotate_exciton_pair(self, iq, iQ):
        """
        Rotate exciton wavefunctions for a given phonon index and exciton Q index.

        This function computes the rotated exciton wavefunctions for the given q-point 
        and Q-point indices. It returns the wavefunctions corresponding to the 
        A^{S,Q}_{cvk} and A^{S,Q+q}_{cvk} states.

        Parameters
        ----------
        iq : int
            The index of the q-point.
        iQ : int
            The index of the Q-point.

        Returns
        -------
        tuple
            A tuple containing two rotated exciton wavefunctions:
            - Ak : ndarray
                The wavefunction for A^{S,Q}_{cvk}.
            - Akq : ndarray
                The wavefunction for A^{S,Q+q}_{cvk}.
        """
        # For q-point k
        ik_ibz, isym = self.kmap[self.qidx_in_kpts[iq]]
        Ak = self._rotate_exciton_wavefunction(ik_ibz, isym)

        # For q + Q point
        ik_ibz_qplusQ, isym_qplusQ = self.kmap[self.qplusQ_in_kpts[iQ, iq]]
        Akq = self._rotate_exciton_wavefunction(ik_ibz_qplusQ, isym_qplusQ)

        return Ak, Akq

    def _rotate_exciton_wavefunction(self, ik_ibz, isym):
        """
        Rotate exciton wavefunction from IBZ point using symmetry operation.

        Parameters
        ----------
        ik_ibz : int
            Index of the k-point in the IBZ.
        isym : int
            Index of the symmetry operation.

        Returns
        -------
        Ak : ndarray
            The rotated exciton wavefunction A^{S,RQ}_{cvk}, where RQ is the rotated Q-point.

        Notes
        -----
        The rotation is done as follows:
        A^{S,RQ}_{cvk} = \sum_{k'i'j'} \mathcal{U}^{Q}_{k'i'j'kij}(g)A^{S,Q}_{k'i'j'}
        """
        key = (ik_ibz, isym)
        # The same time a key pair appears we save time and do not call rotate_exc_wf
        if key in self._Ak_cache:
            return self._Ak_cache[key]

        is_sym_time_rev = isym >= self.symm_mats.shape[0] // (int(self.ele_time_rev) + 1)

        Ak = rotate_exc_wf(
            self.BS_wfcs[ik_ibz],
            self.sym_red[isym],
            self.kpts.data,
            self.wfdb.kpts_iBZ[ik_ibz],
            self.Dmats[isym],
            is_sym_time_rev,
            ktree=self.kpt_tree
        )

        self._Ak_cache[key] = Ak
        return Ak

    def _save_or_load_exph(self, ex_ph):
        """
        Save or load exciton-phonon matrix elements to/from a file.

        This function either saves the provided exciton-phonon matrix elements to a 
        file named 'Ex-ph.npy' or loads these elements from the file, depending on 
        the value of the 'save_files' attribute. If saving, the matrix elements are 
        converted to the same data type as `BS_wfcs` before saving. If loading, it 
        checks for the file's existence and raises an error if it is not found.

        Parameters
        ----------
        ex_ph : ndarray
            The exciton-phonon matrix elements to be saved or loaded.

        Returns
        -------
        ex_ph : ndarray
            The saved or loaded exciton-phonon matrix elements.

        Raises
        ------
        FileNotFoundError
            If loading is attempted and 'Ex-ph.npy' does not exist.
        """
        from time import time
        import os

        if self.save_files:
            print('Saving exciton-phonon matrix elements...')
            t0 = time()
            np.save('Ex-ph.npy', ex_ph.astype(self.BS_wfcs.dtype))
            self.timings['exph_io'] += time() - t0
            return ex_ph
        else:
            print('Loading exciton-phonon matrix elements...')
            t0 = time()
            if not os.path.exists('Ex-ph.npy'):
                raise FileNotFoundError("Cannot load 'Ex-ph.npy' - file does not exist.")
            ex_ph_loaded = np.load('Ex-ph.npy')
            self.timings['exph_io'] += time() - t0
            return ex_ph_loaded

