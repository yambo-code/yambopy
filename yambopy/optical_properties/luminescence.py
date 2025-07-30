#
# License-Identifier: GPL
#
# Copyright (C) 2024 The Yambo Team
#
# Authors: RR, MN
#
# This file is part of the yambopy project
#
import warnings
from numba import njit, prange
import os
from netCDF4 import Dataset
from yambopy.letzelphc_interface.lelphcdb import LetzElphElectronPhononDB
from yambopy.dbs.latticedb import YamboLatticeDB
from yambopy.dbs.wfdb import YamboWFDB
from yambopy.dbs.excitondb import YamboExcitonDB
from yambopy.dbs.dipolesdb import YamboDipolesDB
from yambopy.units import *
from yambopy.optical_properties.ex_dipole import ExcitonDipole
from yambopy.optical_properties.ex_phonon import ExcitonPhonon
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

class Luminescence(object):
    def __init__(self, path=None, save='SAVE', lelph_db=None, latdb=None, wfdb=None, \
                 ydipdb=None, bands_range=[], BSE_dir='bse', LELPH_dir='lelph', \
                 DIP_dir='gw',save_files=True):
        """
        Initialize the Luminescence class.

        Args:
            path (str, optional): Path to the directory containing the wavefunction files. Defaults to the current directory.
            save (str, optional): Subdirectory containing the wavefunction files. Defaults to 'SAVE'.
            lelph_db (LetzElphElectronPhononDB, optional): Electron-Phonon coupling database. Defaults to None.
            latdb (YamboLatticeDB, optional): Yambo Lattice database. Defaults to None.
            wfdb (YamboWFDB, optional): Yambo wavefunction database. Defaults to None.
            ydipdb (YamboDipolesDB, optional): Yambo Dipoles database. Defaults to None.
            bands_range (list, optional): Range of bands to load. Defaults to all bands. Python indexing. Right one is excluded.
            BSE_dir (str, optional): Subdirectory containing the BSE output. Defaults to 'bse'.
            LELPH_dir (str, optional): Subdirectory containing the LELPH output. Defaults to 'lelph'.
            DIP_dir (str, optional): Subdirectory containing the dipoles output. Defaults to 'gw'.
            save_files (bool, optional): Whether to save files in .npy database. Defaults to True.
        """
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
        """
        Read Yambo databases and wavefunction files.

        Reads the Yambo Lattice database, Yambo wavefunction database, Yambo Exciton
        database, Yambo Dipoles database, and LetzElPhC database. Additionally, reads
        the wavefunction files and computes the electron-phonon coupling matrix elements
        using the save_Dmat method of the Yambo wavefunction database.

        Parameters
        ----------
        lelph_db : LetzElphElectronPhononDB
            Electron-Phonon coupling database.
        latdb : YamboLatticeDB
            Yambo Lattice database.
        wfdb : YamboWFDB
            Yambo wavefunction database.
        ydipdb : YamboDipolesDB
            Yambo Dipoles database.
        bands_range : list
            Range of bands to load. Defaults to all bands. Python indexing. Right one is excluded.
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
        Dmat_data = Dataset(f'{self.LELPH_dir}/ndb.Dmats', 'r')
        Dmats = Dmat_data['Dmats'][:, :, 0, start_bnd_idx:end_bnd,
                                    start_bnd_idx:end_bnd, :].data
        Dmats = Dmats[...,0] + 1j * Dmats[..., 1]
        Dmats = Dmats.astype(self.wfdb.wf.dtype)
        self.Dmats    = Dmats #self.wfdb.save_Dmat
        Dmat_data.close()
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
        """
        Reads YamboExcitonDB objects for all q-points, extracting the number of bands involved in the BSE,
        the eigenenergies of the excitons, their wavefunctions, and the Q-point of each BSE calculation.
        
        Parameters
        ----------
        BSE_dir : str
            Path to the directory containing the BSE calculations.
        
        Returns
        -------
        bs_bands : list
            List of the number of bands involved in the BSE for each q-point.
        BS_eigs : list
            List of the eigenenergies of the excitons for each q-point.
        BS_wfcs : list
            List of the exciton wavefunctions for each q-point.
        excQpt : list
            List of the Q-point of each BSE calculation.
        """
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
    
    def compute_luminescence(self, 
                             ome_range,
                             temp = 20,
                             broadening = 0.00124, 
                             npol = 2, 
                             ph_thr = 1e-9                             
                             ):   
        """
        Compute luminescence intensities.

        Parameters
        ----------
        ome_range : tuple (start, end, num)
            Range of energies to compute luminescence intensities.
        temp : float, optional
            Temperature in Kelvin. Default is 20.
        broadening : float, optional
            Broadening of the luminescence peaks. Default is 0.00124 Ha.
        npol : int, optional
            Number of polarizations. Default is 2.
        ph_thr : float, optional
            Threshold for phonon frequencies. Default is 1e-9 Ry.

        Returns
        -------
        ome_range : array
            Array of energies.
        self_inten : array
            Array of luminescence intensities.
        """
        ome_range = np.linspace(ome_range[0], ome_range[1], num=ome_range[2])
        exe_ene = self.BS_eigs[self.kmap[self.qidx_in_kpts, 0], :]
        self_inten = []

        Exe_min = np.min(exe_ene)
        print(f'Minimum energy of the exciton is        : {Exe_min*ha2ev:.4f} eV')
        print(
            f'Minimum energy of the Direct exciton is : {min(self.BS_eigs[0]*ha2ev):.4f} eV'
        )
        print('Computing luminescence intensities ...')
        try: 
            if hasattr(self, 'ex_dip'):
                print('exciton-photon matrix elements founds')
                #self.ex_dip = np.load(f'ex_dip')
            else:             
                print('Computing exciton-photon matrix elements')
                ExDipole = ExcitonDipole(self.path, self.SAVE_dir, self.latdb, self.wfdb, \
                self.ydipdb, self.bands_range, self.BSE_dir, \
                self.DIP_dir,self.save_files)
                ExDipole.compute_Exdipole()
                self.ex_dip = ExDipole.ex_dip
        except Exception as e:
            raise IOError(f'Cannot compute exciton-photon matrix elements')

        try: 
            if hasattr(self, 'ex_ph'):
                print('exciton-phonon matrix elements founds')
            else:
                print('Computing exciton-phonon matrix elements')
                ExPhonon =  ExcitonPhonon(self.path, self.SAVE_dir, self.lelph_db, self.latdb, self.wfdb, \
                self.ydipdb, self.bands_range, self.BSE_dir, self.LELPH_dir, \
                self.DIP_dir,self.save_files)
                ExPhonon.compute_Exph(gamma_only = True)
                self.ex_ph  = ExPhonon.ex_ph[0]
                
        except Exception as e:
            raise IOError(f'Cannot compute exciton-phonon matrix elements')
        for i in tqdm(range(ome_range.shape[0]), desc="Luminescence "):
            inte_tmp = compute_luminescence_per_freq(ome_range[i], self.ph_freq, exe_ene, \
                Exe_min, self.ex_dip, self.ex_ph, npol=npol, ph_thr=ph_thr,broadening=broadening, temp=temp)
            self_inten.append(inte_tmp)
        ## save intensties
        if self.save_files:
            np.savetxt('luminescence_intensities.dat', np.c_[ome_range,
                                                        np.array(self_inten)].real)
            np.save('Intensties_self', np.array(self_inten))
        return ome_range,np.array(self_inten).real
    
@njit(cache=True, nogil=True, parallel=True)
def compute_luminescence_per_freq(ome_light,
                        ph_freq,
                        ex_ene,
                        exe_low_energy,
                        ex_dip,
                        ex_ph,
                        temp=20,
                        broadening=0.00124,
                        npol=2,
                        ph_thr = 1e-9):
    ## We need exciton dipoles for light emission (<0|r|S>)
    ## and exciton phonon matrix elements for phonon absorption <S',Q|dV_Q|S,0>
    ## energy of the lowest energy energy exe_low_energy
    """
    Compute the luminescence intensity per frequency.

    Parameters
    ----------
    ome_light : float
        Photon energy in eV
    ph_freq : array_like
        Phonon frequencies in Ha
    ex_ene : array_like
        Exciton energies in Ha
    exe_low_energy : float
        Lowest energy of the exciton in Ha
    ex_dip : array_like
        Exciton-photon matrix elements in a.u.
    ex_ph : array_like
        Exciton-phonon matrix elements in a.u.
    temp : float, optional
        Temperature in Kelvin
    broadening : float, optional
        Broadening in eV
    npol : int, optional
        Number of polarizations
    ph_thr : float, optional
        Threshold for negative frequencies

    Returns
    -------
    luminescence : float
        Luminescence intensity per frequency
    """
    Nqpts, nmode, nbnd_i, nbnd_f = ex_ph.shape
    broadening = (broadening / 27.211 / 2)
    ome_light_Ha = (ome_light / 27.211)
    KbT = (3.1726919127302026e-06 * temp)  ## Ha
    bolt_man_fac = -(ex_ene - exe_low_energy) / KbT
    bolt_man_fac = np.exp(bolt_man_fac)  ##(iq,nexe)
    sum_out = 0.0
    for iq in prange(Nqpts):
        for iv in range(nmode):
            ome_fac = ome_light_Ha * (ome_light_Ha + 2 * ph_freq[iq, iv])**2
            if ph_freq[iq, iv] < ph_thr:
                bose_ph_fac = 1.0
                Warning('Negative frequencies set to zero')
            else:
                bose_ph_fac = 1 + 1.0 / (np.exp(ph_freq[iq, iv] / KbT) - 1.0)
            E_f_omega = ex_ene[iq, :] - ph_freq[iq, iv]
            Tmu = np.zeros((npol, nbnd_f), dtype=np.complex64)  # D*G
            ## compute scattering matrix
            for ipol in range(npol):
                for ii in range(nbnd_i):
                    Tmu[ipol,:] = Tmu[ipol,:] + np.conj(ex_ph[iq,iv,ii,:]) * ex_dip[ipol,ii] \
                        /(ex_ene[0,ii] - E_f_omega + (1j*broadening)).astype(np.complex64)
            ## abs and sum over initial states and pols
            Gamma_mu = bose_ph_fac * np.sum(np.abs(Tmu)**2,axis=0) * ome_fac * bolt_man_fac[iq,:] \
                        /E_f_omega/((ome_light_Ha-E_f_omega)**2 + broadening
**2)
            sum_out = sum_out + np.sum(Gamma_mu)
    return sum_out * broadening/ np.pi / Nqpts