import warnings
from numba import njit, prange
import os
import numpy as np
from netCDF4 import Dataset
from yambopy.units import *
from yambopy.optical_properties.base_optical import BaseOpticalProperties
from yambopy.optical_properties.utils import (
    validate_path, setup_directories, safe_file_operation
)

from tqdm import tqdm
warnings.filterwarnings('ignore')

class ExcitonDipole(BaseOpticalProperties):
    def __init__(self, path=None, save='SAVE', latdb=None, wfdb=None, 
                 ydipdb=None, bands_range=None, BSE_dir='bse', 
                 DIP_dir='gw', save_files=True):
        """
        Initialize ExcitonDipole class.
        
        Parameters
        ----------
        path : str, optional
            Path to calculation directory. Defaults to current directory.
        save : str, optional
            SAVE directory name. Defaults to 'SAVE'.
        latdb : YamboLatticeDB, optional
            Pre-loaded lattice database.
        wfdb : YamboWFDB, optional
            Pre-loaded wavefunction database.
        ydipdb : YamboDipolesDB, optional
            Pre-loaded dipoles database.
        bands_range : list, optional
            Range of bands to load.
        BSE_dir : str, optional
            BSE directory name. Defaults to 'bse'.
        DIP_dir : str, optional
            Dipoles directory name. Defaults to 'gw'.
        save_files : bool, optional
            Whether to save files in .npy database. Defaults to True.
        """
        # Initialize base class
        super().__init__(path=path, save=save, latdb=latdb, wfdb=wfdb, 
                        bands_range=bands_range, BSE_dir=BSE_dir, save_files=save_files)
        
        # Setup additional directories
        self._setup_directories(DIP_dir=DIP_dir)
        
        # Store dipoles database
        self.ydipdb = ydipdb
        
        # Read all necessary databases
        self.read(latdb=latdb, wfdb=wfdb, ydipdb=ydipdb, bands_range=bands_range)

    def read(self, latdb=None, wfdb=None, ydipdb=None, bands_range=None):
        """
        Read all necessary databases for exciton-dipole calculations.
        
        Parameters
        ----------
        latdb : YamboLatticeDB, optional
            Pre-loaded lattice database.
        wfdb : YamboWFDB, optional
            Pre-loaded wavefunction database.
        ydipdb : YamboDipolesDB, optional
            Pre-loaded dipoles database.
        bands_range : list, optional
            Range of bands to load.
        """
        # Read common databases using base class method
        self.read_common_databases(latdb=latdb, wfdb=wfdb, bands_range=bands_range)
        
        # Read dipoles database
        self._read_dipoles_db(ydipdb, dip_dir=self.DIP_dir, bands_range=bands_range)

    def compute(self):
        """
        Main computation method - computes exciton-photon coupling matrix elements.
        
        Returns
        -------
        np.ndarray
            Exciton-dipole matrix elements.
        """
        return self.compute_Exdipole()
    
    def compute_Exdipole(self):
        """
        Compute exciton-photon coupling matrix elements.
        
        Returns
        -------
        np.ndarray
            Exciton-dipole matrix elements.
        """
        from time import time
        
        start_time = time()
        print('Computing Exciton-photon matrix elements')
        
        # Compute exciton-dipole matrix elements
        self.ex_dip = self.exe_dipoles(
            self.ele_dips, self.BS_wfcs[0],
            self.kmap, self.symm_mats, self.ele_time_rev
        )
        
        # Save results if requested
        if self.save_files:
            np.save('Ex-dipole', self.ex_dip)
        
        computation_time = time() - start_time
        print(f'Exciton-dipole computation completed in {computation_time:.4f} s')
        print('*' * 60)
        
        return self.ex_dip

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