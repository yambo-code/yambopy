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
import numpy as np
from netCDF4 import Dataset
from yambopy.units import *
from yambopy.optical_properties.base_optical import BaseOpticalProperties
from yambopy.optical_properties.ex_dipole import ExcitonDipole
from yambopy.optical_properties.ex_phonon import ExcitonPhonon
from yambopy.optical_properties.utils import (
    read_lelph_database, create_progress_bar
)

from tqdm import tqdm
warnings.filterwarnings('ignore')

class Luminescence(BaseOpticalProperties):
    def __init__(self, path=None, save='SAVE', lelph_db=None, latdb=None, wfdb=None, 
                 ydipdb=None, bands_range=None, BSE_dir='bse', LELPH_dir='lelph', 
                 DIP_dir='gw', save_files=True):
        """
        Initialize the Luminescence class.

        Parameters
        ----------
        path : str, optional
            Path to calculation directory. Defaults to current directory.
        save : str, optional
            SAVE directory name. Defaults to 'SAVE'.
        lelph_db : LetzElphElectronPhononDB, optional
            Pre-loaded electron-phonon database.
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
        LELPH_dir : str, optional
            LELPH directory name. Defaults to 'lelph'.
        DIP_dir : str, optional
            Dipoles directory name. Defaults to 'gw'.
        save_files : bool, optional
            Whether to save files in .npy database. Defaults to True.
        """
        # Initialize base class
        super().__init__(path=path, save=save, latdb=latdb, wfdb=wfdb, 
                        bands_range=bands_range, BSE_dir=BSE_dir, save_files=save_files)
        
        # Setup additional directories
        self._setup_directories(LELPH_dir=LELPH_dir, DIP_dir=DIP_dir)
        
        # Store specific parameters
        self.lelph_db = lelph_db
        self.ydipdb = ydipdb
        
        # Read all necessary databases
        self.read(lelph_db=lelph_db, latdb=latdb, wfdb=wfdb, 
                  ydipdb=ydipdb, bands_range=bands_range)

    def read(self, lelph_db=None, latdb=None, wfdb=None, ydipdb=None, bands_range=None):
        """
        Read all necessary databases for luminescence calculations.
        
        Parameters
        ----------
        lelph_db : LetzElphElectronPhononDB, optional
            Pre-loaded electron-phonon database.
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
        
        # Read LetzElPhC database
        self.lelph_db = read_lelph_database(self.LELPH_dir, lelph_db)
        self.kpts = self.lelph_db.kpoints
        self.qpts = self.lelph_db.qpoints
        self.elph_bnds_range = self.lelph_db.bands
        self.ph_freq = self.lelph_db.ph_energies / ha2ev  # Convert to Hartree
        
        # Read dipoles database
        self._read_dipoles_db(ydipdb, dip_dir='gw', bands_range=bands_range)

    def compute(self):
        """
        Main computation method - computes luminescence properties.
        
        Returns
        -------
        dict
            Luminescence calculation results.
        """
        return self.compute_luminescence()
    
    def compute_luminescence(self):
        """
        Compute luminescence properties using exciton-dipole and exciton-phonon coupling.
        
        Returns
        -------
        dict
            Dictionary containing luminescence results.
        """
        from time import time
        
        start_time = time()
        print('Computing Luminescence properties')
        
        # Initialize ExcitonDipole and ExcitonPhonon objects using existing data
        ex_dipole = ExcitonDipole(
            path=self.path, latdb=self.ydb, wfdb=self.wfdb, 
            ydipdb=self.ydipdb, bands_range=self.bands_range
        )
        
        ex_phonon = ExcitonPhonon(
            path=self.path, lelph_db=self.lelph_db, latdb=self.ydb, 
            wfdb=self.wfdb, ydipdb=self.ydipdb, bands_range=self.bands_range
        )
        
        # Compute exciton-dipole matrix elements
        dipole_matrix = ex_dipole.compute()
        
        # Placeholder for luminescence calculation
        # The actual implementation would combine dipole and phonon matrix elements
        print("Luminescence computation method needs full implementation")
        print("This would combine exciton-dipole and exciton-phonon matrix elements")
        
        computation_time = time() - start_time
        print(f'Luminescence computation completed in {computation_time:.4f} s')
        print('*' * 60)
        
        return {
            'dipole_matrix': dipole_matrix,
            'computation_time': computation_time
        }


    
    def compute_luminescence_spectrum(self, 
                                     ome_range,
                                     temp=20,
                                     broadening=0.00124, 
                                     npol=2, 
                                     ph_thr=1e-9):   
        """
        Compute luminescence spectrum intensities.

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
        tuple
            (ome_range, luminescence_intensities) arrays.
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