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
                 DIP_dir='gw', save_files=True, field_dir = [1,1,1]):
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
        field_dir : list, optional
            Direction of the electric field. Defaults to [1,1,1].            
        """
        # Initialize base class
        super().__init__(path=path, save=save, latdb=latdb, wfdb=wfdb, 
                        bands_range=bands_range, BSE_dir=BSE_dir, save_files=save_files)
        
        # Setup additional directories
        self._setup_directories(LELPH_dir=LELPH_dir, DIP_dir=DIP_dir)
        
        # Store specific parameters
        self.lelph_db = lelph_db
        self.ydipdb = ydipdb
        self.field_dir = field_dir
        
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
        
        # Read LetzElPhC database (gracefully handles missing files)
        self.lelph_db = read_lelph_database(self.LELPH_dir, lelph_db)
        if self.lelph_db:
            self.kpts = self.lelph_db.kpoints
            self.qpts = self.lelph_db.qpoints
            self.elph_bnds_range = self.lelph_db.bands
            self.ph_freq = self.lelph_db.ph_energies / ha2ev  # Convert to Hartree
        else:
            # Fall back to geometry manager k-points
            self.kpts = self.red_kpoints
            self.qpts = None
            self.elph_bnds_range = None
            self.ph_freq = None
        
        # Read dipoles database
        self._read_dipoles_db(ydipdb = ydipdb, DIP_dir=self.DIP_dir, bands_range=bands_range)
        
        # Build k-point tree and find q-point indices (needed for luminescence)
        from yambopy.kpoints import build_ktree, find_kpt
        print('Building kD-tree for kpoints')
        self.kpt_tree = build_ktree(self.kpts)
        if self.qpts is not None:
            self.qidx_in_kpts = find_kpt(self.kpt_tree, self.qpts)
        else:
            self.qidx_in_kpts = None  # No q-points available

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
        Compute luminescence matrix elements (exciton-dipole and exciton-phonon).
        
        Returns
        -------
        dict
            Dictionary containing computed matrix elements.
        """
        from time import time
        
        start_time = time()
        print('Computing Luminescence matrix elements')
        
        # Compute exciton-dipole matrix elements
        try: 
            if hasattr(self, 'ex_dip'):
                print('exciton-photon matrix elements found')
            else:             
                print('Computing exciton-photon matrix elements')
                ex_dipole = ExcitonDipole(
                    path=self.path, latdb=self.ydb, wfdb=self.wfdb, 
                    ydipdb=self.ydipdb, bands_range=self.bands_range, 
                    BSE_dir=self.BSE_dir, DIP_dir=self.BSE_dir, save_files=self.save_files
                )
                self.ex_dip = ex_dipole.compute()
        except Exception as e:
            raise IOError(f'Cannot compute exciton-photon matrix elements: {e}')

        # Compute exciton-phonon matrix elements
        try: 
            if hasattr(self, 'ex_ph'):
                print('exciton-phonon matrix elements found')
            else:
                print('Computing exciton-phonon matrix elements')
                ex_phonon = ExcitonPhonon(
                    path=self.path, lelph_db=self.lelph_db, latdb=self.ydb, 
                    wfdb=self.wfdb, ydipdb=self.ydipdb, bands_range=self.bands_range,
                    BSE_dir=self.BSE_dir, LELPH_dir=self.LELPH_dir, DIP_dir=self.BSE_dir,
                )
                ex_phonon.compute(gamma_only=True)
                self.ex_ph = ex_phonon.ex_ph[0]
                
        except Exception as e:
            raise IOError(f'Cannot compute exciton-phonon matrix elements: {e}')
        
        computation_time = time() - start_time
        print(f'Luminescence matrix elements computed in {computation_time:.4f} s')
        print('*' * 60)
        
        return {
            'ex_dip': self.ex_dip,
            'ex_ph': self.ex_ph,
            'computation_time': computation_time
        }


    
    def compute_luminescence_spectrum(self, 
                                     ome_range,
                                     temp=20,
                                     broadening=0.00124, 
                                     npol=2, 
                                     ph_thr=1e-9,
                                     direct_pl = False,
                                     ):   
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
                print('exciton-photon matrix elements found')
            else:             
                print('Computing exciton-photon matrix elements')
                ExDipole = ExcitonDipole(
                    path=self.path, latdb=self.ydb, wfdb=self.wfdb, 
                    ydipdb=self.ydipdb, bands_range=self.bands_range, BSE_dir = self.BSE_dir,
                    DIP_dir=self.BSE_dir, save_files=self.save_files
                )
                self.ex_dip = ExDipole.compute()
        except Exception as e:
            raise IOError(f'Cannot compute exciton-photon matrix elements: {e}')

        try: 
            if hasattr(self, 'ex_ph'):
                print('exciton-phonon matrix elements found')
            else:
                print('Computing exciton-phonon matrix elements')
                ExPhonon = ExcitonPhonon(
                    path=self.path, lelph_db=self.lelph_db, latdb=self.ydb, 
                    wfdb=self.wfdb, bands_range=self.bands_range,
                    BSE_dir=self.BSE_dir, LELPH_dir=self.LELPH_dir, DIP_dir=self.BSE_dir,
                    save_files=self.save_files
                )
                ExPhonon.compute(gamma_only=True)
                self.ex_ph = ExPhonon.ex_ph[0]
                
        except Exception as e:
            raise IOError(f'Cannot compute exciton-phonon matrix elements: {e}')
        for i in tqdm(range(ome_range.shape[0]), desc="Luminescence "):
            inte_tmp = compute_luminescence_per_freq(ome_range[i], self.ph_freq, exe_ene, \
                Exe_min, self.ex_dip, self.ex_ph, npol=npol, ph_thr=ph_thr,broadening=broadening, temp=temp)
            if direct_pl == True :
                inte_direct = compute_dir_luminescence_per_freq(ome_range[i], self.ph_freq, exe_ene, \
                    Exe_min, self.ex_dip, self.ex_ph, npol=npol, ph_thr=ph_thr,broadening=broadening,temp=temp)
                inte_tmp += inte_direct  # Add direct PL contribution
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
                        temp=20.0,
                        broadening=0.00124,
                        npol=2,
                        ph_thr=1e-9,
                        ):
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
    # Extract shape information
    Nqpts = ex_ph.shape[0]
    nmode = ex_ph.shape[1]
    nbnd_i = ex_ph.shape[2]
    nbnd_f = ex_ph.shape[3]
    
    # Convert units and compute constants
    broadening_Ha = broadening / 27.211 / 2.0
    ome_light_Ha = ome_light / 27.211
    KbT = 3.1726919127302026e-06 * temp  ## Ha
    
    # Compute Boltzmann factors
    bolt_man_fac = np.exp(-(ex_ene - exe_low_energy) / KbT)  ##(iq,nexe)
        
    # Initialize output array for parallel reduction
    sum_array = np.zeros(Nqpts, dtype=np.float64)
    
    for iq in prange(Nqpts):
        local_sum = 0.0
        for iv in range(nmode):  # mu
            # Compute frequency factor
            ome_fac = ome_light_Ha * (ome_light_Ha - 2.0 * ph_freq[iq, iv])**2
            
            # Compute Bose-Einstein factors
            if ph_freq[iq, iv] < ph_thr:
                bose_ph_fac = 1.0
                bos_occ = 1.0
                # Note: Warning removed - not compatible with numba
            else:
                bos_occ = 1.0 / (np.exp(ph_freq[iq, iv] / KbT) - 1.0)
                bose_ph_fac = 1.0 + bos_occ
            
            # Compute final state energy
            E_f_omega = ex_ene[iq, :] - ph_freq[iq, iv]
            
            # Initialize scattering matrix
            Tmu = np.zeros((npol, nbnd_f), dtype=np.complex128)  # D*G
            
            ## Compute scattering matrix
            for ipol in range(npol): 
                for ii in range(nbnd_i):  # lambda
                    # Compute denominator with complex broadening
                    denom = ex_ene[0, ii] - E_f_omega + 1j * broadening_Ha
                    # Compute contribution
                    contrib = np.conj(ex_ph[iq, iv, ii, :]) * ex_dip[ipol, ii] / denom
                    Tmu[ipol, :] = Tmu[ipol, :] + contrib
            
            ## Compute Gamma_mu: abs and sum over initial states and pols
            Tmu_abs_sq = np.abs(Tmu)**2
            Tmu_sum = np.sum(Tmu_abs_sq, axis=0)
            
            # Compute denominator for Gamma_mu
            denom_gamma = E_f_omega * ((ome_light_Ha - E_f_omega)**2 + broadening_Ha**2)
            
            Gamma_mu = bose_ph_fac * Tmu_sum * ome_fac * bolt_man_fac[iq, :] / denom_gamma
            
            # Compute direct contribution
            # Use squared denominator to match Lorentzian form
            local_sum = local_sum + np.sum(Gamma_mu).real
        
        sum_array[iq] = local_sum
    
    sum_out = np.sum(sum_array)
    return sum_out * broadening_Ha / np.pi / Nqpts

@njit(cache=True, nogil=True, parallel=True)
def compute_dir_luminescence_per_freq(ome_light,
                        ph_freq,
                        ex_ene,
                        exe_low_energy,
                        ex_dip,
                        ex_ph,
                        temp=20.0,
                        broadening=0.00124,
                        npol=2,
                        ph_thr=1e-9,
                        ):
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
    # Extract shape information
    Nqpts = ex_ph.shape[0]
    nmode = ex_ph.shape[1]
    nbnd_i = ex_ph.shape[2]
    nbnd_f = ex_ph.shape[3]
    
    # Convert units and compute constants
    broadening_Ha = broadening / 27.211 / 2.0
    ome_light_Ha = ome_light / 27.211
    KbT = 3.1726919127302026e-06 * temp  ## Ha
    
    # Compute Boltzmann factors
    bolt_man_fac = np.exp(-(ex_ene - exe_low_energy) / KbT)  ##(iq,nexe)
    
    # Pre-compute ex_dip absolute square sum (independent of q-point loop)
    ex_dip_pol_sum = np.sum(ex_dip,axis=0)
    ex_dip_sum = np.abs(ex_dip_pol_sum)**2  # Sum over all polarizations and bands
    
    # Step 1: Compute R_lambda by summing over all q, mu, beta
    # R_lambda has shape (nbnd_f,) representing the scattering rate for each final state lambda
    R_lambda_array = np.zeros((Nqpts, nbnd_f), dtype=np.complex128)
    
    for iq in prange(Nqpts):
        R_lambda_local = np.zeros(nbnd_f, dtype=np.complex128)
        
        for iv in range(nmode):  # mu (phonon modes)
            if ph_freq[iq, iv] < ph_thr:
                bose_ph_fac = 1.0
                bos_occ = 1.0
            else:
                bos_occ = 1.0 / (np.exp(ph_freq[iq, iv] / KbT) - 1.0)
                bose_ph_fac = 1.0 + bos_occ
                     
            reg = broadening_Ha**2  # Use broadening as regularization 
                   
            ## Sum over beta (initial bands) and polarizations
            for ipol in range(npol): 
                for ii in range(nbnd_i):  # beta (initial exciton states)
                    # |G|^2 term
                    G_squared = np.abs(ex_ph[iq, iv, ii, :])**2
                    
                    # Compute denominators: (E_lambda - E_beta +/- omega_q,mu)^2
                    denom1 = (ex_ene[0, :] - ex_ene[iq, ii] + ph_freq[iq, iv])**2 + reg
                    denom2 = (ex_ene[0, :] - ex_ene[iq, ii] - ph_freq[iq, iv])**2 + reg
                    
                    # Bose factors: (n+1)/denom1 + n/denom2
                    term1 = bose_ph_fac / denom1
                    term2 = bos_occ / denom2
                    
                    # Accumulate R_lambda
                    R_lambda_local[:] += G_squared * (term1 + term2)
        
        R_lambda_array[iq, :] = R_lambda_local
    
    # Step 2: Sum R_lambda over all q-points
    R_lambda = np.sum(R_lambda_array, axis=0)
    
    # Step 3: Compute I^PL(omega) = sum over lambda
    # Denominator: (omega - E_lambda + i*eta)
    denom_direct = (ome_light_Ha - ex_ene[0, :])**2 + broadening_Ha**2
    
    # Direct PL intensity: |T_lambda|^2 * omega^3 * (1 - R_lambda) / (omega - E_lambda + i*eta) * exp(...)
    Direct_lambda = (ex_dip_sum * ome_light_Ha**3 * (1.0 - R_lambda) / denom_direct * bolt_man_fac[0, :]).real
    
    # Sum over all lambda states
    sum_out = np.sum(Direct_lambda).real
    return sum_out * broadening_Ha / np.pi / Nqpts