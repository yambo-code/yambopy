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
import numpy as np
import os
from netCDF4 import Dataset
from yambopy.letzelphc_interface.lelphcdb import LetzElphElectronPhononDB
from yambopy.units import *
from yambopy.bse.rotate_excitonwf import rotate_exc_wf
from yambopy.optical_properties.base_optical import BaseOpticalProperties
from yambopy.optical_properties.utils import (
    read_lelph_database, compute_symmetry_matrices
)
from tqdm import tqdm
import re

warnings.filterwarnings('ignore')

class ExcitonGroupTheory(BaseOpticalProperties):
    """
    This class performs group theory analysis of exciton states.
    
    It analyzes the irreducible representations of exciton states under the 
    little group of the exciton momentum, providing insight into the symmetry
    properties of excitonic states.
    
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
    bands_range : list or tuple, optional
        The range of bands for which the analysis will be performed.
        Default is all bands.
    BSE_dir : str, optional
        The name of the folder which contains the BSE calculation. Default is 'bse'.
    LELPH_dir : str, optional
        The name of the folder which contains the electron-phonon matrix elements.
        Default is 'lelph'.
    read_symm_from_ns_db_file : bool, optional
        If True, will read symmetry matrices from ns.db1 file else from ndb.elph.
        Default is False.
    
    Attributes
    ----------
    SAVE_dir : str
        The path of the SAVE folder.
    BSE_dir : str
        The path of the BSE folder.
    LELPH_dir : str
        The path of the folder which contains the electron-phonon matrix elements.
    latdb : YamboLatticeDB
        The YamboLatticeDB object which contains the lattice information.
    lelph_db : LetzElphElectronPhononDB
        The LetzElphElectronPhononDB object which contains the electron-phonon matrix
        elements.
    wfdb : YamboWFDB
        The YamboWFDB object which contains the wavefunction information.
    """
    
    def __init__(self, path=None, save='SAVE', lelph_db=None, latdb=None, wfdb=None, 
                 bands_range=None, BSE_dir='bse', LELPH_dir='lelph', 
                 read_symm_from_ns_db_file=True):
        """
        Initialize ExcitonGroupTheory class.
        
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
        bands_range : list, optional
            Range of bands to load.
        BSE_dir : str, optional
            BSE directory name. Defaults to 'bse'.
        LELPH_dir : str, optional
            LELPH directory name. Defaults to 'lelph'.
        read_symm_from_ns_db_file : bool, optional
            Whether to read symmetry from ns.db1 file. Defaults to True.
        """
        # Initialize base class
        super().__init__(path=path, save=save, latdb=latdb, wfdb=wfdb, 
                        bands_range=bands_range, BSE_dir=BSE_dir)
        
        # Setup additional directories
        self._setup_directories(LELPH_dir=LELPH_dir)
        
        # Store specific parameters
        self.lelph_db = lelph_db
        self.read_symm_from_ns_db_file = read_symm_from_ns_db_file
        
        # Read all necessary databases
        self.read(lelph_db=lelph_db, latdb=latdb, wfdb=wfdb, bands_range=bands_range)

    def read(self, lelph_db=None, latdb=None, wfdb=None, bands_range=None):
        """
        Read all necessary databases for group theory analysis.
        
        Parameters
        ----------
        lelph_db : LetzElphElectronPhononDB, optional
            Pre-loaded electron-phonon database.
        latdb : YamboLatticeDB, optional
            Pre-loaded lattice database.
        wfdb : YamboWFDB, optional
            Pre-loaded wavefunction database.
        bands_range : list, optional
            Range of bands to load.
        """
        # Read common databases using base class method
        self.read_common_databases(latdb=latdb, wfdb=wfdb, bands_range=bands_range)
        
        # Read LetzElPhC database
        self.lelph_db = read_lelph_database(self.LELPH_dir, lelph_db)
        self.qpts = self.lelph_db.qpoints
        self.elph_bnds_range = self.lelph_db.bands

        # Read D-matrices
        if bands_range:
            self.Dmats = self.wfdb.Dmat()[:,:,0,:,:]
        else:
            self.Dmats = self.wfdb.Dmat()[:,:,0,:,:]

        # Handle symmetry matrices and kmap specific to group theory
        self._setup_symmetry_data()

    def _setup_symmetry_data(self):
        """Setup symmetry data specific to group theory analysis."""
        ndb_lelph_fname = os.path.join(self.LELPH_dir, 'ndb.elph')
        
        if not self.read_symm_from_ns_db_file:
            try:
                elph_file = Dataset(ndb_lelph_fname, 'r')
                self.kpts = elph_file['kpoints'][...].data
                self.kmap = elph_file['kmap'][...].data
                self.symm_mats = elph_file['symmetry_matrices'][...].data
                self.ele_time_rev = elph_file['time_reversal_phonon'][...].data
                self.frac_trans = elph_file['fractional_translation'][...].data
                # Convert to crystal coordinates
                self.frac_trans = np.einsum('ij,nj->ni', self.blat_vecs.T, self.frac_trans)
                elph_file.close()
            except Exception as e:
                print(f"Warning: Could not read symmetry from elph file: {e}")
                self.frac_trans = np.zeros((self.symm_mats.shape[0], 3))
                # Fallback to lattice database kmap
                self.kmap = self.kmap
        else:
            self.frac_trans = np.zeros((self.symm_mats.shape[0], 3))
            # Use lattice database kmap and kpoints from lelph_db
            self.kpts = self.lelph_db.kpoints

        # Use existing k-point tree from base class
        if hasattr(self.wfdb, 'ktree'):
            self.kpt_tree = self.wfdb.ktree
        else:
            self._build_kpoint_tree(self.kpts)
        
        # Compute symmetry matrices in reduced coordinates using utility function
        self.sym_red = compute_symmetry_matrices(self.symm_mats, self.blat_vecs, self.lat_vecs)

        # Construct kpts_iBZ (following original algorithm exactly)
        self.kpts_iBZ = np.zeros((len(np.unique(self.kmap[:, 0])), 3))
        for i in range(self.kmap.shape[0]):
            ik_ibz, isym = self.kmap[i]
            if isym == 0:
                self.kpts_iBZ[ik_ibz, :] = self.kpts[i]

    def read_excdb_single(self, BSE_dir, iQ, nstates):
        """
        Read yambo exciton database for a specific Q-point.

        Parameters
        ----------
        BSE_dir : str
            The directory containing the BSE calculation data.
        iQ : int
            The Q-point index (1-based indexing as in Yambo).
        nstates : int
            Number of exciton states to read.

        Returns
        -------
        tuple
            (bands_range, BS_eigs, BS_wfcs) for the specific Q-point.
        """
        from yambopy.dbs.excitondb import YamboExcitonDB
        
        try:
            bse_db_iq = YamboExcitonDB.from_db_file(self.ydb, folder=BSE_dir,
                                                   filename=f'ndb.BS_diago_Q{iQ+1}')
        except Exception as e:
            raise IOError(f'Cannot read ndb.BS_diago_Q{iQ+1} file: {e}')
            
        bands_range = bse_db_iq.nbands
        BS_eigs = bse_db_iq.eigenvalues[:nstates]
        BS_wfcs = bse_db_iq.get_Akcv()[:nstates]
        
        # Convert to Hartree units
        BS_eigs = BS_eigs / ha2ev
        
        return bands_range, BS_eigs, BS_wfcs

    def compute(self):
        """
        Main computation method - placeholder for group theory analysis.
        
        Returns
        -------
        dict
            Results of group theory analysis.
        """
        print("ExcitonGroupTheory compute method called.")
        print("Use analyze_exciton_symmetry() method for specific Q-point analysis.")
        #For now it is a dummy method
        return {}

    def analyze_exciton_symmetry(self, iQ, nstates, degen_thres=0.001):
        """
        Perform group theory analysis for exciton states at a given Q-point.
        This implementation follows the algorithm in exe_rep_program.py exactly.

        Parameters
        ----------
        iQ : int
            The Q-point index (1-based indexing as in Yambo).
        nstates : int
            Number of exciton states to analyze.
        degen_thres : float, optional
            Degeneracy threshold in eV. Default is 0.001 eV.

        Returns
        -------
        results : dict
            Dictionary containing the analysis results including:
            - 'little_group': Little group symmetries
            - 'point_group_label': Point group label
            - 'unique_energies': Unique energy levels
            - 'degeneracies': Degeneracy of each level
            - 'irrep_decomposition': Irreducible representation decomposition
        """
        print('Reading BSE eigen vectors')
        bands_range, BS_eigs, BS_wfcs = self.read_excdb_single(self.BSE_dir, iQ-1, nstates)
        
        # Convert energies to eV for analysis (following original algorithm exactly)
        BS_eigs_eV = BS_eigs * ha2ev
        
        # Get unique values up to threshold (following original algorithm exactly)
        uni_eigs, degen_eigs = np.unique((BS_eigs_eV / degen_thres).astype(int),
                                        return_counts=True)
        uni_eigs = uni_eigs * degen_thres
        
        print('=' * 40)
        print('Group theory analysis for Q point : ', self.kpts_iBZ[iQ - 1])
        print('*' * 40)

        # Find little group (following original algorithm exactly)
        trace_all_real = []
        trace_all_imag = []
        little_group = []
        # Loop over symmetries (excluding time reversal operations)
        for isym in range(int(self.sym_red.shape[0] / (self.ele_time_rev + 1))):
            # Check if Sq = q (following original algorithm exactly)
            Sq_minus_q = np.einsum('ij,j->i', self.sym_red[isym],
                                  self.kpts_iBZ[iQ - 1]) - self.kpts_iBZ[iQ - 1]
            Sq_minus_q = Sq_minus_q - np.rint(Sq_minus_q)
            
            # Check if Sq = q (within tolerance)
            if np.linalg.norm(Sq_minus_q) > 1e-5:
                continue
            little_group.append(isym + 1)
            # Phase factor from fractional translations
            tau_dot_k = np.exp(1j * 2 * np.pi *
                              np.dot(self.kpts_iBZ[iQ - 1], self.frac_trans[isym]))
            
            # Rotate exciton wavefunction
            wfc_tmp = rotate_exc_wf(
                BS_wfcs,
                self.sym_red[isym],
                self.kpts,
                self.kpts_iBZ[iQ - 1],
                self.Dmats[isym],
                False,
                ktree=self.kpt_tree
            )
            
            # Compute representation matrix (following original algorithm exactly)
            rep = np.einsum('n...,m...->nm', wfc_tmp, BS_wfcs.conj(),
                           optimize=True) * tau_dot_k
            
            # Compute traces for each degenerate subspace (following original algorithm)
            irrep_sum = 0
            real_trace = []
            imag_trace = []
            for iirepp in range(len(uni_eigs)):
                idegen = degen_eigs[iirepp]
                idegen2 = irrep_sum + idegen
                trace_tmp = np.trace(rep[irrep_sum:idegen2, irrep_sum:idegen2])
                real_trace.append(trace_tmp.real.round(4))
                imag_trace.append(trace_tmp.imag.round(4))
                irrep_sum = idegen2
                
            trace_all_real.append(real_trace)
            trace_all_imag.append(imag_trace)

        little_group = np.array(little_group, dtype=int)
        
        # Get point group information (following original algorithm)
        try:
            from .point_group_ops import get_pg_info, decompose_rep2irrep
            little_group_mats = self.symm_mats[little_group - 1]
            # Debug info removed for cleaner output
            
            # Use original matrices without rounding (like the reference script)
            pg_label, classes, class_dict, char_tab, irreps = get_pg_info(little_group_mats)
        except ImportError:
            print("Warning: Point group analysis module not available")
            pg_label = "Unknown"
            classes = []
            class_dict = {}
            char_tab = None
            irreps = []
        except Exception as e:
            print(f"Warning: Point group analysis failed due to numerical precision issues.")
            print("Continuing with basic symmetry analysis...")
            pg_label = "Unknown"
            classes = []
            class_dict = {}
            char_tab = None
            irreps = []

        print('Little group : ', pg_label)
        print('Little group symmetries : ', little_group)

        # Print class information (following original algorithm exactly)
        irrep_decompositions = []
        if classes:
            print('Classes (symmetry indices in each class): ')
            req_sym_characters = np.zeros(len(classes), dtype=int)
            class_orders = np.zeros(len(classes), dtype=int)
            for ilab, iclass in class_dict.items():
                if ilab < len(classes):  # Safety check
                    print("%16s    : " % (classes[ilab]), little_group[np.array(iclass)])
                    req_sym_characters[ilab] = min(iclass)
                    class_orders[ilab] = len(iclass)
                else:
                    print(f"Warning: Class index {ilab} out of range for classes list (len={len(classes)})")
            print()

            # Process traces (following original algorithm exactly)
            trace_all_real = np.array(trace_all_real)
            trace_all_imag = np.array(trace_all_imag)
            trace = trace_all_real + 1j * trace_all_imag
            trace_req = trace[req_sym_characters, :].T

            print("====== Exciton representations ======")
            print("Energy (eV),  degeneracy  : representation")
            print('-' * 40)
            
            # Decompose representations (following original algorithm exactly)
            for i in range(len(trace_req)):
                if char_tab is not None:
                    rep_str_tmp = decompose_rep2irrep(trace_req[i], char_tab, 
                                                     len(little_group),
                                                     class_orders, irreps)
                else:
                    rep_str_tmp = "Analysis not available"
                print('%.4f        %9d  : ' % (uni_eigs[i], degen_eigs[i]), rep_str_tmp)
                irrep_decompositions.append(rep_str_tmp)
        else:
            # Fallback when point group analysis fails
            print("====== Exciton representations ======")
            print("Energy (eV),  degeneracy  : representation")
            print('-' * 40)
            for i in range(len(uni_eigs)):
                rep_str_tmp = "Point group analysis failed"
                print('%.4f        %9d  : ' % (uni_eigs[i], degen_eigs[i]), rep_str_tmp)
                irrep_decompositions.append(rep_str_tmp)

        print('*' * 40)

        # Return results
        results = {
            'q_point': self.kpts_iBZ[iQ - 1],
            'little_group': little_group,
            'point_group_label': pg_label,
            'unique_energies': uni_eigs,
            'degeneracies': degen_eigs,
            'irrep_decomposition': irrep_decompositions,
            'exciton_energies': BS_eigs_eV,
            'classes': classes,
            'class_dict': class_dict,
            'trace_characters': np.array(trace_all_real) + 1j * np.array(trace_all_imag) if len(trace_all_real) > 0 else None
        }
        
        return results

    def save_analysis_results(self, results, filename=None):
        """
        Save the group theory analysis results to a file.

        Parameters
        ----------
        results : dict
            Results dictionary from analyze_exciton_symmetry.
        filename : str, optional
            Output filename. If None, uses default naming.
        """
        if filename is None:
            q_str = '_'.join([f'{q:.3f}' for q in results['q_point']])
            filename = f'exciton_group_theory_Q{q_str}.txt'
        
        with open(filename, 'w') as f:
            f.write("Exciton Group Theory Analysis\n")
            f.write("=" * 40 + "\n")
            f.write(f"Q-point: {results['q_point']}\n")
            f.write(f"Little group: {results['point_group_label']}\n")
            f.write(f"Little group symmetries: {results['little_group']}\n\n")
            
            if results['classes']:
                f.write("Classes:\n")
                for class_name in results['classes']:
                    f.write(f"  {class_name}\n")
                f.write("\n")
            
            f.write("Exciton representations:\n")
            f.write("Energy (eV)    Degeneracy    Representation\n")
            f.write("-" * 50 + "\n")
            
            for i, (energy, degen, irrep) in enumerate(zip(
                results['unique_energies'], 
                results['degeneracies'],
                results['irrep_decomposition'])):
                f.write(f"{energy:8.4f}    {degen:8d}    {irrep}\n")
        
        print(f"Analysis results saved to {filename}")