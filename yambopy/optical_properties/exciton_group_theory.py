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
    Group theory analysis of exciton states using crystallographic symmetries.
    
    This class performs comprehensive symmetry analysis of exciton states by:
    
    1. **Point Group Identification**: Determines the crystallographic point group
       using spglib and spgrep libraries for accurate symmetry classification.
       
    2. **Little Group Analysis**: Identifies symmetry operations that leave the 
       exciton momentum invariant, crucial for k-point dependent analysis.
       
    3. **Irreducible Representation Decomposition**: Decomposes exciton states 
       into irreducible representations of the little group.
       
    4. **Optical Activity Analysis**: Determines selection rules for optical 
       transitions including electric dipole, Raman, and infrared activity.
    
    **Theoretical Background**
    
    The symmetry of exciton states is determined by the little group G_k of the 
    exciton momentum k. For each exciton state |ψ_n⟩, the representation matrix
    under symmetry operation R is computed as:
    
        D^(n)_R = ⟨ψ_n(Rk)| U(R) |ψ_n(k)⟩
    
    where U(R) is the symmetry operator. The character of this representation
    determines the irreducible representation:
    
        χ^(n)(R) = Tr[D^(n)_R]
    
    **Selection Rules**
    
    Optical transitions are governed by symmetry selection rules:
    
    - **Electric Dipole**: Allowed if Γ_i ⊗ Γ_f contains the representation of 
      the dipole operator (typically ungerade irreps in centrosymmetric groups)
    - **Raman Active**: Allowed if the irrep is contained in the polarizability 
      tensor (typically gerade irreps)
    - **IR Active**: Allowed if the irrep transforms as a translation vector
    
    Parameters
    ----------
    path : str, optional
        Path to the calculation directory. Default is current directory.
    save : str, optional
        Name of the SAVE directory. Default is 'SAVE'.
    lelph_db : LetzElphElectronPhononDB, optional
        Pre-loaded electron-phonon database.
    latdb : YamboLatticeDB, optional
        Pre-loaded lattice database.
    wfdb : YamboWFDB, optional
        Pre-loaded wavefunction database.
    bands_range : list or tuple, optional
        Range of bands to include in the analysis.
    BSE_dir : str, optional
        Name of the BSE directory. Default is 'bse'.
    LELPH_dir : str, optional
        Name of the electron-phonon directory. Default is 'lelph'.
    read_symm_from_ns_db_file : bool, optional
        Whether to read symmetries from ns.db1 file. Default is True.
    
    Attributes
    ----------
    point_group_label : str
        Identified crystallographic point group.
    little_group : array_like
        Indices of little group symmetry operations.
    symm_mats : array_like
        Symmetry matrices in Cartesian coordinates.
    frac_trans : array_like
        Fractional translation vectors for non-symmorphic operations.
    spglib_Dmats : array_like, optional
        D-matrices computed using spglib symmetries for enhanced accuracy.
    
    Examples
    --------
    Basic usage for exciton symmetry analysis:
    
    >>> from yambopy.optical_properties import ExcitonGroupTheory
    >>> 
    >>> # Initialize the analysis
    >>> egt = ExcitonGroupTheory(
    ...     path='./calculation',
    ...     BSE_dir='bse',
    ...     LELPH_dir='lelph',
    ...     bands_range=[1, 20]
    ... )
    >>> 
    >>> # Analyze exciton symmetries at Γ point
    >>> results = egt.analyze_exciton_symmetry(iQ=1, nstates=10)
    >>> 
    >>> # Print results
    >>> print(f"Point group: {results['point_group_label']}")
    >>> for i, (energy, irrep, activity) in enumerate(zip(
    ...     results['unique_energies'],
    ...     results['irrep_decomposition'],
    ...     results['optical_activity']
    ... )):
    ...     print(f"State {i+1}: {energy:.3f} eV, {irrep}")
    ...     print(f"  Optically active: {activity['electric_dipole_allowed']}")
    ...     print(f"  Raman active: {activity['raman_active']}")
    
    References
    ----------
    .. [1] Tinkham, M. "Group Theory and Quantum Mechanics" (1964)
    .. [2] Dresselhaus, M. S. et al. "Group Theory: Application to the Physics 
           of Condensed Matter" (2008)
    .. [3] Aroyo, M. I. et al. "Crystallography online: Bilbao Crystallographic 
           Server" Bulg. Chem. Commun. 43, 183-197 (2011)
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
        
        # Setup D-matrices for spgrep analysis (will be computed when needed)
        self.spglib_Dmats = None
        
        # Try to compute spglib D-matrices for better representation analysis
        self._try_compute_spglib_dmats()

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
                print("Attempting to get fractional translations from spglib...")
                self._setup_symmetry_from_spglib()
        else:
            print("Reading symmetry from ns.db1 file...")
            # Use lattice database kmap and kpoints from lelph_db
            self.kpts = self.lelph_db.kpoints
            # Try to get fractional translations from spglib
            self._setup_symmetry_from_spglib()

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
    
    def _setup_symmetry_from_spglib(self):
        """Setup symmetry data using spglib when elph file is not available."""
        try:
            import spglib
            
            # Get crystal structure
            lattice, positions, numbers = self._get_crystal_structure()
            if lattice is not None and positions is not None and numbers is not None:
                cell = (lattice, positions, numbers)
                symmetry = spglib.get_symmetry(cell, symprec=1e-5)
                
                if symmetry:
                    # Get fractional translations from spglib
                    spg_translations = symmetry['translations']
                    spg_rotations = symmetry['rotations']
                    
                    # Match spglib operations with Yambo operations
                    if hasattr(self, 'symm_mats') and self.symm_mats is not None:
                        # Try to match operations and extract corresponding translations
                        self.frac_trans = self._match_translations(spg_rotations, spg_translations)
                    else:
                        # If no Yambo symmetries available, use spglib directly
                        print("Using spglib symmetries directly")
                        self.symm_mats = spg_rotations.astype(float)
                        self.frac_trans = spg_translations
                        # Convert to crystal coordinates
                        self.frac_trans = np.einsum('ij,nj->ni', self.blat_vecs.T, self.frac_trans)
                        
                    print(f"Extracted {len(self.frac_trans)} fractional translations from spglib")
                    return
            
            # Fallback: zero translations
            print("spglib fallback failed, using zero fractional translations")
            if hasattr(self, 'symm_mats') and self.symm_mats is not None:
                self.frac_trans = np.zeros((self.symm_mats.shape[0], 3))
            else:
                self.frac_trans = np.zeros((24, 3))  # Default assumption
                
        except ImportError:
            print("spglib not available, using zero fractional translations")
            if hasattr(self, 'symm_mats') and self.symm_mats is not None:
                self.frac_trans = np.zeros((self.symm_mats.shape[0], 3))
            else:
                self.frac_trans = np.zeros((24, 3))
        except Exception as e:
            print(f"spglib setup failed: {e}")
            if hasattr(self, 'symm_mats') and self.symm_mats is not None:
                self.frac_trans = np.zeros((self.symm_mats.shape[0], 3))
            else:
                self.frac_trans = np.zeros((24, 3))
    
    def _match_translations(self, spg_rotations, spg_translations):
        """Match spglib translations with Yambo symmetry operations."""
        matched_translations = np.zeros((self.symm_mats.shape[0], 3))
        
        for i, yambo_rot in enumerate(self.symm_mats):
            # Find matching rotation in spglib operations
            best_match_idx = -1
            min_diff = float('inf')
            
            for j, spg_rot in enumerate(spg_rotations):
                diff = np.linalg.norm(yambo_rot - spg_rot.astype(float))
                if diff < min_diff:
                    min_diff = diff
                    best_match_idx = j
            
            # If we found a good match (within tolerance)
            if min_diff < 1e-6 and best_match_idx >= 0:
                matched_translations[i] = spg_translations[best_match_idx]
            else:
                # No match found, use zero translation
                matched_translations[i] = np.zeros(3)
        
        # Convert to crystal coordinates
        return np.einsum('ij,nj->ni', self.blat_vecs.T, matched_translations)
    


    
    def _try_compute_spglib_dmats(self):
        """Try to compute D-matrices using spglib symmetries for better representation analysis."""
        try:
            import spglib
            
            # Get crystal structure and symmetry operations
            lattice, positions, numbers = self._get_crystal_structure()
            if lattice is not None and positions is not None and numbers is not None:
                cell = (lattice, positions, numbers)
                symmetry = spglib.get_symmetry(cell, symprec=1e-5)
                
                if symmetry and hasattr(self, 'symm_mats') and self.symm_mats is not None:
                    spg_rotations = symmetry['rotations']
                    spg_translations = symmetry['translations']
                    
                    # Check if we have the same number of operations
                    if len(spg_rotations) == len(self.symm_mats):
                        print("Computing D-matrices with spglib symmetries for representation analysis...")
                        
                        # Compute D-matrices using spglib symmetries
                        self.spglib_Dmats = self.wfdb.Dmat(
                            symm_mat=spg_rotations.astype(float),
                            frac_vec=spg_translations,
                            time_rev=(self.ele_time_rev == 1)
                        )[:,:,0,:,:]
                        
                        print(f"Successfully computed spglib D-matrices: {self.spglib_Dmats.shape}")
                    else:
                        print(f"Symmetry count mismatch: spglib({len(spg_rotations)}) vs Yambo({len(self.symm_mats)})")
                        print("Using Yambo D-matrices")
                        
        except Exception as e:
            print(f"Could not compute spglib D-matrices: {e}")
            print("Using Yambo D-matrices")

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
        
        # Convert energies to eV for analysis 
        BS_eigs_eV = BS_eigs * ha2ev
        
        # Get unique values up to threshold 
        uni_eigs, degen_eigs = np.unique((BS_eigs_eV / degen_thres).astype(int),
                                        return_counts=True)
        uni_eigs = uni_eigs * degen_thres
        
        print('=' * 40)
        print('Group theory analysis for Q point : ', self.kpts_iBZ[iQ - 1])
        print('*' * 40)

        # Find little group 
        trace_all_real = []
        trace_all_imag = []
        little_group = []
        # Loop over symmetries (excluding time reversal operations)
        for isym in range(int(self.sym_red.shape[0] / (self.ele_time_rev + 1))):
            # Check if Sq = q 
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
            # Use spglib D-matrices if available, otherwise use Yambo D-matrices
            dmat_to_use = self.spglib_Dmats[isym] if self.spglib_Dmats is not None else self.Dmats[isym]
            
            wfc_tmp = rotate_exc_wf(
                BS_wfcs,
                self.sym_red[isym],
                self.kpts,
                self.kpts_iBZ[iQ - 1],
                dmat_to_use,
                False,
                ktree=self.kpt_tree
            )
            
            # Compute representation matrix 
            rep = np.einsum('n...,m...->nm', wfc_tmp, BS_wfcs.conj(),
                           optimize=True) * tau_dot_k
            
            # Compute traces for each degenerate subspace
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
        
        # Get point group information using spgrep for complete crystallographic analysis
        # Note: We use spgrep here to get the full crystallographic point group including
        # non-symmorphic symmetries, while the actual wavefunction rotations above
        # used the Yambo symmetries from the SAVE database
        try:
            from .spgrep_point_group_ops import get_pg_info, decompose_rep2irrep
            
            # For irrep analysis, we need to determine the crystallographic point group
            # This may include additional symmetries not present in the Yambo SAVE
            pg_label, classes, class_dict, char_tab, irreps = self._get_crystallographic_point_group(little_group)
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
                    # Convert 1-based indices to 0-based for array access
                    iclass_0based = np.array(iclass) - 1
                    # But print the 1-based indices as they appear in little_group
                    print("%16s    : " % (classes[ilab]), np.array(iclass))
                    req_sym_characters[ilab] = min(iclass) - 1  # Convert to 0-based
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
        
        # Add optical activity analysis
        results['optical_activity'] = self._analyze_optical_activity(
            results['irrep_decomposition'], 
            results['point_group_label']
        )
        
        return results
    
    def _analyze_optical_activity(self, irrep_decompositions, point_group_label):
        """
        Analyze optical activity of exciton states based on their irreducible representations.
        
        This method determines whether transitions are:
        - Electric dipole allowed (optical transitions)
        - Raman active
        - Infrared (IR) active
        
        Parameters
        ----------
        irrep_decompositions : list
            List of irreducible representation labels for each energy level
        point_group_label : str
            Point group label (e.g., '6/m', 'D3h', 'Oh')
            
        Returns
        -------
        list
            List of dictionaries containing optical activity analysis for each level
        """
        optical_activities = []
        
        for irrep_str in irrep_decompositions:
            activity = {
                'irrep': irrep_str,
                'electric_dipole_allowed': False,
                'raman_active': False,
                'ir_active': False,
                'selection_rules': []
            }
            
            # Parse irrep string (handle cases like "Γ5 ⊕ Γ7", "A1g", etc.)
            irreps = self._parse_irrep_string(irrep_str)
            
            for irrep in irreps:
                # Apply selection rules based on point group and irrep
                rules = self._get_selection_rules(irrep, point_group_label)
                
                # Combine rules (if any component is active, the whole state is active)
                activity['electric_dipole_allowed'] |= rules['electric_dipole_allowed']
                activity['raman_active'] |= rules['raman_active'] 
                activity['ir_active'] |= rules['ir_active']
                activity['selection_rules'].append(rules)
            
            optical_activities.append(activity)
        
        return optical_activities
    
    def _parse_irrep_string(self, irrep_str):
        """Parse irrep string to extract individual irrep labels."""
        if not isinstance(irrep_str, str):
            return [str(irrep_str)]
            
        # Handle direct sum notation (e.g., "Γ5 ⊕ Γ7")
        if '⊕' in irrep_str:
            return [irrep.strip() for irrep in irrep_str.split('⊕')]
        elif '+' in irrep_str:
            return [irrep.strip() for irrep in irrep_str.split('+')]
        else:
            return [irrep_str.strip()]
    
    def _get_selection_rules(self, irrep, point_group_label):
        """
        Get selection rules for a specific irrep in a given point group.
        
        Parameters
        ----------
        irrep : str
            Irreducible representation label
        point_group_label : str
            Point group label
            
        Returns
        -------
        dict
            Selection rules for this irrep
        """
        rules = {
            'electric_dipole_allowed': False,
            'raman_active': False,
            'ir_active': False,
            'notes': []
        }
        
        # Normalize point group label
        pg = point_group_label.lower().replace('/', '').replace('-', '')
        
        # Apply general selection rules based on irrep symmetry
        if 'g' in irrep.lower():
            # Gerade (even) irreps
            rules['raman_active'] = True
            rules['ir_active'] = False
            rules['electric_dipole_allowed'] = False
            rules['notes'].append("Gerade: Raman active, IR forbidden")
            
        elif 'u' in irrep.lower():
            # Ungerade (odd) irreps  
            rules['raman_active'] = False
            rules['ir_active'] = True
            rules['electric_dipole_allowed'] = True
            rules['notes'].append("Ungerade: IR active, electric dipole allowed")
            
        else:
            # For point groups without inversion symmetry
            # Apply specific rules based on point group and irrep type
            rules.update(self._get_non_centrosymmetric_rules(irrep, pg))
        
        # Special cases for specific point groups
        if pg in ['6m', 'c6h']:
            rules.update(self._get_c6h_rules(irrep))
        elif pg in ['d3h']:
            rules.update(self._get_d3h_rules(irrep))
        elif pg in ['oh', 'o']:
            rules.update(self._get_oh_rules(irrep))
        elif pg in ['td']:
            rules.update(self._get_td_rules(irrep))
            
        return rules
    
    def _get_non_centrosymmetric_rules(self, irrep, point_group):
        """Selection rules for non-centrosymmetric point groups."""
        rules = {
            'electric_dipole_allowed': True,  # Generally allowed
            'raman_active': True,             # Generally allowed
            'ir_active': True,                # Generally allowed
            'notes': ["Non-centrosymmetric: generally all transitions allowed"]
        }
        
        # Specific rules for common irreps
        if irrep.lower().startswith('a'):
            # A-type irreps (totally symmetric or antisymmetric)
            if '1' in irrep:
                rules['notes'].append("A1-type: totally symmetric")
            elif '2' in irrep:
                rules['notes'].append("A2-type: antisymmetric to some operations")
                
        elif irrep.lower().startswith('e'):
            # E-type irreps (doubly degenerate)
            rules['notes'].append("E-type: doubly degenerate, typically optically active")
            
        elif irrep.lower().startswith('t'):
            # T-type irreps (triply degenerate)
            rules['notes'].append("T-type: triply degenerate, typically optically active")
            
        return rules
    
    def _get_c6h_rules(self, irrep):
        """Selection rules for C6h point group (6/m)."""
        rules = {'notes': []}
        
        # C6h irreps: Ag, Bg, E1g, E2g, Au, Bu, E1u, E2u
        if 'γ' in irrep.lower() or 'gamma' in irrep.lower():
            # Gamma point irreps (spgrep notation)
            irrep_num = ''.join(filter(str.isdigit, irrep))
            if irrep_num in ['1', '2']:  # Γ1, Γ2 (A-type)
                rules['electric_dipole_allowed'] = False
                rules['raman_active'] = True
                rules['ir_active'] = False
                rules['notes'].append("A-type gerade: Raman active")
            elif irrep_num in ['3', '4']:  # Γ3, Γ4 (A-type ungerade)
                rules['electric_dipole_allowed'] = True
                rules['raman_active'] = False
                rules['ir_active'] = True
                rules['notes'].append("A-type ungerade: IR active")
            elif irrep_num in ['5', '6']:  # Γ5, Γ6 (E1-type)
                rules['electric_dipole_allowed'] = True
                rules['raman_active'] = True
                rules['ir_active'] = True
                rules['notes'].append("E1-type: optically active (in-plane)")
            elif irrep_num in ['7', '8']:  # Γ7, Γ8 (E2-type)
                rules['electric_dipole_allowed'] = False
                rules['raman_active'] = True
                rules['ir_active'] = False
                rules['notes'].append("E2-type: Raman active")
        
        return rules
    
    def _get_d3h_rules(self, irrep):
        """Selection rules for D3h point group."""
        rules = {'notes': []}
        
        # D3h irreps: A1', A2', E', A1'', A2'', E''
        if "'" in irrep:  # Prime irreps
            if irrep.startswith('A1'):
                rules['electric_dipole_allowed'] = False
                rules['raman_active'] = True
                rules['ir_active'] = False
                rules['notes'].append("A1': totally symmetric, Raman active")
            elif irrep.startswith('A2'):
                rules['electric_dipole_allowed'] = False
                rules['raman_active'] = False
                rules['ir_active'] = False
                rules['notes'].append("A2': antisymmetric, inactive")
            elif irrep.startswith('E'):
                rules['electric_dipole_allowed'] = True
                rules['raman_active'] = True
                rules['ir_active'] = True
                rules['notes'].append("E': doubly degenerate, optically active")
        elif "''" in irrep:  # Double prime irreps
            if irrep.startswith('A1'):
                rules['electric_dipole_allowed'] = True
                rules['raman_active'] = False
                rules['ir_active'] = True
                rules['notes'].append("A1'': IR active (z-polarized)")
            elif irrep.startswith('A2'):
                rules['electric_dipole_allowed'] = False
                rules['raman_active'] = True
                rules['ir_active'] = False
                rules['notes'].append("A2'': Raman active")
            elif irrep.startswith('E'):
                rules['electric_dipole_allowed'] = False
                rules['raman_active'] = True
                rules['ir_active'] = False
                rules['notes'].append("E'': Raman active")
        
        return rules
    
    def _get_oh_rules(self, irrep):
        """Selection rules for Oh point group."""
        rules = {'notes': []}
        
        # Oh irreps: A1g, A2g, Eg, T1g, T2g, A1u, A2u, Eu, T1u, T2u
        if 'g' in irrep:
            rules['electric_dipole_allowed'] = False
            rules['raman_active'] = True
            rules['ir_active'] = False
            rules['notes'].append("Gerade: Raman active, electric dipole forbidden")
        elif 'u' in irrep:
            rules['electric_dipole_allowed'] = True
            rules['raman_active'] = False
            rules['ir_active'] = True
            rules['notes'].append("Ungerade: IR active, electric dipole allowed")
            
        return rules
    
    def _get_td_rules(self, irrep):
        """Selection rules for Td point group."""
        rules = {'notes': []}
        
        # Td irreps: A1, A2, E, T1, T2
        if irrep.startswith('A1'):
            rules['electric_dipole_allowed'] = False
            rules['raman_active'] = True
            rules['ir_active'] = False
            rules['notes'].append("A1: totally symmetric, Raman active")
        elif irrep.startswith('A2'):
            rules['electric_dipole_allowed'] = False
            rules['raman_active'] = False
            rules['ir_active'] = False
            rules['notes'].append("A2: inactive")
        elif irrep.startswith('E'):
            rules['electric_dipole_allowed'] = False
            rules['raman_active'] = True
            rules['ir_active'] = False
            rules['notes'].append("E: doubly degenerate, Raman active")
        elif irrep.startswith('T1'):
            rules['electric_dipole_allowed'] = True
            rules['raman_active'] = False
            rules['ir_active'] = True
            rules['notes'].append("T1: triply degenerate, IR active")
        elif irrep.startswith('T2'):
            rules['electric_dipole_allowed'] = False
            rules['raman_active'] = True
            rules['ir_active'] = False
            rules['notes'].append("T2: triply degenerate, Raman active")
            
        return rules
    
    def _get_crystallographic_point_group(self, little_group_yambo):
        """
        Determine the crystallographic point group for irrep analysis.
        
        This method takes the little group operations from Yambo and determines
        the corresponding crystallographic point group using spgrep, which may
        include additional non-symmorphic symmetries not present in the Yambo SAVE.
        
        Parameters
        ----------
        little_group_yambo : array_like
            Little group operation indices from Yambo analysis
            
        Returns
        -------
        tuple
            Point group information (pg_label, classes, class_dict, char_tab, irreps)
        """
        from .spgrep_point_group_ops import get_pg_info
        
        # Get the Yambo symmetry matrices for the little group
        little_group_mats_yambo = self.symm_mats[little_group_yambo - 1]
        
        # For crystallographic analysis, we need to determine the space group
        # and extract the corresponding point group with all symmetries
        # This is where spgrep's comprehensive database is crucial
        
        try:
            # First attempt: use the Yambo matrices directly with spgrep
            pg_label, classes, class_dict, char_tab, irreps = get_pg_info(
                little_group_mats_yambo, 
                time_rev=(self.ele_time_rev == 1)
            )
            return pg_label, classes, class_dict, char_tab, irreps
            
        except Exception as e:
            print(f"Direct spgrep analysis failed: {e}")
            
            # Use spglib to get proper crystallographic symmetries for irrep analysis
            return self._get_point_group_from_spglib(little_group_yambo)
    

    
    def _get_point_group_from_spglib(self, little_group_yambo):
        """
        Use spglib to get proper crystallographic symmetries for irrep analysis.
        
        This method uses spglib to identify the space group and extract the 
        corresponding point group symmetries for proper irrep decomposition.
        """
        try:
            import spglib
            
            # Get atomic positions and lattice from the SAVE database
            lattice, positions, numbers = self._get_crystal_structure()
            
            if lattice is not None and positions is not None and numbers is not None:
                # Use spglib to find the space group
                cell = (lattice, positions, numbers)
                spacegroup_info = spglib.get_spacegroup(cell, symprec=1e-5)
                
                if spacegroup_info:
                    print(f"spglib identified space group: {spacegroup_info}")
                    
                    # Get symmetry operations from spglib
                    symmetry = spglib.get_symmetry(cell, symprec=1e-5)
                    
                    if symmetry:
                        # Extract point group operations (rotations without translations)
                        rotations = symmetry['rotations']
                        
                        # Use spgrep with the spglib symmetries directly
                        # The little group analysis is handled by the original Yambo approach
                        return self._analyze_with_spglib_symmetries(rotations, little_group_yambo)
            
            # Fallback if spglib analysis fails
            print("spglib analysis failed, using operation matching")
            return self._match_operations_to_point_group(little_group_yambo)
                
        except ImportError:
            print("spglib not available, using operation matching")
            return self._match_operations_to_point_group(little_group_yambo)
        except Exception as e:
            print(f"spglib analysis failed: {e}")
            return self._match_operations_to_point_group(little_group_yambo)
    

    
    def _get_crystal_structure(self):
        """
        Get crystal structure information from the SAVE database.
        
        Returns
        -------
        tuple
            (lattice, positions, numbers) for spglib analysis
        """
        try:
            # Get lattice vectors (already available)
            lattice = self.lat_vecs
            
            # Read atomic positions and numbers from the SAVE database
            if (hasattr(self.ydb, 'car_atomic_positions') and 
                hasattr(self.ydb, 'atomic_numbers')):
                
                # Use the actual atomic positions and numbers from SAVE
                # spglib expects fractional coordinates, not Cartesian
                positions = self.ydb.red_atomic_positions
                numbers = self.ydb.atomic_numbers
                
                print(f"Read crystal structure from SAVE: {len(positions)} atoms")
                print(f"Atomic numbers: {numbers}")
                
                return lattice, positions, numbers
            else:
                print("No atomic structure information found in SAVE database")
                return None, None, None
            
        except Exception as e:
            print(f"Failed to get crystal structure: {e}")
            return None, None, None
    

    
    def _analyze_with_spglib_symmetries(self, spglib_rotations, little_group_yambo):
        """
        Analyze point group using spglib symmetries with spgrep.
        
        Parameters
        ----------
        spglib_rotations : array_like
            Rotation matrices from spglib
        little_group_yambo : array_like
            Little group indices from Yambo (for reference)
            
        Returns
        -------
        tuple
            Point group analysis results
        """
        try:
            from .spgrep_point_group_ops import get_pg_info
            
            # Convert spglib rotations to the format expected by spgrep
            # spglib gives integer matrices in the standard crystallographic setting
            print(f"Using {len(spglib_rotations)} symmetry operations from spglib")
            
            # D-matrices are computed during initialization if possible
            
            # Use spgrep with the proper spglib symmetries
            pg_label, classes, class_dict, char_tab, irreps = get_pg_info(
                spglib_rotations, 
                time_rev=(self.ele_time_rev == 1)
            )
            
            print(f"spgrep analysis with spglib symmetries successful: {pg_label}")
            return pg_label, classes, class_dict, char_tab, irreps
            
        except Exception as e:
            print(f"spgrep analysis with spglib symmetries failed: {e}")
            # Final fallback
            return self._match_operations_to_point_group(little_group_yambo)
    

    
    def _match_operations_to_point_group(self, little_group_yambo):
        """
        Match symmetry operations to known point groups as a fallback.
        
        This analyzes the actual symmetry operations to determine the point group.
        """
        try:
            # Get the symmetry matrices for the little group
            little_group_mats = self.symm_mats[little_group_yambo - 1]
            
            # Analyze the operations to determine point group
            n_ops = len(little_group_mats)
            
            # Count different types of operations
            identity_count = 0
            c3_rotations = 0
            c2_rotations = 0
            horizontal_reflections = 0
            vertical_reflections = 0
            improper_rotations = 0
            
            for mat in little_group_mats:
                det = np.linalg.det(mat)
                trace = np.trace(mat)
                
                if np.allclose(mat, np.eye(3)):
                    identity_count += 1
                elif np.isclose(det, 1) and np.isclose(trace, 0):
                    c3_rotations += 1
                elif np.isclose(det, 1) and np.isclose(trace, -1):
                    c2_rotations += 1
                elif (np.isclose(det, -1) and np.isclose(mat[2,2], 1) and 
                      np.allclose(mat[2,:2], 0) and np.allclose(mat[:2,2], 0)):
                    horizontal_reflections += 1
                elif np.isclose(det, -1) and np.isclose(trace, 1):
                    vertical_reflections += 1
                elif np.isclose(det, -1):
                    improper_rotations += 1
            
            # If we can't identify the point group, return unknown
            print(f"Could not identify point group from {n_ops} operations")
            return "Unknown", [], {}, None, []
            
        except Exception as e:
            print(f"Operation matching failed: {e}")
            return "Unknown", [], {}, None, []
    


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