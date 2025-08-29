#
# License-Identifier: GPL
#
# Copyright (C) 2024 The Yambo Team
#
# Authors: RR, MN
#
# This file is part of the yambopy project
#
"""
Base class for optical properties calculations to reduce code duplication.
"""

import warnings
import os
import numpy as np
from abc import ABC, abstractmethod
from tqdm import tqdm

from yambopy import YamboLatticeDB
from yambopy.dbs.wfdb import YamboWFDB
from yambopy.dbs.excitondb import YamboExcitonDB
from yambopy.dbs.dipolesdb import YamboDipolesDB
from yambopy.units import ha2ev
from yambopy.kpoints import build_ktree, find_kpt
from yambopy.geometry_manager import get_geometry_manager

try:
    from pykdtree.kdtree import KDTree 
    # pykdtree is much faster and is recommended
    # pip install pykdtree
    # useful in Dmat computation
except ImportError:
    from scipy.spatial import KDTree

warnings.filterwarnings('ignore')


class BaseOpticalProperties(ABC):
    """
    Base class for optical properties calculations.
    
    This class provides common functionality for reading Yambo databases,
    setting up k-point trees, and handling common data structures.
    """
    
    def __init__(self, path=None, save='SAVE', latdb=None, wfdb=None, 
                 bands_range=None, BSE_dir='bse', save_files=True, use_geometry_manager=True):
        """
        Initialize base optical properties class.
        
        Parameters
        ----------
        path : str, optional
            Path to the directory containing the calculation files. 
            Defaults to current directory.
        save : str, optional
            Subdirectory containing the SAVE files. Defaults to 'SAVE'.
        latdb : YamboLatticeDB, optional
            Pre-loaded Yambo Lattice database. Defaults to None.
        wfdb : YamboWFDB, optional
            Pre-loaded Yambo wavefunction database. Defaults to None.
        bands_range : list, optional
            Range of bands to load. Python indexing. Defaults to all bands.
        BSE_dir : str, optional
            Subdirectory containing BSE output. Defaults to 'bse'.
        save_files : bool, optional
            Whether to save files in .npy database. Defaults to True.
        use_geometry_manager : bool, optional
            Whether to use the geometry manager for shared geometry data. Defaults to True.
        """
        self.path = path if path is not None else os.getcwd()
        self.SAVE_dir = os.path.join(self.path, save)
        self.BSE_dir = os.path.join(self.path, BSE_dir)
        self.save_files = save_files
        self.bands_range = bands_range if bands_range is not None else []
        self.use_geometry_manager = use_geometry_manager
        
        # Initialize databases
        self.latdb = latdb
        self.wfdb = wfdb
        
        # Initialize geometry manager if requested
        if self.use_geometry_manager:
            self.geometry_manager = get_geometry_manager(self.path, save, latdb)
            # Initialize common attributes as None - they'll be accessed via properties
            self.ydb = None
            self.kpt_tree = None
            self.qidx_in_kpts = None
        else:
            self.geometry_manager = None
            # Initialize common attributes the traditional way
            self.ydb = None
            self.kpt_tree = None
            self.qidx_in_kpts = None
        
    def _setup_directories(self, **additional_dirs):
        """
        Setup additional directories based on provided keyword arguments.
        
        Parameters
        ----------
        **additional_dirs : dict
            Additional directories to setup (e.g., DIP_dir='gw', LELPH_dir='lelph')
        """
        for dir_name, dir_value in additional_dirs.items():
            setattr(self, dir_name, os.path.join(self.path, dir_value))
    
    # Properties for geometry manager integration
    @property
    def lat_vecs(self):
        """Get lattice vectors - from geometry manager if enabled, otherwise from ydb."""
        if self.use_geometry_manager and self.geometry_manager:
            return self.geometry_manager.get_lattice_info()['lat_vecs']
        elif hasattr(self, '_lat_vecs'):
            return self._lat_vecs
        elif self.ydb:
            return self.ydb.lat
        return None
    
    @lat_vecs.setter
    def lat_vecs(self, value):
        """Set lattice vectors."""
        if not self.use_geometry_manager:
            self._lat_vecs = value
    
    @property
    def blat_vecs(self):
        """Get reciprocal lattice vectors - from geometry manager if enabled, otherwise from ydb."""
        if self.use_geometry_manager and self.geometry_manager:
            return self.geometry_manager.get_lattice_info()['blat_vecs']
        elif hasattr(self, '_blat_vecs'):
            return self._blat_vecs
        elif self.ydb:
            return self.ydb.rlat.T
        return None
    
    @blat_vecs.setter
    def blat_vecs(self, value):
        """Set reciprocal lattice vectors."""
        if not self.use_geometry_manager:
            self._blat_vecs = value
    
    @property
    def nibz(self):
        """Get number of IBZ k-points - from geometry manager if enabled, otherwise from ydb."""
        if self.use_geometry_manager and self.geometry_manager:
            return self.geometry_manager.get_lattice_info()['nibz']
        elif hasattr(self, '_nibz'):
            return self._nibz
        elif self.ydb:
            return self.ydb.ibz_nkpoints
        return None
    
    @nibz.setter
    def nibz(self, value):
        """Set number of IBZ k-points."""
        if not self.use_geometry_manager:
            self._nibz = value
    
    @property
    def symm_mats(self):
        """Get symmetry matrices - from geometry manager if enabled, otherwise from ydb."""
        if self.use_geometry_manager and self.geometry_manager:
            return self.geometry_manager.get_symmetry_info()['symm_mats']
        elif hasattr(self, '_symm_mats'):
            return self._symm_mats
        elif self.ydb:
            return self.ydb.sym_car
        return None
    
    @symm_mats.setter
    def symm_mats(self, value):
        """Set symmetry matrices."""
        if not self.use_geometry_manager:
            self._symm_mats = value
    
    @property
    def ele_time_rev(self):
        """Get time reversal symmetry - from geometry manager if enabled, otherwise from ydb."""
        if self.use_geometry_manager and self.geometry_manager:
            return self.geometry_manager.get_symmetry_info()['ele_time_rev']
        elif hasattr(self, '_ele_time_rev'):
            return self._ele_time_rev
        elif self.ydb:
            return self.ydb.time_rev
        return None
    
    @ele_time_rev.setter
    def ele_time_rev(self, value):
        """Set time reversal symmetry."""
        if not self.use_geometry_manager:
            self._ele_time_rev = value
    
    @property
    def red_kpoints(self):
        """Get reduced k-points - from geometry manager if enabled, otherwise from ydb."""
        if self.use_geometry_manager and self.geometry_manager:
            return self.geometry_manager.get_kpoint_info()['red_kpoints']
        elif hasattr(self, '_red_kpoints'):
            return self._red_kpoints
        elif self.ydb:
            return self.ydb.red_kpoints
        return None
    
    @red_kpoints.setter
    def red_kpoints(self, value):
        """Set reduced k-points."""
        if not self.use_geometry_manager:
            self._red_kpoints = value
    
    def _read_lattice_db(self, latdb=None):
        """
        Read Yambo Lattice database.
        
        Parameters
        ----------
        latdb : YamboLatticeDB, optional
            Pre-loaded database. If None, reads from ns.db1 file.
        """
        if self.use_geometry_manager:
            # When using geometry manager, just ensure it's initialized
            if not self.geometry_manager.is_initialized():
                self.geometry_manager.initialize()
            # Set ydb for compatibility
            self.ydb = self.geometry_manager.get_lattice_info()['ydb']
        else:
            # Traditional approach
            try:
                ns_db1_fname = os.path.join(self.SAVE_dir, 'ns.db1')
                if latdb:
                    if not hasattr(latdb, 'ibz_kpoints'):
                        latdb.expand_kpoints()
                    self.ydb = latdb
                else:
                    self.ydb = YamboLatticeDB.from_db_file(ns_db1_fname, Expand=True)
            except Exception as e:
                raise IOError(f'Cannot read ns.db1 file: {e}')
            
            # Set common lattice properties using setters (which will store in private attributes)
            self.lat_vecs = self.ydb.lat
            self.nibz = self.ydb.ibz_nkpoints
            self.symm_mats = self.ydb.sym_car
            self.ele_time_rev = self.ydb.time_rev
            self.blat_vecs = self.ydb.rlat.T
    
    def _read_wavefunction_db(self, wfdb=None, bands_range=None):
        """
        Read Yambo wavefunction database.
        
        Parameters
        ----------
        wfdb : YamboWFDB, optional
            Pre-loaded database. If None, reads from ns.wf file.
        bands_range : list, optional
            Range of bands to load.
        """
        bands_range = bands_range or self.bands_range
        
        try:
            if wfdb:
                if not hasattr(wfdb, 'save_Dmat'):
                    wfdb.Dmat()
                self.wfdb = wfdb
            else:
                self.wfdb = YamboWFDB(path=self.path,save='SAVE', filename='ns.wf', latdb=self.ydb, 
                                    bands_range=bands_range)
        except Exception as e:
            raise IOError(f'Cannot read ns.wf file: {e}')
        
        # Set common wavefunction properties
        self.nkpoints = self.wfdb.nkpoints
        self.nspin = self.wfdb.nspin
        self.nspinor = self.wfdb.nspinor
        self.nbands = self.wfdb.nbands
    
    def _setup_kpoint_mapping(self):
        """Setup k-point mapping and symmetry operations."""
        if self.use_geometry_manager:
            # When using geometry manager, this is handled during initialization
            if not self.geometry_manager.is_initialized():
                self.geometry_manager.initialize()
            # Set kmap and sym_red for compatibility
            symmetry_info = self.geometry_manager.get_symmetry_info()
            self.kmap = symmetry_info.get('kmap')
            self.sym_red = symmetry_info.get('sym_red')
        else:
            # Traditional approach
            # Set kmap
            kmap = np.zeros((self.wfdb.nkBZ, 2), dtype=int)
            kmap[:, 0] = self.ydb.kpoints_indexes
            kmap[:, 1] = self.ydb.symmetry_indexes
            self.kmap = kmap
            
            # Setup symmetry operations in reduced coordinates
            temp = np.matmul(self.symm_mats, self.blat_vecs)
            sym_red = np.matmul(self.lat_vecs[None, :, :], temp)
            self.sym_red = np.rint(sym_red).astype(int)
    
    def _build_kpoint_tree(self, kpts=None):
        """
        Build k-point tree for efficient k-point searching.
        
        Parameters
        ----------
        kpts : array_like, optional
            K-points to build tree from. If None, uses red_kpoints.
        """
        if self.use_geometry_manager:
            # When using geometry manager, this is handled during initialization
            if not self.geometry_manager.is_initialized():
                self.geometry_manager.initialize()
            
            # If custom k-points are provided, rebuild the tree
            if kpts is not None:
                print('Building custom kD-tree for kpoints')
                self.kpts = kpts
                self.kpt_tree = build_ktree(kpts)
                self.qidx_in_kpts = find_kpt(self.kpt_tree, kpts)
            else:
                # Use the tree from geometry manager
                kpoint_info = self.geometry_manager.get_kpoint_info()
                self.kpts = kpoint_info['red_kpoints']
                self.kpt_tree = kpoint_info['kpt_tree']
                self.qidx_in_kpts = kpoint_info['qidx_in_kpts']
        else:
            # Traditional approach
            if kpts is None:
                kpts = self.red_kpoints
            
            print('Building kD-tree for kpoints')
            self.kpts = kpts
            self.kpt_tree = build_ktree(kpts)
            self.qidx_in_kpts = find_kpt(self.kpt_tree, kpts)
    
    def _read_dipoles_db(self, ydipdb=None, dip_dir='gw', bands_range=None):
        """
        Read Yambo dipoles database.
        
        Parameters
        ----------
        ydipdb : YamboDipolesDB, optional
            Pre-loaded database. If None, reads from ndb.dipoles file.
        dip_dir : str, optional
            Directory containing dipoles file. Defaults to 'gw'.
        bands_range : list, optional
            Range of bands to load.
        """
        bands_range = bands_range or self.bands_range
        dip_dir_path = os.path.join(self.path, dip_dir)
        
        try:
            ndb_dipoles_fname = os.path.join(dip_dir_path, 'ndb.dipoles')
            if ydipdb:
                self.ydipdb = ydipdb
            else:
                self.ydipdb = YamboDipolesDB(
                    self.ydb, save='', filename=ndb_dipoles_fname, 
                    dip_type='iR', field_dir=[1, 1, 1], project=False, 
                    expand=False, bands_range=bands_range
                )
        except Exception as e:
            print(f"Warning: Could not read dipoles database: {e}")
            print("Continuing without dipoles data - some functionality may be limited")
            self.ydipdb = None
            self.ele_dips = None
            return
        
        # Process dipoles based on spin
        if self.ydipdb.spin == 2:
            self.ele_dips = self.ydipdb.dipoles.conjugate().transpose(1, 2, 3, 4, 0)
        elif self.ydipdb.spin == 1:
            self.ele_dips = self.ydipdb.dipoles.conjugate().transpose(0, 2, 3, 1)
    
    def read_excdb(self, BSE_dir=None):
        """
        Read yambo exciton database for each Q-point.
        
        Parameters
        ----------
        BSE_dir : str, optional
            Directory containing BSE files. If None, uses self.BSE_dir.
        
        Returns
        -------
        tuple
            (bs_bands, BS_eigs, BS_wfcs, excQpt) containing bands involved in BSE,
            eigenenergies, exciton wavefunctions, and Q-points.
        """
        BSE_dir = BSE_dir or self.BSE_dir
        
        bs_bands = []  # bands involved in BSE
        BS_eigs = []   # eigenenergies BSE
        BS_wfcs = []   # exciton wavefunctions
        excQpt = []    # Q-point of BSE -> The q of A^{\lambda Q}_{cvk}
        
        for iq in tqdm(range(self.nibz), desc="Loading Ex-wfcs "):
            try:
                bse_db_iq = YamboExcitonDB.from_db_file(
                    self.ydb, folder=BSE_dir, filename=f'ndb.BS_diago_Q{iq+1}'
                )
            except Exception as e:
                raise IOError(f'Cannot read ndb.BS_diago_Q{iq+1} file: {e}')
            
            bs_bands = bse_db_iq.nbands
            tmp_eigs = bse_db_iq.eigenvalues
            tmp_wfcs = bse_db_iq.get_Akcv()
            tmp_qpt = self.ydb.lat @ bse_db_iq.car_qpoint
            
            BS_eigs.append(tmp_eigs)
            BS_wfcs.append(tmp_wfcs)
            excQpt.append(tmp_qpt)
        
        return (bs_bands, 
                (np.array(BS_eigs) / ha2ev).astype(self.wfdb.wf.dtype),
                np.array(BS_wfcs).astype(self.wfdb.wf.dtype), 
                excQpt)
    
    def read_common_databases(self, latdb=None, wfdb=None, bands_range=None):
        """
        Read common databases used by most optical properties calculations.
        
        Parameters
        ----------
        latdb : YamboLatticeDB, optional
            Pre-loaded lattice database.
        wfdb : YamboWFDB, optional
            Pre-loaded wavefunction database.
        bands_range : list, optional
            Range of bands to load.
        """
        # Read lattice database
        self._read_lattice_db(latdb)
        
        # Read wavefunction database
        self._read_wavefunction_db(wfdb, bands_range)
        
        # Setup k-point mapping
        self._setup_kpoint_mapping()
        
        # Read exciton database
        self.bs_bands, self.BS_eigs, self.BS_wfcs, self.excQpt = self.read_excdb()
        
        # Build k-point tree
        self._build_kpoint_tree()
    
    def enable_geometry_manager(self):
        """
        Enable geometry manager for this instance.
        
        This method switches the instance to use the geometry manager for
        shared geometry data, providing memory and performance benefits.
        
        Returns
        -------
        dict
            Information about the geometry manager benefits.
        """
        if self.use_geometry_manager:
            print("Geometry manager is already enabled for this instance.")
            return {'status': 'already_enabled'}
        
        print("Enabling geometry manager...")
        
        # Store current data if any
        current_ydb = self.ydb
        
        # Enable geometry manager
        self.use_geometry_manager = True
        self.geometry_manager = get_geometry_manager(self.path, 
                                                   os.path.basename(self.SAVE_dir), 
                                                   current_ydb)
        
        # Clear private attributes to force use of geometry manager
        for attr in ['_lat_vecs', '_blat_vecs', '_nibz', '_symm_mats', '_ele_time_rev', '_red_kpoints']:
            if hasattr(self, attr):
                delattr(self, attr)
        
        migration_info = {
            'status': 'enabled',
            'benefits': [
                'Reduced memory usage through shared geometry data',
                'Faster initialization when multiple classes use same calculation',
                'Centralized geometry management',
                'Singleton pattern prevents redundant database reads'
            ],
            'usage': [
                'Geometry data now accessed via properties (lat_vecs, symm_mats, etc.)',
                'Data is shared with other instances using the same calculation path',
                'No code changes needed - same API maintained'
            ],
            'compatibility': 'Full backward compatibility maintained'
        }
        
        print("Geometry Manager Enabled:")
        print("=" * 40)
        for key, value in migration_info.items():
            if isinstance(value, list):
                print(f"{key.replace('_', ' ').title()}:")
                for item in value:
                    print(f"  - {item}")
            else:
                print(f"{key.replace('_', ' ').title()}: {value}")
        
        return migration_info
    
    def get_geometry_manager_info(self):
        """
        Get information about geometry manager usage.
        
        Returns
        -------
        dict
            Information about geometry manager status and benefits.
        """
        info = {
            'geometry_manager_enabled': self.use_geometry_manager,
            'geometry_manager_initialized': False,
            'shared_instances': 0,
            'memory_benefits': 'Not using geometry manager'
        }
        
        if self.use_geometry_manager and self.geometry_manager:
            info['geometry_manager_initialized'] = self.geometry_manager.is_initialized()
            info['memory_benefits'] = 'Sharing geometry data across instances'
            
            # Count shared instances (approximate)
            from yambopy.geometry_manager import GeometryManagerRegistry
            info['shared_instances'] = len(GeometryManagerRegistry._instances)
        
        return info
    
    @abstractmethod
    def compute(self):
        """
        Abstract method for main computation.
        Must be implemented by subclasses.
        """
        pass