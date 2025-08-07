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
                 bands_range=None, BSE_dir='bse', save_files=True):
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
        """
        self.path = path if path is not None else os.getcwd()
        self.SAVE_dir = os.path.join(self.path, save)
        self.BSE_dir = os.path.join(self.path, BSE_dir)
        self.save_files = save_files
        self.bands_range = bands_range if bands_range is not None else []
        
        # Initialize databases
        self.latdb = latdb
        self.wfdb = wfdb
        
        # Initialize common attributes
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
    
    def _read_lattice_db(self, latdb=None):
        """
        Read Yambo Lattice database.
        
        Parameters
        ----------
        latdb : YamboLatticeDB, optional
            Pre-loaded database. If None, reads from ns.db1 file.
        """
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
        
        # Set common lattice properties
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
            K-points to build tree from. If None, uses self.ydb.red_kpoints.
        """
        if kpts is None:
            kpts = self.ydb.red_kpoints
        
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
            raise IOError(f'Cannot read ndb.dipoles file: {e}')
        
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
    
    @abstractmethod
    def compute(self):
        """
        Abstract method for main computation.
        Must be implemented by subclasses.
        """
        pass