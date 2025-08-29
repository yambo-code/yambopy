#
# License-Identifier: GPL
#
# Copyright (C) 2024 The Yambo Team
#
# Authors: RR
#
# This file is part of the yambopy project
#
"""
Abstract geometry manager for centralizing lattice and k-point information.

This module provides a singleton-based geometry manager that centralizes the reading
and management of lattice information, k-points, and related geometric data to 
minimize redundant database reads across different classes.
"""

import os
import numpy as np
from abc import ABC, abstractmethod
from typing import Optional, Dict, Tuple, Any
import warnings

from yambopy import YamboLatticeDB
from yambopy.kpoints import build_ktree, find_kpt

try:
    from pykdtree.kdtree import KDTree 
    # pykdtree is much faster and is recommended
    # pip install pykdtree
    # useful in Dmat computation
except ImportError:
    from scipy.spatial import KDTree

warnings.filterwarnings('ignore')


class GeometryManagerRegistry:
    """
    Registry for managing geometry manager instances.
    
    This class implements a singleton pattern to ensure that geometry data
    is shared across instances when they refer to the same calculation path
    and configuration.
    """
    
    _instances: Dict[str, 'BaseGeometryManager'] = {}
    
    @classmethod
    def get_instance(cls, path: str, save: str = 'SAVE', 
                    latdb: Optional[YamboLatticeDB] = None) -> 'BaseGeometryManager':
        """
        Get or create a geometry manager instance.
        
        Parameters
        ----------
        path : str
            Path to the calculation directory.
        save : str, optional
            SAVE directory name. Defaults to 'SAVE'.
        latdb : YamboLatticeDB, optional
            Pre-loaded lattice database.
            
        Returns
        -------
        BaseGeometryManager
            Geometry manager instance.
        """
        # Create a unique key for this configuration
        key = f"{os.path.abspath(path)}:{save}"
        
        if key not in cls._instances:
            cls._instances[key] = StandardGeometryManager(path, save, latdb)
        
        return cls._instances[key]
    
    @classmethod
    def clear_cache(cls):
        """Clear all cached instances."""
        cls._instances.clear()


class BaseGeometryManager(ABC):
    """
    Abstract base class for geometry managers.
    
    This class defines the interface for geometry managers that handle
    lattice information, k-points, and related geometric data.
    """
    
    def __init__(self, path: str, save: str = 'SAVE', 
                 latdb: Optional[YamboLatticeDB] = None):
        """
        Initialize geometry manager.
        
        Parameters
        ----------
        path : str
            Path to the calculation directory.
        save : str, optional
            SAVE directory name. Defaults to 'SAVE'.
        latdb : YamboLatticeDB, optional
            Pre-loaded lattice database.
        """
        self.path = os.path.abspath(path)
        self.save_dir = os.path.join(self.path, save)
        self._initialized = False
        self._latdb = latdb
        
        # Initialize all geometry-related attributes
        self._init_attributes()
    
    def _init_attributes(self):
        """Initialize all geometry-related attributes."""
        # Lattice information
        self.ydb: Optional[YamboLatticeDB] = None
        self.lat_vecs: Optional[np.ndarray] = None
        self.blat_vecs: Optional[np.ndarray] = None
        self.nibz: Optional[int] = None
        self.symm_mats: Optional[np.ndarray] = None
        self.ele_time_rev: Optional[int] = None
        
        # K-point information
        self.red_kpoints: Optional[np.ndarray] = None
        self.car_kpoints: Optional[np.ndarray] = None
        self.iku_kpoints: Optional[np.ndarray] = None
        self.ibz_kpoints: Optional[np.ndarray] = None
        self.kpt_tree: Optional[Any] = None
        self.qidx_in_kpts: Optional[np.ndarray] = None
        
        # Symmetry information
        self.sym_red: Optional[np.ndarray] = None
        self.kmap: Optional[np.ndarray] = None
    
    @abstractmethod
    def initialize(self) -> None:
        """Initialize the geometry manager by reading necessary databases."""
        pass
    
    @abstractmethod
    def get_lattice_info(self) -> Dict[str, Any]:
        """Get lattice information as a dictionary."""
        pass
    
    @abstractmethod
    def get_kpoint_info(self) -> Dict[str, Any]:
        """Get k-point information as a dictionary."""
        pass
    
    @abstractmethod
    def get_symmetry_info(self) -> Dict[str, Any]:
        """Get symmetry information as a dictionary."""
        pass
    
    def is_initialized(self) -> bool:
        """Check if the geometry manager is initialized."""
        return self._initialized


class StandardGeometryManager(BaseGeometryManager):
    """
    Standard implementation of geometry manager.
    
    This class provides the standard implementation for managing geometry
    information from Yambo databases.
    """
    
    def initialize(self) -> None:
        """Initialize the geometry manager by reading necessary databases."""
        if self._initialized:
            return
        
        print(f"Initializing geometry manager for path: {self.path}")
        
        # Read lattice database
        self._read_lattice_db()
        
        # Setup k-point information
        self._setup_kpoint_info()
        
        # Setup symmetry information
        self._setup_symmetry_info()
        
        # Build k-point tree
        self._build_kpoint_tree()
        
        self._initialized = True
        print("Geometry manager initialized successfully")
    
    def _read_lattice_db(self):
        """Read Yambo Lattice database."""
        try:
            ns_db1_fname = os.path.join(self.save_dir, 'ns.db1')
            if self._latdb:
                if not hasattr(self._latdb, 'ibz_kpoints'):
                    self._latdb.expand_kpoints()
                self.ydb = self._latdb
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
    
    def _setup_kpoint_info(self):
        """Setup k-point information."""
        # Set different k-point representations
        self.iku_kpoints = self.ydb.iku_kpoints
        self.red_kpoints = self.ydb.red_kpoints
        self.car_kpoints = self.ydb.car_kpoints
        self.ibz_kpoints = self.ydb.ibz_kpoints
    
    def _setup_symmetry_info(self):
        """Setup symmetry information."""
        # Setup symmetry operations in reduced coordinates
        temp = np.matmul(self.symm_mats, self.blat_vecs)
        sym_red = np.matmul(self.lat_vecs[None, :, :], temp)
        self.sym_red = np.rint(sym_red).astype(int)
        
        # Set kmap if we have the necessary information
        if hasattr(self.ydb, 'nkBZ') and hasattr(self.ydb, 'kpoints_indexes'):
            kmap = np.zeros((self.ydb.nkBZ, 2), dtype=int)
            kmap[:, 0] = self.ydb.kpoints_indexes
            kmap[:, 1] = self.ydb.symmetry_indexes
            self.kmap = kmap
    
    def _build_kpoint_tree(self, kpts: Optional[np.ndarray] = None):
        """
        Build k-point tree for efficient k-point searching.
        
        Parameters
        ----------
        kpts : array_like, optional
            K-points to build tree from. If None, uses self.red_kpoints.
        """
        if kpts is None:
            kpts = self.red_kpoints
        
        if kpts is not None:
            print('Building kD-tree for kpoints')
            self.kpt_tree = build_ktree(kpts)
            self.qidx_in_kpts = find_kpt(self.kpt_tree, kpts)
    
    def get_lattice_info(self) -> Dict[str, Any]:
        """Get lattice information as a dictionary."""
        if not self._initialized:
            self.initialize()
        
        return {
            'lat_vecs': self.lat_vecs,
            'blat_vecs': self.blat_vecs,
            'nibz': self.nibz,
            'ydb': self.ydb
        }
    
    def get_kpoint_info(self) -> Dict[str, Any]:
        """Get k-point information as a dictionary."""
        if not self._initialized:
            self.initialize()
        
        return {
            'red_kpoints': self.red_kpoints,
            'car_kpoints': self.car_kpoints,
            'iku_kpoints': self.iku_kpoints,
            'ibz_kpoints': self.ibz_kpoints,
            'kpt_tree': self.kpt_tree,
            'qidx_in_kpts': self.qidx_in_kpts
        }
    
    def get_symmetry_info(self) -> Dict[str, Any]:
        """Get symmetry information as a dictionary."""
        if not self._initialized:
            self.initialize()
        
        return {
            'symm_mats': self.symm_mats,
            'ele_time_rev': self.ele_time_rev,
            'sym_red': self.sym_red,
            'kmap': self.kmap
        }
    
    def get_all_info(self) -> Dict[str, Any]:
        """Get all geometry information as a dictionary."""
        if not self._initialized:
            self.initialize()
        
        info = {}
        info.update(self.get_lattice_info())
        info.update(self.get_kpoint_info())
        info.update(self.get_symmetry_info())
        
        return info


def get_geometry_manager(path: str, save: str = 'SAVE', 
                        latdb: Optional[YamboLatticeDB] = None) -> BaseGeometryManager:
    """
    Get a geometry manager instance.
    
    This is the main entry point for getting a geometry manager. It uses
    the registry to ensure that instances are shared when appropriate.
    
    Parameters
    ----------
    path : str
        Path to the calculation directory.
    save : str, optional
        SAVE directory name. Defaults to 'SAVE'.
    latdb : YamboLatticeDB, optional
        Pre-loaded lattice database.
        
    Returns
    -------
    BaseGeometryManager
        Geometry manager instance.
        
    Examples
    --------
    >>> from yambopy.geometry_manager import get_geometry_manager
    >>> 
    >>> # Get geometry manager for a calculation
    >>> geom_mgr = get_geometry_manager('/path/to/calculation')
    >>> 
    >>> # Get lattice information
    >>> lattice_info = geom_mgr.get_lattice_info()
    >>> lat_vecs = lattice_info['lat_vecs']
    >>> 
    >>> # Get k-point information
    >>> kpoint_info = geom_mgr.get_kpoint_info()
    >>> red_kpoints = kpoint_info['red_kpoints']
    """
    return GeometryManagerRegistry.get_instance(path, save, latdb)


def clear_geometry_cache():
    """
    Clear the geometry manager cache.
    
    This function clears all cached geometry manager instances, which can
    be useful for freeing memory or when working with different calculations.
    """
    GeometryManagerRegistry.clear_cache()