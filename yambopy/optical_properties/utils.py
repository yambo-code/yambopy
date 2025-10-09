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
Utility functions for optical properties calculations.
"""

import os
import numpy as np
from typing import Optional, Tuple, List, Union

from yambopy.letzelphc_interface.lelphcdb import LetzElphElectronPhononDB
from yambopy.units import ha2ev


def validate_path(path: Optional[str] = None) -> str:
    """
    Validate and return absolute path.
    
    Parameters
    ----------
    path : str, optional
        Input path. If None, returns current directory.
    
    Returns
    -------
    str
        Validated absolute path.
    """
    if path is None:
        return os.getcwd()
    return os.path.abspath(path)


def setup_directories(base_path: str, **dir_configs) -> dict:
    """
    Setup multiple directories based on configuration.
    
    Parameters
    ----------
    base_path : str
        Base path for all directories.
    **dir_configs : dict
        Directory configurations as key-value pairs.
    
    Returns
    -------
    dict
        Dictionary of directory paths.
    
    Examples
    --------
    >>> dirs = setup_directories('/path/to/calc', 
    ...                          SAVE='SAVE', BSE='bse', DIP='gw')
    >>> print(dirs['SAVE'])  # '/path/to/calc/SAVE'
    """
    directories = {}
    for dir_name, dir_value in dir_configs.items():
        directories[dir_name] = os.path.join(base_path, dir_value)
    return directories


def process_bands_range(bands_range: Optional[List[int]] = None, 
                       total_bands: Optional[int] = None) -> Tuple[List[int], int, int]:
    """
    Process and validate bands range.
    
    Parameters
    ----------
    bands_range : list, optional
        Range of bands [start, end). If None, uses all bands.
    total_bands : int, optional
        Total number of bands available.
    
    Returns
    -------
    tuple
        (processed_bands_range, start_idx, num_bands)
    """
    if not bands_range:
        if total_bands is None:
            return [], 0, 0
        bands_range = [0, total_bands]
    
    start_idx = min(bands_range)
    end_idx = max(bands_range)
    num_bands = end_idx - start_idx
    
    return bands_range, start_idx, num_bands


def read_lelph_database(lelph_dir: str, 
                       lelph_db: Optional[LetzElphElectronPhononDB] = None) -> Optional[LetzElphElectronPhononDB]:
    """
    Read LetzElPhC database with graceful error handling.
    
    Parameters
    ----------
    lelph_dir : str
        Directory containing LELPH files.
    lelph_db : LetzElphElectronPhononDB, optional
        Pre-loaded database. If None, reads from file.
    
    Returns
    -------
    LetzElphElectronPhononDB or None
        Loaded database, or None if reading failed.
    """
    if lelph_db is not None:
        return lelph_db
    
    try:
        ndb_lelph_fname = os.path.join(lelph_dir, 'ndb.elph')
        return LetzElphElectronPhononDB(filename=ndb_lelph_fname)
    except Exception as e:
        print(f"Warning: Could not read LELPH database: {e}")
        print("Continuing without LELPH data - some functionality may be limited")
        return None


def compute_symmetry_matrices(symm_mats: np.ndarray, 
                             blat_vecs: np.ndarray, 
                             lat_vecs: np.ndarray) -> np.ndarray:
    """
    Compute symmetry matrices in reduced coordinates.
    
    Parameters
    ----------
    symm_mats : np.ndarray
        Symmetry matrices in Cartesian coordinates.
    blat_vecs : np.ndarray
        Reciprocal lattice vectors.
    lat_vecs : np.ndarray
        Real space lattice vectors.
    
    Returns
    -------
    np.ndarray
        Symmetry matrices in reduced coordinates.
    """
    temp = np.matmul(symm_mats, blat_vecs)
    sym_red = np.matmul(lat_vecs[None, :, :], temp)
    return np.rint(sym_red).astype(int)


def create_kpoint_mapping(nkBZ: int, 
                         kpoints_indexes: np.ndarray, 
                         symmetry_indexes: np.ndarray) -> np.ndarray:
    """
    Create k-point mapping array.
    
    Parameters
    ----------
    nkBZ : int
        Number of k-points in Brillouin zone.
    kpoints_indexes : np.ndarray
        K-point indices.
    symmetry_indexes : np.ndarray
        Symmetry operation indices.
    
    Returns
    -------
    np.ndarray
        K-point mapping array of shape (nkBZ, 2).
    """
    kmap = np.zeros((nkBZ, 2), dtype=int)
    kmap[:, 0] = kpoints_indexes
    kmap[:, 1] = symmetry_indexes
    return kmap


def process_dipoles_by_spin(dipoles: np.ndarray, spin: int) -> np.ndarray:
    """
    Process dipole matrix elements based on spin configuration.
    
    Parameters
    ----------
    dipoles : np.ndarray
        Raw dipole matrix elements.
    spin : int
        Spin configuration (1 or 2).
    
    Returns
    -------
    np.ndarray
        Processed dipole matrix elements.
    """
    if spin == 2:
        return dipoles.conjugate().transpose(1, 2, 3, 4, 0)
    elif spin == 1:
        return dipoles.conjugate().transpose(0, 2, 3, 1)
    else:
        raise ValueError(f"Unsupported spin configuration: {spin}")


def safe_file_operation(operation_func, filename: str, error_msg: str):
    """
    Safely perform file operations with error handling.
    
    Parameters
    ----------
    operation_func : callable
        Function to perform the file operation.
    filename : str
        Name of the file being operated on.
    error_msg : str
        Error message to display if operation fails.
    
    Returns
    -------
    Any
        Result of the operation function.
    
    Raises
    ------
    IOError
        If the file operation fails.
    """
    try:
        return operation_func()
    except Exception as e:
        raise IOError(f'{error_msg}: {e}')


def convert_energy_units(energies: np.ndarray, 
                        from_unit: str = 'Ry', 
                        to_unit: str = 'Ha',
                        dtype: Optional[np.dtype] = None) -> np.ndarray:
    """
    Convert energy units.
    
    Parameters
    ----------
    energies : np.ndarray
        Energy values to convert.
    from_unit : str, optional
        Source unit ('Ry', 'Ha', 'eV'). Defaults to 'Ry'.
    to_unit : str, optional
        Target unit ('Ry', 'Ha', 'eV'). Defaults to 'Ha'.
    dtype : np.dtype, optional
        Target data type. If None, preserves original dtype.
    
    Returns
    -------
    np.ndarray
        Converted energy values.
    """
    conversion_factors = {
        ('Ry', 'Ha'): 1.0 / 2.0,
        ('Ha', 'Ry'): 2.0,
        ('Ry', 'eV'): ha2ev / 2.0,
        ('Ha', 'eV'): ha2ev,
        ('eV', 'Ha'): 1.0 / ha2ev,
        ('eV', 'Ry'): 2.0 / ha2ev,
    }
    
    if from_unit == to_unit:
        factor = 1.0
    else:
        factor = conversion_factors.get((from_unit, to_unit))
        if factor is None:
            raise ValueError(f"Unsupported conversion: {from_unit} -> {to_unit}")
    
    converted = energies * factor
    
    if dtype is not None:
        converted = converted.astype(dtype)
    
    return converted


def validate_matrix_dimensions(matrix: np.ndarray, 
                              expected_shape: Optional[Tuple[int, ...]] = None,
                              min_dims: Optional[int] = None,
                              max_dims: Optional[int] = None) -> None:
    """
    Validate matrix dimensions.
    
    Parameters
    ----------
    matrix : np.ndarray
        Matrix to validate.
    expected_shape : tuple, optional
        Expected shape. If None, no shape validation.
    min_dims : int, optional
        Minimum number of dimensions.
    max_dims : int, optional
        Maximum number of dimensions.
    
    Raises
    ------
    ValueError
        If matrix dimensions are invalid.
    """
    if expected_shape is not None and matrix.shape != expected_shape:
        raise ValueError(f"Matrix shape {matrix.shape} does not match expected {expected_shape}")
    
    if min_dims is not None and matrix.ndim < min_dims:
        raise ValueError(f"Matrix has {matrix.ndim} dimensions, minimum required: {min_dims}")
    
    if max_dims is not None and matrix.ndim > max_dims:
        raise ValueError(f"Matrix has {matrix.ndim} dimensions, maximum allowed: {max_dims}")


def create_progress_bar(iterable, desc: str = "Processing", **kwargs):
    """
    Create a progress bar with consistent formatting.
    
    Parameters
    ----------
    iterable : iterable
        Iterable to wrap with progress bar.
    desc : str, optional
        Description for the progress bar.
    **kwargs
        Additional arguments for tqdm.
    
    Returns
    -------
    tqdm
        Progress bar iterator.
    """
    from tqdm import tqdm
    return tqdm(iterable, desc=desc, **kwargs)