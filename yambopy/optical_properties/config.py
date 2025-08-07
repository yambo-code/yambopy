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
Configuration constants and settings for optical properties calculations.
"""

import numpy as np

# Default tolerances
DEFAULT_TOLERANCE = 1e-6
SYMMETRY_TOLERANCE = 1e-4
ENERGY_TOLERANCE = 1e-9

# Default file names
DEFAULT_FILES = {
    'lattice_db': 'ns.db1',
    'wavefunction_db': 'ns.wf',
    'dipoles_db': 'ndb.dipoles',
    'elph_db': 'ndb.elph',
    'dmats_db': 'ndb.Dmats',
}

# Default directory names
DEFAULT_DIRECTORIES = {
    'SAVE': 'SAVE',
    'BSE': 'bse',
    'DIP': 'gw',
    'LELPH': 'lelph',
}

# Default computation parameters
DEFAULT_PARAMS = {
    'temperature': 20.0,  # Kelvin
    'broadening': 0.00124,  # Ha
    'npol': 2,
    'ph_threshold': 1e-9,  # Ry
    'degen_threshold': 0.001,  # eV
}

# Default dipole parameters
DIPOLE_PARAMS = {
    'dip_type': 'iR',
    'field_dir': [1, 1, 1],
    'project': False,
    'expand': False,
}

# Progress bar settings
PROGRESS_BAR_SETTINGS = {
    'desc_loading': "Loading Ex-wfcs",
    'desc_computing': "Computing",
    'desc_luminescence': "Luminescence",
    'desc_processing': "Processing",
}

# Data type preferences
DTYPE_PREFERENCES = {
    'complex': np.complex128,
    'real': np.float64,
    'integer': np.int32,
}

# File extensions
FILE_EXTENSIONS = {
    'numpy': '.npy',
    'text': '.dat',
    'hdf5': '.h5',
    'netcdf': '.nc',
}

# Unit conversion factors (relative to Hartree atomic units)
UNIT_CONVERSIONS = {
    'Ha_to_eV': 27.211386245988,
    'Ry_to_Ha': 0.5,
    'eV_to_Ha': 1.0 / 27.211386245988,
    'Ha_to_Ry': 2.0,
}

# Symmetry operation types
SYMMETRY_TYPES = {
    'identity': 'E',
    'rotation': 'C',
    'reflection': 'σ',
    'inversion': 'i',
    'improper_rotation': 'S',
}

# Common validation patterns
VALIDATION_PATTERNS = {
    'positive_integer': lambda x: isinstance(x, int) and x > 0,
    'non_negative_real': lambda x: isinstance(x, (int, float)) and x >= 0,
    'positive_real': lambda x: isinstance(x, (int, float)) and x > 0,
    'valid_range': lambda x: isinstance(x, (list, tuple)) and len(x) == 2 and x[0] < x[1],
}

# Error messages
ERROR_MESSAGES = {
    'file_not_found': "Cannot read {filename}: {error}",
    'invalid_shape': "Matrix shape {actual} does not match expected {expected}",
    'invalid_range': "Invalid range: {range}. Expected [start, end) with start < end",
    'missing_attribute': "Required attribute '{attr}' not found in {obj}",
    'computation_failed': "Computation failed: {error}",
}

# Logging configuration
LOGGING_CONFIG = {
    'format': '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    'level': 'INFO',
    'file_level': 'DEBUG',
}

# Performance settings
PERFORMANCE_SETTINGS = {
    'use_numba': True,
    'parallel_threshold': 1000,  # Minimum size for parallel operations
    'chunk_size': 100,  # Default chunk size for batch operations
    'memory_limit_gb': 8,  # Memory limit for large arrays
}

# Plotting defaults (if matplotlib is available)
PLOT_DEFAULTS = {
    'figsize': (10, 6),
    'dpi': 300,
    'style': 'seaborn-v0_8',
    'colors': ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd'],
    'linewidth': 2,
    'markersize': 6,
}

def get_default_config():
    """
    Get default configuration dictionary.
    
    Returns
    -------
    dict
        Default configuration settings.
    """
    return {
        'tolerances': {
            'default': DEFAULT_TOLERANCE,
            'symmetry': SYMMETRY_TOLERANCE,
            'energy': ENERGY_TOLERANCE,
        },
        'files': DEFAULT_FILES,
        'directories': DEFAULT_DIRECTORIES,
        'parameters': DEFAULT_PARAMS,
        'dipole': DIPOLE_PARAMS,
        'progress': PROGRESS_BAR_SETTINGS,
        'dtypes': DTYPE_PREFERENCES,
        'extensions': FILE_EXTENSIONS,
        'units': UNIT_CONVERSIONS,
        'validation': VALIDATION_PATTERNS,
        'errors': ERROR_MESSAGES,
        'logging': LOGGING_CONFIG,
        'performance': PERFORMANCE_SETTINGS,
        'plotting': PLOT_DEFAULTS,
    }

def validate_config(config):
    """
    Validate configuration dictionary.
    
    Parameters
    ----------
    config : dict
        Configuration dictionary to validate.
    
    Returns
    -------
    bool
        True if configuration is valid.
    
    Raises
    ------
    ValueError
        If configuration is invalid.
    """
    required_keys = ['tolerances', 'files', 'directories', 'parameters']
    
    for key in required_keys:
        if key not in config:
            raise ValueError(f"Missing required configuration key: {key}")
    
    # Validate tolerance values
    for tol_name, tol_value in config['tolerances'].items():
        if not isinstance(tol_value, (int, float)) or tol_value <= 0:
            raise ValueError(f"Invalid tolerance value for {tol_name}: {tol_value}")
    
    return True