"""
Comprehensive optical activity database for all crystallographic point groups.

This module contains selection rules for IR, Raman, and electric dipole transitions
for all 32 crystallographic point groups. The data is organized by point group
symbol and includes the irreducible representations that are active for each
type of optical transition.

References:
- Tinkham, M. "Group Theory and Quantum Mechanics" (1964)
- Cotton, F.A. "Chemical Applications of Group Theory" (1990)
- Dresselhaus, M.S. "Group Theory: Application to the Physics of Condensed Matter" (2008)
- International Tables for Crystallography, Volume A (2016)
"""

# Comprehensive optical activity database for all 32 crystallographic point groups
OPTICAL_ACTIVITY_DATABASE = {
    
    # ========================================================================
    # TRICLINIC SYSTEM
    # ========================================================================
    
    'C1': {  # Point group 1
        'crystal_system': 'triclinic',
        'ir_active': ['A'],
        'raman_active': ['A'],
        'electric_dipole': ['A'],
        'notes': 'All transitions allowed - no symmetry restrictions'
    },
    
    'Ci': {  # Point group -1
        'crystal_system': 'triclinic',
        'ir_active': ['Au'],
        'raman_active': ['Ag'],
        'electric_dipole': ['Au'],
        'notes': 'Mutual exclusion rule: IR and Raman active modes are different'
    },
    
    # ========================================================================
    # MONOCLINIC SYSTEM
    # ========================================================================
    
    'C2': {  # Point group 2
        'crystal_system': 'monoclinic',
        'ir_active': ['A', 'B'],
        'raman_active': ['A', 'B'],
        'electric_dipole': ['A', 'B'],
        'notes': 'All irreps are both IR and Raman active'
    },
    
    'Cs': {  # Point group m
        'crystal_system': 'monoclinic',
        'ir_active': ['A\'', 'A\'\''],
        'raman_active': ['A\'', 'A\'\''],
        'electric_dipole': ['A\'', 'A\'\''],
        'notes': 'All irreps are both IR and Raman active'
    },
    
    'C2h': {  # Point group 2/m
        'crystal_system': 'monoclinic',
        'ir_active': ['Au', 'Bu'],
        'raman_active': ['Ag', 'Bg'],
        'electric_dipole': ['Au', 'Bu'],
        'notes': 'Mutual exclusion rule applies'
    },
    
    # ========================================================================
    # ORTHORHOMBIC SYSTEM
    # ========================================================================
    
    'C2v': {  # Point group mm2
        'crystal_system': 'orthorhombic',
        'ir_active': ['A1', 'B1', 'B2'],
        'raman_active': ['A1', 'A2', 'B1', 'B2'],
        'electric_dipole': ['A1', 'B1', 'B2'],
        'notes': 'A2 is only Raman active'
    },
    
    'D2': {  # Point group 222
        'crystal_system': 'orthorhombic',
        'ir_active': ['B1', 'B2', 'B3'],
        'raman_active': ['A', 'B1', 'B2', 'B3'],
        'electric_dipole': ['B1', 'B2', 'B3'],
        'notes': 'A is only Raman active'
    },
    
    'D2h': {  # Point group mmm
        'crystal_system': 'orthorhombic',
        'ir_active': ['B1u', 'B2u', 'B3u'],
        'raman_active': ['Ag', 'B1g', 'B2g', 'B3g'],
        'electric_dipole': ['B1u', 'B2u', 'B3u'],
        'notes': 'Mutual exclusion rule applies'
    },
    
    # ========================================================================
    # TETRAGONAL SYSTEM
    # ========================================================================
    
    'C4': {  # Point group 4
        'crystal_system': 'tetragonal',
        'ir_active': ['A', 'E'],
        'raman_active': ['A', 'B', 'E'],
        'electric_dipole': ['A', 'E'],
        'notes': 'B is only Raman active'
    },
    
    'S4': {  # Point group -4
        'crystal_system': 'tetragonal',
        'ir_active': ['E'],
        'raman_active': ['A', 'B', 'E'],
        'electric_dipole': ['E'],
        'notes': 'A and B are only Raman active'
    },
    
    'C4h': {  # Point group 4/m
        'crystal_system': 'tetragonal',
        'ir_active': ['Au', 'Eu'],
        'raman_active': ['Ag', 'Bg', 'Eg'],
        'electric_dipole': ['Au', 'Eu'],
        'notes': 'Mutual exclusion rule applies'
    },
    
    'C4v': {  # Point group 4mm
        'crystal_system': 'tetragonal',
        'ir_active': ['A1', 'E'],
        'raman_active': ['A1', 'A2', 'B1', 'B2', 'E'],
        'electric_dipole': ['A1', 'E'],
        'notes': 'A2, B1, B2 are only Raman active'
    },
    
    'D4': {  # Point group 422
        'crystal_system': 'tetragonal',
        'ir_active': ['E'],
        'raman_active': ['A1', 'A2', 'B1', 'B2', 'E'],
        'electric_dipole': ['E'],
        'notes': 'A1, A2, B1, B2 are only Raman active'
    },
    
    'D2d': {  # Point group -42m
        'crystal_system': 'tetragonal',
        'ir_active': ['E'],
        'raman_active': ['A1', 'A2', 'B1', 'B2', 'E'],
        'electric_dipole': ['E'],
        'notes': 'A1, A2, B1, B2 are only Raman active'
    },
    
    'D4h': {  # Point group 4/mmm
        'crystal_system': 'tetragonal',
        'ir_active': ['A2u', 'Eu'],
        'raman_active': ['A1g', 'B1g', 'B2g', 'Eg'],
        'electric_dipole': ['A2u', 'Eu'],
        'notes': 'Mutual exclusion rule applies'
    },
    
    # ========================================================================
    # TRIGONAL SYSTEM
    # ========================================================================
    
    'C3': {  # Point group 3
        'crystal_system': 'trigonal',
        'ir_active': ['A', 'E'],
        'raman_active': ['A', 'E'],
        'electric_dipole': ['A', 'E'],
        'notes': 'All irreps are both IR and Raman active'
    },
    
    'C3i': {  # Point group -3
        'crystal_system': 'trigonal',
        'ir_active': ['Au', 'Eu'],
        'raman_active': ['Ag', 'Eg'],
        'electric_dipole': ['Au', 'Eu'],
        'notes': 'Mutual exclusion rule applies'
    },
    
    'C3v': {  # Point group 3m
        'crystal_system': 'trigonal',
        'ir_active': ['A1', 'E'],
        'raman_active': ['A1', 'A2', 'E'],
        'electric_dipole': ['A1', 'E'],
        'notes': 'A2 is only Raman active'
    },
    
    'D3': {  # Point group 32
        'crystal_system': 'trigonal',
        'ir_active': ['E'],
        'raman_active': ['A1', 'A2', 'E'],
        'electric_dipole': ['E'],
        'notes': 'A1, A2 are only Raman active'
    },
    
    'D3d': {  # Point group -3m
        'crystal_system': 'trigonal',
        'ir_active': ['A2u', 'Eu'],
        'raman_active': ['A1g', 'Eg'],
        'electric_dipole': ['A2u', 'Eu'],
        'notes': 'Mutual exclusion rule applies'
    },
    
    # ========================================================================
    # HEXAGONAL SYSTEM
    # ========================================================================
    
    'C6': {  # Point group 6
        'crystal_system': 'hexagonal',
        'ir_active': ['A', 'E1'],
        'raman_active': ['A', 'B', 'E1', 'E2'],
        'electric_dipole': ['A', 'E1'],
        'notes': 'B and E2 are only Raman active'
    },
    
    'C3h': {  # Point group -6
        'crystal_system': 'hexagonal',
        'ir_active': ['A\'\'', 'E\'\''],
        'raman_active': ['A\'', 'E\'', 'E\'\''],
        'electric_dipole': ['A\'\'', 'E\'\''],
        'notes': 'A\' and E\' are only Raman active'
    },
    
    'C6h': {  # Point group 6/m
        'crystal_system': 'hexagonal',
        'ir_active': ['A2u', 'E1u'],
        'raman_active': ['A1g', 'E1g', 'E2g'],
        'electric_dipole': ['A2u', 'E1u'],
        'notes': 'Mutual exclusion rule applies'
    },
    
    'C6v': {  # Point group 6mm
        'crystal_system': 'hexagonal',
        'ir_active': ['A1', 'E1'],
        'raman_active': ['A1', 'A2', 'B1', 'B2', 'E1', 'E2'],
        'electric_dipole': ['A1', 'E1'],
        'notes': 'A2, B1, B2, E2 are only Raman active'
    },
    
    'D6': {  # Point group 622
        'crystal_system': 'hexagonal',
        'ir_active': ['E1'],
        'raman_active': ['A1', 'A2', 'B1', 'B2', 'E1', 'E2'],
        'electric_dipole': ['E1'],
        'notes': 'A1, A2, B1, B2, E2 are only Raman active'
    },
    
    'D3h': {  # Point group -6m2
        'crystal_system': 'hexagonal',
        'ir_active': ['A2\'\'', 'E\'\''],
        'raman_active': ['A1\'', 'E\''],
        'electric_dipole': ['A2\'\'', 'E\'\''],
        'notes': 'Mutual exclusion rule applies'
    },
    
    'D6h': {  # Point group 6/mmm
        'crystal_system': 'hexagonal',
        'ir_active': ['A2u', 'E1u'],
        'raman_active': ['A1g', 'E1g', 'E2g'],
        'electric_dipole': ['A2u', 'E1u'],
        'notes': 'Mutual exclusion rule applies - validated with hBN'
    },
    
    # ========================================================================
    # CUBIC SYSTEM
    # ========================================================================
    
    'T': {  # Point group 23
        'crystal_system': 'cubic',
        'ir_active': ['T'],
        'raman_active': ['A', 'E', 'T'],
        'electric_dipole': ['T'],
        'notes': 'A and E are only Raman active'
    },
    
    'Th': {  # Point group m-3
        'crystal_system': 'cubic',
        'ir_active': ['Tu'],
        'raman_active': ['Ag', 'Eg', 'Tg'],
        'electric_dipole': ['Tu'],
        'notes': 'Mutual exclusion rule applies'
    },
    
    'Td': {  # Point group -43m
        'crystal_system': 'cubic',
        'ir_active': ['T2'],
        'raman_active': ['A1', 'E', 'T2'],
        'electric_dipole': ['T2'],
        'notes': 'A1 and E are only Raman active'
    },
    
    'O': {  # Point group 432
        'crystal_system': 'cubic',
        'ir_active': ['T1'],
        'raman_active': ['A1', 'A2', 'E', 'T2'],
        'electric_dipole': ['T1'],
        'notes': 'A1, A2, E, T2 are only Raman active'
    },
    
    'Oh': {  # Point group m-3m
        'crystal_system': 'cubic',
        'ir_active': ['T1u'],
        'raman_active': ['A1g', 'Eg', 'T2g'],
        'electric_dipole': ['T1u'],
        'notes': 'Mutual exclusion rule applies'
    },
}

# Alternative notation mappings for common point group symbols
POINT_GROUP_ALIASES = {
    # Hermann-Mauguin to Schoenflies mapping
    '1': 'C1',
    '-1': 'Ci',
    '2': 'C2',
    'm': 'Cs',
    '2/m': 'C2h',
    'mm2': 'C2v',
    '222': 'D2',
    'mmm': 'D2h',
    '4': 'C4',
    '-4': 'S4',
    '4/m': 'C4h',
    '4mm': 'C4v',
    '422': 'D4',
    '-42m': 'D2d',
    '4/mmm': 'D4h',
    '3': 'C3',
    '-3': 'C3i',
    '3m': 'C3v',
    '32': 'D3',
    '-3m': 'D3d',
    '6': 'C6',
    '-6': 'C3h',
    '6/m': 'C6h',
    '6mm': 'C6v',
    '622': 'D6',
    '-6m2': 'D3h',
    '6/mmm': 'D6h',
    '23': 'T',
    'm-3': 'Th',
    '-43m': 'Td',
    '432': 'O',
    'm-3m': 'Oh',
}

def get_optical_activity(point_group_symbol):
    """
    Get optical activity information for a given point group.
    
    Parameters
    ----------
    point_group_symbol : str
        Point group symbol in either Hermann-Mauguin or Schoenflies notation
        
    Returns
    -------
    dict or None
        Dictionary containing optical activity information, or None if not found
        
    Examples
    --------
    >>> activity = get_optical_activity('6/mmm')  # D6h
    >>> print(activity['ir_active'])
    ['A2u', 'E1u']
    
    >>> activity = get_optical_activity('D6h')  # Alternative notation
    >>> print(activity['raman_active'])
    ['A1g', 'E1g', 'E2g']
    """
    # Try direct lookup first
    if point_group_symbol in OPTICAL_ACTIVITY_DATABASE:
        return OPTICAL_ACTIVITY_DATABASE[point_group_symbol].copy()
    
    # Try alias lookup
    if point_group_symbol in POINT_GROUP_ALIASES:
        schoenflies_symbol = POINT_GROUP_ALIASES[point_group_symbol]
        if schoenflies_symbol in OPTICAL_ACTIVITY_DATABASE:
            return OPTICAL_ACTIVITY_DATABASE[schoenflies_symbol].copy()
    
    # Try reverse alias lookup (Schoenflies to Hermann-Mauguin)
    for hm_symbol, sf_symbol in POINT_GROUP_ALIASES.items():
        if point_group_symbol == sf_symbol and hm_symbol in OPTICAL_ACTIVITY_DATABASE:
            return OPTICAL_ACTIVITY_DATABASE[hm_symbol].copy()
    
    return None

def analyze_optical_activity(point_group_symbol, irrep_multiplicities):
    """
    Analyze optical activity based on irreducible representations.
    
    Parameters
    ----------
    point_group_symbol : str
        Point group symbol
    irrep_multiplicities : list
        List of (irrep_label, multiplicity) tuples
        
    Returns
    -------
    str
        Activity description
        
    Examples
    --------
    >>> activity = analyze_optical_activity('6/mmm', [('A1g', 1), ('E1u', 2)])
    >>> print(activity)
    'IR active, electric dipole allowed'
    """
    activity_data = get_optical_activity(point_group_symbol)
    
    if activity_data is None:
        return "activity unknown (point group not in database)"
    
    present_irreps = [label for label, mult in irrep_multiplicities]
    activities = []
    
    # Check IR activity
    if any(irrep in present_irreps for irrep in activity_data['ir_active']):
        activities.append("IR active")
    
    # Check Raman activity
    if any(irrep in present_irreps for irrep in activity_data['raman_active']):
        activities.append("Raman active")
    
    # Check electric dipole activity
    if any(irrep in present_irreps for irrep in activity_data['electric_dipole']):
        activities.append("electric dipole allowed")
    
    if not activities:
        activities.append("optically inactive")
    
    return ", ".join(activities)

def get_crystal_system_info():
    """
    Get summary information about crystal systems and their point groups.
    
    Returns
    -------
    dict
        Dictionary with crystal system information
    """
    crystal_systems = {}
    
    for pg_symbol, data in OPTICAL_ACTIVITY_DATABASE.items():
        crystal_system = data['crystal_system']
        if crystal_system not in crystal_systems:
            crystal_systems[crystal_system] = []
        crystal_systems[crystal_system].append(pg_symbol)
    
    return crystal_systems

def validate_database():
    """
    Validate the optical activity database for consistency.
    
    Returns
    -------
    dict
        Validation results
    """
    results = {
        'total_point_groups': len(OPTICAL_ACTIVITY_DATABASE),
        'crystal_systems': {},
        'issues': []
    }
    
    # Count point groups per crystal system
    for pg_symbol, data in OPTICAL_ACTIVITY_DATABASE.items():
        crystal_system = data['crystal_system']
        if crystal_system not in results['crystal_systems']:
            results['crystal_systems'][crystal_system] = 0
        results['crystal_systems'][crystal_system] += 1
        
        # Check for required keys
        required_keys = ['crystal_system', 'ir_active', 'raman_active', 'electric_dipole']
        for key in required_keys:
            if key not in data:
                results['issues'].append(f"Missing key '{key}' in {pg_symbol}")
    
    # Check if we have all 32 point groups
    if results['total_point_groups'] != 32:
        results['issues'].append(f"Expected 32 point groups, found {results['total_point_groups']}")
    
    return results

