#
# This file is part of yambopy
#
import numpy as np
from qepy import xml

"""
This file contains a list of functions to deal with Bravais lattices
"""

def get_ibrav(datafile_xml):
    """
    Get ibrav number (as string) from QE data-file-schema.xml
    
    [NB] datafile_xml must be in the form:
    
        datafile_xml = ET.parse( "PREFIX.save/data-file-schema.xml" )
    """
    
    ibrav = xml.get_xml_attrib(datafile_xml,'atomic_structure','bravais_index',repeated=True)
    return ibrav

def lattice_dictionary(ibrav):
    """
    Dictionary with text descriptions for ibrav numbers
    """    
    bravais_lattices = {
    0 : 'free',
    1 : 'cubic P (sc)',
    2 : 'cubic F (fcc)',
    3 : 'cubic I (bcc)',
    4 : 'Hexagonal and Trigonal P',
    5 : 'Trigonal R, 3fold axis c',
    6 : 'Tetragonal P (st)',
    7 : 'Tet ragonal I (bct)',
    8 : 'Orthorhombic P',
    9 : 'Orthorhombic base-centered(bco)',
    10 : 'Orthorhombic face-centered',
    11 : 'Orthorhombic body-centered',
    12 : 'Monoclinic P, unique axis c',
    13 : 'Monoclinic base-centered',
    14 : 'Triclinic'
    }
    try: bravais_lattices[ibrav]
    except KeyError: print("This lattice type is not implemented ):")
    else: return bravais_lattices[ibrav]
    
def lattice_type(ibrav,cell):
    """
    Input: specific lattice type (ibrav) and cell parameters (as list)
    Output: lattice vectors 
    """
    ibrav = int(ibrav)
    lattice = lattice_dictionary(ibrav)
    sqrt = np.sqrt
    
    if lattice=='Hexagonal and Trigonal P':
        a    = cell[0]
        c_a  = cell[2] #c/a
        vecs = a*np.array([ [  1.,          0.,  0. ], 
                            [-0.5, sqrt(3.)/2.,  0. ],
                            [  0.,          0., c_a ] ])
    if lattice=='free':
        vecs = np.array([ datafile_xml.find("CELL/DIRECT_LATTICE_VECTORS/a1").text.split(),
                          datafile_xml.find("CELL/DIRECT_LATTICE_VECTORS/a2").text.split(),
                          datafile_xml.find("CELL/DIRECT_LATTICE_VECTORS/a3").text.split() ]).astype(float)
    if lattice=='cubic P (sc)':
        a    = cell[0]
        vecs = a*np.array([ [ 1., 0., 0. ],
                            [ 0., 1., 0. ],
                            [ 0., 0., 1. ] ])
    if lattice== 'cubic F (fcc)':
        a    = cell[0]
        vecs = a/2.*np.array([ [-1., 0., 1. ],
                               [ 0., 1., 1. ],
                               [-1., 1., 0. ] ]) 
    if lattice=='cubic I (bcc)':
        a    = cell[0]
        vecs = a/2.*np.array([ [ 1., 1., 1. ],
                               [-1., 1., 1. ],
                               [-1.,-1., 0. ] ])
    if lattice=='Trigonal R, 3fold axis c':
        a    = cell[0]
        c_g  = cell[3] #cos(gamma)
        tx, ty, tz = [ sqrt((1.-c_g)/2.), sqrt((1.-c_g)/6.), sqrt((1.+2*c_g)/3.) ]
        vecs = a*np.array([ [ tx,  -ty, tz ],
                            [ 0., 2*ty, tz ],
                            [-tx,  -ty, tz ] ]) 
    if lattice=='Tetragonal P (st)':
        a    = cell[0]
        c_a  = cell[2]
        vecs = a*np.array([ [ 1., 0.,  0. ],
                            [ 0., 1.,  0. ],
                            [ 0., 0., c_a ] ])
    if lattice=='Tetragonal I (bct)':
        a    = cell[0]
        c_a  = cell[2]
        vecs = a/2.*np.array([ [ 1.,-1., c_a ],
                               [ 1., 1., c_a ],
                               [-1.,-1., c_a ] ])
    if lattice=='Orthorhombic P':
        a    = cell[0]
        b_a  = cell[1] #b/a
        c_a  = cell[2]
        vecs = a*np.array([ [ 1.,  0.,  0. ],
                            [ 0., b_a,  0. ],
                            [ 0.,  0., c_a ] ])
    if lattice=='Orthorhombic base-centered(bco)':
        a    = cell[0]
        b_a  = cell[1]
        c_a  = cell[2]
        vecs = a*np.array([ [ 1./2., b_a/2.,  0. ],
                            [-1./2., b_a/2.,  0. ],
                            [    0.,     0., c_a ] ])
    if lattice=='Orthorhombic face-centered':
        a    = cell[0]
        b_a  = cell[1]
        c_a  = cell[2]
        vecs = a/2.*np.array([ [ 1.,  0., c_a ],
                               [ 1., b_a,  0. ],
                               [ 0., b_a, c_a ] ])
    if lattice=='Orthorhombic body-centered':
        a    = cell[0]
        b_a  = cell[1]
        c_a  = cell[2]
        vecs = a/2.*np.array([ [ 1., b_a, c_a ],
                               [-1., b_a, c_a ],
                               [-1.,-b_a, c_a ] ])
    if lattice=='Monoclinic P, unique axis c':
        a    = cell[0]
        b_a  = cell[1]
        c_a  = cell[2]
        c_g  = cell[3]
        s_g  =+sqrt(1.-c_g*c_g)
        vecs = a*np.array([ [      1.,      0.,  0. ],
                            [ b_a*c_g, b_a*s_g,  0. ],
                            [      0.,      0., c_a ] ])
    if lattice=='Monoclinic base-centered':
        a    = cell[0]
        b_a  = cell[1]
        c_a  = cell[2]
        c_g  = cell[3]
        s_g  =+sqrt(1.-c_g*c_g)
        vecs = a*np.array([ [   1./2.,      0.,-c_a/2. ],
                            [ b_a*c_g, b_a*s_g,     0. ],
                            [   1./2.,      0., c_a/2. ] ])
    if lattice=='Triclinic':
        a    = cell[0]
        b_a  = cell[1]
        c_a  = cell[2]
        c_a  = cell[3] #cos(alpha)
        c_b  = cell[4] #cos(beta) 
        c_g  = cell[5] #cos(gamma)
        s_g  =+sqrt(1.-c_g*c_g)
        V2   = a**6.*b_a**2.*c_a**2.*(1.+2.*c_a*c_b*c_g-c_a*c_a-c_b*c_b-c_g*c_g)
        vecs = a*np.array([ [      1.,      0., 0. ],
                            [ b_a*c_g, b_a*s_g, 0. ],
                            [ c_a*c_b, c_a*(c_a-c_b*c_g)/s_g, c_a*sqrt(V2)/(a*b_a*a*c_a*a)/s_g ] ])  
    #
    try: vecs
    except NameError: print("This lattice type is not implemented ):")
    else: return vecs

def crys_to_car(lattice_vectors,reciprocal_space=False):
    """ Calculation of the metric tensor in real or reciprocal space
    
        Input:
            - Lattice vectors (real space)
    """
    
    #Metric tensor
    G = np.array([ [np.dot(a1,a2) for a1 in lattice_vectors] for a2 in lattice_vectors ])
    #Metric tensor in reciprocal space
    if reciprocal_space: return np.linalg.inv(G)
    else: return G
