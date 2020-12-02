import yaml
import argparse
import numpy as np
import xml.etree.ElementTree as ET

"""
Script to produce band structure data and visualization from QE.

It reads data-file-schema.xml found in the save folders of newer pw versions.

Arguments are:
    -bi, --band_input  -> Path to yaml input file

- INPUT 
This is a yaml file in the following format, to be copy-pasted and edited:

 :: yaml

    ---
     save_dir: "PREFIX.save"
     KPTs: [[0.,0.,0.],[0.,0.5,0.],[0.5,0.5,0.],[0.5,0.,0.]]
     KPTs_labels: [G, Y, S, X]
     shift_Delta_c_v: [0.,1.,1.]

 :: 

The input parameters are:
 - save_dir: path to folder containing 'data-file-schema.xml' (in principle the QE save folder)
 - KPTs: band circuit in reduced coordinates
 - KPTs_labels: labels for the band circuit points
 - shift_Delta_c_v: k-dependent scissor shift as a list of three values (gap shift, cond. stretch, val. stretch)

 - n_valence_bands: AUTOMATICALLY DETERMINE
"""
def get_xml_attrib(xml_root,tag,attrib,repeated=False):
    """ Extract xml attribute
    """
    attribute =  None
    for element in xml_root.iter():
        if element.tag == tag: 
            attribute = element.attrib[attrib]
            if repeated: break

    if attribute is None: raise ValueError('xml attribute %s not found'%attrib)
    else: return attribute

def read_input(inp_file):
    """
    Read input file in yaml format (see above docstring)
    """
    
    # Get input data from yaml file in dictionary form
    stream = open(inp_file, 'r')
    dictionary = yaml.load(stream)
    stream.close()
    
    # Transform input data in shape used by the code
    if 'save_dir' not in dictionary: 
        raise ValueError('QE save path save_dir not found in the input file')
    else: 
        save_dir = dictionary['save_dir']
        prefix   = save_dir.split('/')[-1][:-5]
    
    if 'KPTs' not in dictionary: raise ValueError('Band circuit KPTs not found in the input file.')
    else: KPTs = dictionary['KPTs']
    
    if 'KPTs_labels' not in dictionary: raise ValueError('Kpoint labels KPTs_labels not found in the input file.')
    else: KPTs_labels = dictionary['KPTs_labels']
    
    if 'shift_Delta_c_v' in dictionary: 
        shift_Delta_c_v = np.array( dictionary['shift_Delta_c_v'] ) 
        if ( shift_Delta_c_v == np.array([0.,1.,1.]) ).all(): shift_Delta_c_v = None
    else: 
        shift_Delta_c_v = None
    
    return save_dir, prefix, KPTs, KPTs_labels, shift_Delta_c_v

def setup_BZ_points(input_params):
    """
    Perform preliminary operations such as transform in cartesian
    coordinates and find kpath lengths.
    """
    save_dir, prefix, KPTs, KPTs_labels, shift_Delta_c_v = input_params
    
    # Input files from QE
    eig_xml   = 'eigenval.xml' 
    datafile  = 'data-file-schema.xml'
    datafile_xml = ET.parse( "%s/%s"%(save_dir, datafile)).getroot()

    # Assign path features
    nKPTs = len(KPTs)
    print("=== PATH ===")
    print("directions: %d"%(nKPTs-1))
    if KPTs[-1]==KPTs[0]: print("(Closed path)")
    else: print("(Open path)")

    # Get lattice type 
    ibrav = get_xml_attrib(datafile_xml,'atomic_structure','bravais_index',repeated=True)
    ibrav = int(ibrav)
    lat_type = lattice_dictionary(ibrav)
    exit()
    ## CONTINUE FROM HERE!
    if datafile_xml.find("CELL/CELL_DIMENSIONS") is None: 
        raise ValueError('CELL_PARAMETERS not found')
    else: 
        cell_param = np.array( datafile_xml.find("CELL/CELL_DIMENSIONS").text.split() ).astype(float)

    # Transform edge kpoints from crystal to cartesian coordinates
    G_r = crys_to_car(KPTs,lat_type,cell_param)

    # Length of path lines
    KPTs = np.array(KPTs)
    KPTdist = np.array([ KPTs[i+1]-KPTs[i] for i in range(nKPTs-1) ])

    # ||x||^2 = x^T G_r x
    KPT_lengths = np.array([ sqrt(np.dot(vec, np.dot(G_r, vec))) for vec in KPTdist ])
    KPT_lengths = KPT_lengths/np.max(KPT_lengths)
    print("ratios: ")
    print(KPT_lengths)

def crys_to_car(KPTs,lattice,cell):
    """ Calculation of the metric tensors
    """
    
    lat_vecs = lattice_type(lattice,cell)
    #Metric tensor
    G = np.array([ [np.dot(a1,a2) for a1 in lat_vecs] for a2 in lat_vecs ])
    #Metric tensor in reciprocal space
    G_r = np.linalg.inv(G)
    return G_r

def lattice_dictionary(ibrav):
    """
    Dictionary with text descriptions for ibrav numbers
    """    
    bravais_lattices = {
    0 : 'free',
    1 : 'cubic P (sc)',
    2 : 'cubic F (fcc)',
    3 : 'cubic I (bcc)',
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
    
def lattice_type(lattice,cell):
    """
    This function contains the Bravais lattice basis vector coordinates (ibrav in QE)
    """
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

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Generate band plot')
    parser.add_argument('-bi','--band_input', type=str,help='<Required> Path to input file', required=True)
    args = parser.parse_args()

    inp = args.band_input 
    
    # Read input file
    input_params = read_input(inp)
    
    # BZ setup 
    setup_BZ_points(input_params)
    
    
    


