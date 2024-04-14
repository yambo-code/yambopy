import yaml
import argparse
import numpy as np
from math import sqrt
import xml.etree.ElementTree as ET
from qepy import xml
from qepy import bravais
from copy import deepcopy
from yambocommandline.commands import band_plots

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
"""
def read_input(inp_file):
    """
    Read input file in yaml format (see above docstring)
    """
    
    # Get input data from yaml file in dictionary form
    stream = open(inp_file, 'r')
    dictionary = yaml.safe_load(stream)
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
    
def get_data_from_xml(input_params):
    """
    Read the data-file-schema.xml to obtain all the needed variables
    
    TO BE WRAPPED WITH A DEDICATED PWXML CLASS
    """
    save_dir, prefix, KPTs, KPTs_labels, shift_Delta_c_v = input_params
    
    # Input files from QE
    datafile  = 'data-file-schema.xml'
    datafile_xml = ET.parse( "%s/%s"%(save_dir, datafile))
    
    # Get lattice vectors
    lat_vecs = []
    lat_vecs_tags = ['a1','a2','a3']
    for tag in lat_vecs_tags:
        lat_vecs.append( xml.get_xml_data(datafile_xml,tag,as_type=float,repeated=False) )
    lat_vecs = np.array(lat_vecs)
    
    # Get kpoint number
    nkpoints = xml.get_xml_data(datafile_xml,'nk',as_type=int,repeated=True)
    
    # Get band number
    nbands   = xml.get_xml_data(datafile_xml,'nbnd',as_type=int,repeated=True)
    
    # Get kpoints cartesian coordinates
    kpts_cart, eigen, occs,= xml.get_xml_nk_bands(datafile_xml)
    
    # Check if the system is magnetic
    try:
        spin     = eval( xml.get_xml_data(datafile_xml,'lsda',as_type=str).capitalize() )
        noncolin = eval( xml.get_xml_data(datafile_xml,'noncolin',as_type=str).capitalize() )
        if spin and not noncolin: l_magnetic = True
        else: l_magnetic = False
    except ValueError: l_magnetic = False

    # Get Fermi energy
    try: fermi_e = xml.get_xml_data(datafile_xml,'fermi_energy',as_type=float)
    except ValueError: fermi_e = xml.get_xml_data(datafile_xml,'two_fermi_energies',as_type=float)[0] # majority spin Fermi energy

    # Scissor operator check
    if shift_Delta_c_v is not None:
        occ_kind = xml.get_xml_data(datafile_xml,'occupations_kind',as_type=str)
        if occ_kind!='fixed': 
            print('[WARNING]: if your system is a metal, the scissor shift will NOT work!')
    
    # Get topmost valence band (only valid for gapped systems)    
    if not l_magnetic: n_VBM = int( np.sum(occs[0]) )
    else:              n_VBM = int( np.sum(occs[0,:nbands]) ) # majority spins top valence  

    # Get eigenvalues
    eigen, s_eigen = process_bands(eigen,kpts_cart,fermi_e,shift_Delta_c_v,n_VBM)

    return lat_vecs, nkpoints, nbands, kpts_cart, eigen, s_eigen, n_VBM, l_magnetic
    
def setup_BZ_points(input_params,output_data):
    """
    Perform preliminary operations such as transform in cartesian
    coordinates and find kpath lengths.
    
    ***ASSUMING SAME NUMBER OF POINTS FOR EACH SEGMENT***)
    """
    save_dir, prefix, KPTs, KPTs_labels, shift_Delta_c_v  = input_params
    lat_vecs, nkpoints, nbands, kpts_cart, eigen, s_eigen, n_VBM, l_magnetic = output_data

    # Describe BZ symmetry lines
    nKPTs = len(KPTs)
    nkpt_per_direction = (nkpoints-1)/(nKPTs-1)

    # Metric tensor in reciprocal space
    G_r = bravais.crys_to_car(lat_vecs,reciprocal_space=True)

    # Length of path lines
    KPTs = np.array(KPTs)
    KPTdist = np.array([ KPTs[i+1]-KPTs[i] for i in range(nKPTs-1) ])

    # Transform edge kpoints from crystal to cartesian coordinates:
    #               ||x||^2 = x^T G_r x
    KPT_lengths = np.array([ sqrt(np.dot(vec, np.dot(G_r, vec))) for vec in KPTdist ])
    KPT_lengths = KPT_lengths/np.max(KPT_lengths)
    
    # Generating scaled steps
    kstps = np.zeros(nkpoints)
    points = []
    for i in range(1,nKPTs):
        n_l, n_r = np.array([i-1, i])*nkpt_per_direction
        n_l = int(n_l)
        n_r = int(n_r)
        for ik in range(n_l+1,n_r+1): kstps[ik] = kstps[ik-1]+KPT_lengths[i-1]
        points.append(kstps[n_r])

    # Get energy gaps
    gaps = find_gaps(s_eigen,n_VBM)

    # Print info
    setup_info(nKPTs,nkpt_per_direction,KPTs,KPT_lengths,nkpoints,nbands,points,gaps,shift_Delta_c_v,l_magnetic)
    
    if shift_Delta_c_v is not None: return prefix,nkpoints,nbands,kstps,KPTs_labels,points,gaps,eigen,s_eigen    
    else: return prefix,nkpoints,nbands,kstps,KPTs_labels,points,gaps,eigen
        
    
def setup_info(nKPTs,nkpt_per_direction,KPTs,KPT_lengths,nkpoints,nbands,points,gaps,scissor,l_magnetic):
    """
    Print information about the system
    """
    # BZ info
    print("=== PATH ===")
    print("directions: %d"%(nKPTs-1))
    if np.array_equal(KPTs[-1],KPTs[0]): print("(Closed path)")
    else: print("(Open path)")  
    print("ratios: ")
    print(KPT_lengths) 
    
    # Bands info 
    print("=== BANDS ===")
    print("nkpoints: %d"%nkpoints)
    print("kpoint density per direction: %d"%nkpt_per_direction)
    print("nbands: %d"%nbands)  
    if scissor is None: print("scissor shift: No")
    else: print("scissor shift: Yes")
    if gaps[0] > 1e-6 or gaps[1] >1.e-6:
        print("direct band gap:   %f eV"%gaps[0])
        print("indirect band gap: %f eV"%gaps[1])
    else:
        print("direct band gap:   %f eV"%gaps[0])
        print("indirect band gap: %f eV"%gaps[1])
        print("This is a metal.")
    
    if l_magnetic: print("The system is collinear magnetic.")

    # Symmetry points info
    print("=== PLOT ===")
    print("Internal high-symmetry points at: ")
    print(points)  

def find_gaps(eigen,n_val):
    """
    Find smallest direct and indirect band gaps
    """
    Nk,Nb = eigen.shape
    top_v, bottom_c = eigen[:,n_val-1], eigen[:,n_val]    
    dir_gaps = bottom_c - top_v
    E_dir_gap = np.min(dir_gaps)
    indir_gaps = []
    for ik in range(Nk):
        for ip in range(Nk):
            if ip!=ik: indir_gaps.append(bottom_c[ip]-top_v[ik])
    indir_gaps = np.array(indir_gaps)
    E_ind_gap = np.min(indir_gaps)
    return E_dir_gap, E_ind_gap
            
def process_bands(eigen,kpts_cart,fermi_e,shift_Delta_c_v,n_val):
    """
    Fix Fermi level to zero, apply scissor shift, return useful dictionary for plotting 
    
    Output:
        eigenvalues and shifted eigenvalues (equal if scissor not applied)
    """
    # Rescale in eV and shift Fermi level to zero
    ha2ev  = 27.211396132
    processed_eigen  = ha2ev*(eigen-fermi_e)              

    # Apply scissor shift
    shifted_eigen = deepcopy(processed_eigen)
    if shift_Delta_c_v is not None:
        scissor = shift_Delta_c_v
        top_v, bottom_c = processed_eigen[:,n_val-1], processed_eigen[:,n_val]
        ind_k_dir_gap = np.argmin(bottom_c-top_v)
        ev_max, ec_min = top_v[ind_k_dir_gap], bottom_c[ind_k_dir_gap]
        for ik in range( len(kpts_cart) ):
            for ib in range( len(processed_eigen[0]) ):
                if ib<n_val: shifted_eigen[ik][ib] = ev_max-(ev_max-processed_eigen[ik][ib])*scissor[2]
                else:        shifted_eigen[ik][ib] = ec_min+scissor[0]+(processed_eigen[ik][ib]-ec_min)*scissor[1]

    return processed_eigen, shifted_eigen

def launch_plot(data,plt_type='bands',out_name=None,erange=None,show=True):
    """
    Plot driver. 
    
    Calls band_plots with the appropriate plt_type
    """
    
    band_plots.plot_driver(data,plt_type,out_name=out_name,erange=erange,show=show)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Generate band plot')
    parser.add_argument('-bi','--band_input', type=str,help='<Required> Path to input file', required=True)
    args = parser.parse_args()

    inp = args.band_input 
    
    # Read user input file
    input_params = read_input(inp)
    
    # Read data from .xml file
    output_data = get_data_from_xml(input_params)
    
    # Plot setup
    data_to_plot = setup_BZ_points(input_params,output_data)
    
    # Produce the plot
    launch_plot(data_to_plot,plt_type='bands',out_name=None,erange=4.,show=True)
    
    


