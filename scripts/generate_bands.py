import yaml
import argparse
import numpy as np

"""
Script to produce band structure data and visualization from QE.

It reads data-file-schema.xml found in the save folders of newer pw versions.

- INPUT 
This is a yaml file in the following format, to be copy-pasted and edited:

 :: yaml

 PLACE EXAMPLE

 :: 

The input parameters are:
 - save_dir: path to folder containing 'data-file-schema.xml' (in principle the QE save folder)
 - KPTs: band circuit in reduced coordinates
 - KPTs_labels: labels for the band circuit points
 - shift_Delta_c_v: k-dependent scissor shift as a list of three values (gap shift, cond. stretch, val. stretch)

 - n_valence_bands: AUTOMATICALLY DETERMINE
"""

def read_input(inp_file):
    """
    Read input file in yaml format (see above docstring)
    """
    
    # Get input data from yaml file in dictionary form
    stream = open(inp_file, 'r')
    dictionary = yaml.load(stream)
    stream.close()
    
    global save_dir
    global prefix
    global KPTs
    global KPTs_labels
    global shift_Delta_c_v
    
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

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Generate band plot')
    parser.add_argument('-i','--input', type=str,help='<Required> Path to input file', required=True)
    args = parser.parse_args()

    inp = args.input 
    
    # Read input file
    read_input(inp)
    
    


