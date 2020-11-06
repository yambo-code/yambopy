import os
from yambopy import *
from schedulerpy import *
import argparse

"""
Script to produce band structure data and visualization from QE.

It reads data-file-schema.xml found in the save folders of newer pw versions.

- INPUT 
This is a yaml file in the following format, to be copy-pasted and edited:

 :: yaml

 PLACE EXAMPLE

 :: 

The input parameters are:
 - save_dir: path to QE save folder
 - KPTs: band circuit in reduced coordinates
 - KPTs_labels: labels for the band circuit points
 - shift_Delta_c_v: k-dependent scissor shift as a list of three values (gap shift, cond. stretch, val. stretch)

 - n_valence_bands: AUTOMATICALLY DETERMINE
"""

def some_function():
    """
    Sample function
    """
    return some_function.__doc__

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Generate band plot')
    parser.add_argument('-i','--input', type=str,help='<Required> Path to input file', required=True)
    args = parser.parse_args()

    inp = args.input 


