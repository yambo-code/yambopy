from netCDF4 import *
import argparse
import numpy as np

"""
Script to convert RL number in Ry energy units using ndb.gops.

Inputs:
 1. --ndb_gops='path/to/folder/with/ndb.gops' [i.e., SAVE]
 2. --value= value to convert with units, e.g., '11 RL' or '5 Ry']

The script will read ndb.gops and find the nearest completed G-shell, then give the
converted value in Ry (RL) to the one supplied in input.

"""

def convert(value,ndb_gops):
    """ Perform the conversion
    """ 

    def get_closed_shell(value,shells):
        """ Returns index of closest completed shell to input value
        """
        ind = np.searchsorted(shells, value, side='left')
        if ind==len(shells): return None
        else:                return ind

    Ry2Ha=0.5
    ndb_gops = ndb_gops+"/ndb.gops"
    try:    value, unit = value
    except: raise ValueError("[ERROR] Incorrect format of value to be converted")

    # Read shells from ndb.gops
    db = Dataset(ndb_gops,'r')
    ng_in_shell = np.array(db['ng_in_shell'][:]).astype(int)
    E_of_shell  = db['E_of_shell'][:]
    db.close()

    if unit=='RL': 
        value=int(value)
        shells_i = ng_in_shell
        shells_o = E_of_shell
        unit_o   = 'Ry'
    elif unit=='Ry': 
        value=float(value)*Ry2Ha
        shells_i = E_of_shell
        shells_o = ng_in_shell
        unit_o   = 'RL'
    else: raise ValueError("[ERROR] Unit %s currently not supported"%unit)

    ind = get_closed_shell(value,shells_i)
    if ind is not None:
        close_shell_value = shells_i[ind]
        converted_value   = shells_o[ind]
        print(converted_value,' ',unit_o,' (closest closed shell)')
    else: raise ValueError("[ERROR] value supplied is bigger than maximum shell at %s %s"%(str(shell_i[-1]),unit))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert RL number in energy units (Ry)')
    parser.add_argument('-gops','--ndb_gops', type=str, default="./SAVE", help='<Optional> Path to folder with ndb.gops (default: ./SAVE)')
    parser.add_argument('-v','--value', type=str,help="<Required> Value to be converted along with units, e.g.: '11 RL' or '5 Ry'",nargs=2,required=True)
    args = parser.parse_args()

    ndb_gops = args.ndb_gops
    value    = args.value
    
    convert(value,ndb_gops)
