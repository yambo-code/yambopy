from __future__ import print_function
from builtins import str
from builtins import range
import matplotlib
#matplotlib.use('Agg') # prevents crashes if no X server present (clusters)
from yambopy import *
import matplotlib.pyplot as plt
import sys
import argparse
import numpy as np
import operator

"""
Study the convergence of GW calculations by looking at the change in band-gap value.

The script reads from <folder> all results from <variable> calculations and display them.

Use the band and k-point options (or change default values) according to the size of your k-grid and
the location of the band extrema.
"""

parser = argparse.ArgumentParser(description='Study GW convergence with regards to the band-gap value.')
parser.add_argument('folder'    , help='Folder containing SAVE and convergence runs.')
parser.add_argument('variable'  , help='Variable tested (e.g. FFTGvecs)'             )
parser.add_argument('-bc','--bandc'     , help='Lowest conduction band number'    , default=53, type=int)
parser.add_argument('-kc','--kpointc'   , help='K-point index for conduction band', default=19, type=int)
parser.add_argument('-bv','--bandv'     , help='Highest valence band number'      , default=52, type=int)
parser.add_argument('-kv','--kpointv'   , help='K-point index for valence band'   , default=1, type=int)
parser.add_argument('-np','--nopack'    , help='Skips packing o- files into .json files', action='store_false')
parser.add_argument('-t'  ,'--text'      , help='Also print a text file for reference'   , action='store_true')
args = parser.parse_args()

folder = args.folder
var    = args.variable
bandc  = args.bandc
kpointc= args.kpointc
bandv  = args.bandv
kpointv= args.kpointv
nopack = args.nopack
text   = args.text

print('Valence band: ',bandv,'conduction band: ',bandc)
print('K-point VB: ',kpointv, ' k-point CB: ',kpointc)


# Packing results (o-* files) from the calculations into yambopy-friendly .json files
if nopack: # True by default, False if -np used
    print('Packing ...')
    pack_files_in_folder(folder,mask=var)
    pack_files_in_folder(folder,mask='reference')
    print('Packing done.')
else:
    print('Packing skipped.')

# importing data from .json files in <folder>
print('Importing...')
data = YamboAnalyser(folder)

# extract data according to relevant variable
outvars = data.get_data(var)
invars = data.get_inputfiles_tag(var)
tags = data.get_tags(var)

# Get only files related to the convergence study of the variable,
# ordered to have a smooth plot
keys=[]
sorted_invars = sorted(list(invars.items()), key=operator.itemgetter(1))

for i in range(0,len(sorted_invars)):
    key=sorted_invars[i][0]
    if key.startswith(var) or key=='reference.json':
        keys.append(key)
print('Files detected: ',keys)

print('Preparing output...')
### Output

# Unit of the variable :
unit = invars[keys[0]]['variables'][var][1]

# The following variables are used to make the script compatible with both short and extended output
#kpindex = tags[keys[0]].tolist().index('K-point')
kpindex = tags[keys[0]].tolist().index('Kpoint_index')  # Alejandro
bdindex = tags[keys[0]].tolist().index('Band')
e0index = tags[keys[0]].tolist().index('Eo')
gwindex = tags[keys[0]].tolist().index('E-Eo')


array = np.zeros((len(keys),2))

for i,key in enumerate(keys):
    # input value
    # GbndRnge and BndsRnX_ are special cases
    if var.startswith('GbndRng') or var.startswith('BndsRnX'):
        # format : [1, nband, ...]
        array[i][0] = invars[key]['variables'][var][0][1]
    else:
        array[i][0] = invars[key]['variables'][var][0]

    # Output value (gap energy)
    # First the relevant lines are identified
    valence=[]
    conduction=[]
    for j in range(len(outvars[key]+1)):
        if outvars[key][j][kpindex]==kpointc and outvars[key][j][bdindex]==bandc:
                conduction=outvars[key][j]
        elif outvars[key][j][kpindex]==kpointv and outvars[key][j][bdindex]==bandv:
                valence = outvars[key][j]
    # Then the gap can be calculated
    array[i][1] = conduction[e0index]+conduction[gwindex]-(valence[e0index]+valence[gwindex])

if text:
    filename = folder+'_'+var+'.dat'
    header = +var+'('+str(unit)+'), gap'
    np.savetxt(filename,array,delimiter='\t',header=header)
    print(filename)

plt.plot(array[:,0],array[:,1],'o-')
plt.xlabel(var+' ('+unit+')')
plt.ylabel('E_gw = E_lda + \Delta E')
plt.show()
#plt.savefig(folder+'_'+var+'.png')

# Plot all of the different GW bandstructures in the same plot
#ya = YamboAnalyser(folder)
#ya.plot_gw('qp',cols=(lambda x: x[2],lambda x: x[3]+x[4]))

