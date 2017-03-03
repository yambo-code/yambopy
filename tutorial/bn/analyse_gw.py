import matplotlib
#matplotlib.use('Agg') # prevents crashes if no X server present (clusters)
from yambopy import *
import matplotlib.pyplot as plt
import sys
import argparse
import numpy as np

"""
Study the convergence of GW calculations by looking at the change in band-gap value.

The script reads from <folder> all results from <variable> calculations (skipping the reference run)
and display them. To avoid running the reference run for nothing, use optimize(...,ref_run=False).
Please note that the first value in the convergence dictionnary will thus not be run.

Use the band and k-point options (or change default values) according to the size of your k-grid and
the location of the band extrema.
"""

parser = argparse.ArgumentParser(description='Study GW convergence with regards to the band-gap value.')
#parser.add_argument('-f' ,'--folder'    , help='Folder containing SAVE and convergence runs.', required=True)
#parser.add_argument('-v' ,'--variable'  , help='Variable tested (e.g. FFTGvecs)'             , required=True)
parser.add_argument('folder'    , help='Folder containing SAVE and convergence runs.')
parser.add_argument('variable'  , help='Variable tested (e.g. FFTGvecs)'             )
parser.add_argument('-bc','--bandc'     , help='Lowest conduction band number'    , default=27,type=int)
parser.add_argument('-kc','--kpointc'   , help='K-point index for conduction band', default=19,type=int)
parser.add_argument('-bv','--bandv'     , help='Highest valence band number'      , default=26,type=int)
parser.add_argument('-kv','--kpointv'   , help='K-point index for valence band'   , default=19,type=int)
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

print 'Valence band: ',bandv,'conduction band: ',bandc
print 'K-point VB: ',kpointv, ' k-point CB: ',kpointc


# Packing results (o-* files) from the calculations into yambopy-friendly .json files
if nopack: # True by default, False if -np used
    print 'Packing ...'
    pack_files_in_folder(folder)
    print 'Packing done.'

# importing data from .json files in <folder>
print 'Importing...'
data = YamboAnalyser(folder)

# extract data according to relevant variable
outvars = data.get_data(var)
invars = data.get_inputfiles_tag(var)
tags = data.get_tags(var)

# Get only files related to the convergence study of the variable
keys=[]
for key in invars:
    if key.startswith(var):
        keys.append(key)

# Ordered to help plotting with lines
keys=sorted(keys)

print 'Preparing output...'
### Output

# Unit of the variable :
unit = invars[keys[0]]['variables'][var][1]

# The following variables are used to make the script compatible with both short and extended output
kpindex = tags[keys[0]].tolist().index('K-point')
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
    header = 'Variable: '+var+', unit: '+str(unit)
    np.savetxt(filename,array,delimiter='\t',header=header)
    print filename

plt.plot(array[:,0],array[:,1],'o-')
plt.xlabel(var+' ('+unit+')')
plt.ylabel('$E_{GW} = E_{LDA} + \Delta E$')
plt.show()
#plt.savefig(folder+'_'+var+'.png')
