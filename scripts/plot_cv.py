import matplotlib
matplotlib.use('Agg') # prevents crashes if no X server present
from yambopy import *
import matplotlib.pyplot as plt
import sys

""" Input : python data_grab.py <folder> <variable>
    Other parameters must be changed depending on systems in the first lines of the script.
    
    Uses .json files from pack_files_in_folder routine
    to help plotting the energy gap convergence for the <variable> passed as argument
    <folder> is the folder containing all the different run folders,
    wherein the different o-* files are located.
    NB : The reference run is not used and should be disabled if using the 
    yambopy.optimize() routine using arg "ref_run=False".
    """

## USER-DEFINED VARIABLES
# Variables used to find where the gap is.
# bandc | bandv : number of the conduction | valence band
# kpoint : number of the k-point which is to be evaluated

# Exemple values for 1L-MoS2 on a 12x12 grid
bandc = 27 ; bandv = 26
kpoint = 19
print 'Conduction band: ',bandc,'\nValence band: ',bandv,'\nK-point:',kpoint

folder = sys.argv[1]
var = sys.argv[2]

# <folder> would be the folder containing all .json files
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
# arrays for matplotlib plot
# file for later use
inparray = []
outarray = []
filename = folder+'_'+var+'.dat'
f = open(filename,'w')

# The following variables are used to make the script compatible with both short and extended output
kpindex = tags[keys[0]].tolist().index('K-point')
bdindex = tags[keys[0]].tolist().index('Band')
e0index = tags[keys[0]].tolist().index('Eo')
gwindex = tags[keys[0]].tolist().index('E-Eo')

# Writing the unit of the input value in the first line
unit = invars[keys[0]]['variables'][var][1]
f.write('# Unit of the input value: '+str(unit)+'\n')

for key in keys:
    # input value
    # GbndRnge and BndsRnX_ are special cases
    if var.startswith('GbndRng') or var.startswith('BndsRnX'):
        # format : [1, nband, ...]
        inp = invars[key]['variables'][var][0][1]
    else:
        inp = invars[key]['variables'][var][0]

    # Output value (gap energy)
    # First the relevant lines are identified
    valence=[]
    conduction=[]
    for i in range(len(outvars[key]+1)):
        if outvars[key][i][kpindex]==kpoint and outvars[key][i][bdindex]==bandc:
		conduction=outvars[key][i]
        elif outvars[key][i][kpindex]==kpoint and outvars[key][i][bdindex]==bandv:
		valence = outvars[key][i]
    # Then the gap can be calculated
    out = conduction[e0index]+conduction[gwindex]-(valence[e0index]+valence[gwindex]) 

    #writing value and energy diff in file
    s=str(inp)+'\t'+str(out)+'\n'
    f.write(s)
    inparray.append([inp])
    outarray.append([out])

plt.plot(inparray,outarray,'o-')
plt.xlabel(var+' ('+unit+')')
plt.ylabel('E_gw = E_lda + \Delta E')
#plt.show()
plt.savefig(folder+'_'+var+'.png')
print filename
