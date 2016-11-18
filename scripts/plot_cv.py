import matplotlib
matplotlib.use('Agg') # prevents crashes if no X server present
from yambopy import *
import matplotlib.pyplot as plt
import sys

""" Input : python data_grab.py <folder> <variable>
    
    Uses .json files from pack_files_in_folder routine
    to help plotting convergence for the <variable> passed as argument
    <folder> is the folder containing all the different run folders,
    wherein the different o-* files are located.
    NB : The reference run is not used and should be disabled if using the 
    yambopy.optimize() routine using arg "ref_run=False".
    """

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

# Writing the unit of the input value in the first line
unit = invars[keys[0]]['variables'][var][1]
f.write('# Unit of the input value: '+str(unit)+'\n')

# In the 'short' Yambo output format :
# K-coord | Band | E_LDA | E_GW-E_LDA | Sc|Eo
for key in keys:
    # input value
    # GbndRnge and BndsRnX_ are special cases
    if var.startswith('GbndRng') or var.startswith('BndsRnX'):
        # format : [1, nband, ...]
        inp = invars[key]['variables'][var][0][1]
    else:
        inp = invars[key]['variables'][var][0]

    # output value : E_GW between bands 26 and 27
    # in my particular case, I have band 26 on line 2 and band 27 on line 3 (index 1 and 2)
    # TODO : check band using if [input_key][i][1] == 26|27
    # col 2+3 is E_GW ; bands 27-26 is gap
    out = outvars[key][2][2]+outvars[key][2][3]-( outvars[key][1][2]+outvars[key][1][3] )
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
