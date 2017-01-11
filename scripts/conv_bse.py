import matplotlib
matplotlib.use('Agg') # prevents crashes if no X server present
from yambopy import *
from qepy import *
import json
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt
import sys
import argparse

""" 
Using ypp, you can:
  Create a .png of all absorption spectra relevant to the variable you study
  Look at the eigenvalues of the first n bright excitons.

Pass the folder and variable name using -f and -v.
Inside the script are additional settings.
"""


parser = argparse.ArgumentParser(description='Helps studying convergence on BS calculations using ypp calls.')
parser.add_argument('-f' ,'--folder'    , help='Folder containing SAVE and convergence runs.')
parser.add_argument('-v' ,'--variable'  , help='Variable tested (e.g. FFTGvecs)')
# TODO have options to disable one or the other (text/graph)
args = parser.parse_args()

folder = args.folder
var    = args.var


###################
# USER PARAMETERS #
###################
# Note : Agg backend (line 2) works great for clusters (no X.org server), 
# but if X is available, you might want to comment it and 
# uncomment the plt.show() at the bottom of the script

# Number of excitons to get
exc_n = 2
# Exciton brightness threshhold
exc_int = 0.005 
# Exciton degenerescence threshold
exc_degen = 0.01 # eV
# Max exciton energies before ignoring
exc_max_E = 5 # eV

#####################
#  INPUT AND FILES  #
#####################

print 'Packing ...'
pack_files_in_folder(folder)
print 'Packing done.'

# importing data from .json files in <folder>
print 'Importing...'
data = YamboAnalyser(folder)

# extract data according to relevant var
invars = data.get_inputfiles_tag(var)

# Get only files related to the convergence study of the var
keys=[]
for key in invars:
    if key.startswith(var):
        keys.append(key)

keys=sorted(keys)
print 'Files detected: ',keys

# unit of the input value
unit = invars[keys[0]]['variables'][var][1]

######################

# Output-file filename
outname = folder+'_'+var+'.dat'

# Array that will contain the output
excitons = []

# Loop over all calculations
for key in keys:
    jobname=key.replace('.json','')
    print jobname

        # input value
        # BndsRn__ is a special case
        if var.startswith('BndsRnX'):
        # format : [1, nband, ...]
            inp = invars[key]['variables'][var][0][1]
        else:
            inp = invars[key]['variables'][var][0]

    print 'Preparing JSON file. Calling ypp if necessary.'
    ### Creating the 'absorptionspectra.json' file
    # It will contain the exciton energies
    y = YamboOut(folder=folder,save_folder=folder)
    # Args : name of job, SAVE folder path, folder where job was run path
    a = YamboBSEAbsorptionSpectra(jobname,path=folder)
    # Get excitons values (runs ypp)
    a.get_excitons(min_intensity=exc_int,max_energy=exc_max_E,Degen_Step=exc_degen)
    # Get excitonic WFs (reads eigenvalues)
    a.get_wavefunctions(Degen_Step=exc_degen,repx=range(-1,2),repy=range(-1,2),repz=range(1))
    # Write .json file
    a.write_json(filename=folder+'_'+jobname)

    ### Loading data from .json file
    f = open(folder+'_'+key)
    data = json.load(f)
    f.close()
    print 'JSON file prepared and loaded.'

    ### Plotting the absorption spectra
    # BSE spectra
    plt.plot(data['E/ev[1]'], data['EPS-Im[2]'],label=jobname,lw=2)
#   # Axes : lines for exciton energies, labels
#   for n,exciton in enumerate(data['excitons']):
#       plt.axvline(exciton['energy'])

    ### Creating array with exciton values (according to settings)
    l = [inp]
    for n,exciton in enumerate(data['excitons']):
        if n <= exc_n-1:
            l.append(exciton['energy'])

    excitons.append(l)

header = 'Variable: '+var+', unit: '+unit
np.savetxt(outname,excitons,header=header,fmt='%1f')
print outname

plt.xlabel('$\omega$ (eV)')
plt.gca().yaxis.set_major_locator(plt.NullLocator())
plt.legend()
#plt.draw()
#plt.show()
plotname = folder+'_'+var+'_abs.png'
plt.savefig(plotname, bbox_inches='tight')
print plotname
print 'Done.'
