# Version : 6th dec
""" python plot_BSE.py <PathToFolder> <Variable> 
This script allows to track the energy of the n-first bright exciton (settings in script).
Another way to converge BS calculations in parallel to absorption spectra plotting.

The Folder normally contains the SAVE folder as well as <variable> folders. Path should be at least '.'.
Yambo must be able to run, because ypp is called multiple times.
"""
# Future dev : have both excitonic energies and abs spectra.

import matplotlib
matplotlib.use('Agg') # prevents crashes if no X server present
from yambopy import *
from qepy import *
import json
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt
import sys

path = sys.argv[1]
var = sys.argv[2]

###################
# USER PARAMETERS #
###################

# Number of excitons to get
exc_n = 2
# Exciton brightness threshhold
exc_int = 0.005
# Exciton degenerescence threshold
exc_degen = 0.01
# Max exciton energies before ignoring
exc_max_E = 5 # eV

#####################
#  INPUT AND FILES  #
#####################

print 'Packing ...'
pack_files_in_folder(path)
print 'Packing done.'

# importing data from .json files in <folder>
print 'Importing...'
data = YamboAnalyser(path)

# extract data according to relevant var
invars = data.get_inputfiles_tag(var)

# Get only files related to the convergence study of the var
keys=[]
for key in invars:
    if key.startswith(var):
        keys.append(key)

keys=sorted(keys)
print keys

# unit of the input value
unit = invars[keys[0]]['variables'][var][1]

######################

# Output-file filename
outname = path+'_'+var+'.dat'

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

	print 'Preparing JSON file. Calling ypp ...'
	### Creating the 'absorptionspectra.json' file
	# It will contain the exciton energies
	y = YamboOut(folder=path,save_folder=path)
	# Args : name of job, SAVE folder path, folder where job was run path
	a = YamboBSEAbsorptionSpectra(jobname,path=path)
	# Get excitons values
	a.get_excitons(min_intensity=exc_int,max_energy=exc_max_E,Degen_Step=exc_degen)
	# Get excitonic WFs (contains eigenvalues)
	a.get_wavefunctions(Degen_Step=exc_degen,repx=range(-1,2),repy=range(-1,2),repz=range(1))
	# Write .json file
	a.write_json(filename=path+'_'+jobname)

	### Loading data from .json file
	f = open(path+'_'+jobname+'.json')
	data = json.load(f)
	f.close()
	print 'JSON file prepared and loaded.'
#
#	### Plotting the absorption spectra
#	print 'Absorption spectra ...'
#	# BSE spectra
#	plt.plot(data['E/ev[1]'], data['EPS-Im[2]'],label='BSE',lw=4)
#	# Axes : lines for exciton energies, labels
#	for n,exciton in enumerate(data['excitons']):
#	    plt.axvline(exciton['energy'])
#	plt.xlabel('$\omega$ (eV)')
#	plt.gca().yaxis.set_major_locator(plt.NullLocator())
#	plt.legend()
#	#plt.draw()
#	#plt.show()
#	plt.savefig(path+'_'+jobname+'_abs.png', bbox_inches='tight')

	l = [inp]
	for n,exciton in enumerate(data['excitons']):
		if n <= exc_n-1:
			print exciton['energy']
			l.append(exciton['energy'])

	excitons.append(l)

print excitons
header = 'Variable: '+var+', unit: '+unit
np.savetxt(outname,excitons,header=header,fmt='%1f')
print outname
