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
import operator

"""
Using ypp, you can study the convergence of BSE calculations in 2 ways:
  Create a .png of all absorption spectra relevant to the variable you study
  Look at the eigenvalues of the first n "bright" excitons (given a threshold intensity)

The script reads from <folder> all results from <variable> calculations for processing.
The resulting pictures and data files are saved in the ./analyse_bse/ folder.

By default, the graphical interface is deactivated (assuming you run on a cluster because of ypp calls).
See line 2 inside the script.
"""


parser = argparse.ArgumentParser(description='Study convergence on BS calculations using ypp calls.')
parser.add_argument('folder',           help='Folder containing SAVE and convergence runs.')
parser.add_argument('variable',         help='Variable tested (e.g. FFTGvecs)'             )
parser.add_argument('-ne','--numbexc',  help='Number of excitons to read beyond threshold', default=2,type=int)
parser.add_argument('-ie','--intexc',   help='Minimum intensity for excitons to be considered bright', default=0.05,type=float)
parser.add_argument('-de','--degenexc', help='Energy threshold under which different peaks are merged (eV)', default=0.01,type=float)
parser.add_argument('-me','--maxexc',   help='Energy threshold after which excitons are not read anymore (eV)', default=4.0,type=float)
parser.add_argument('-np','--nopack',   help='Skips packing o- files into .json files', action='store_false')
parser.add_argument('-nt','--notext',   help='Skips writing the .dat file', action='store_false')
parser.add_argument('-nd','--nodraw',   help='Skips drawing (plotting) the abs spectra', action='store_false')

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

folder    = args.folder
var       = args.variable
exc_n     = args.numbexc
exc_int   = args.intexc
exc_degen = args.degenexc
exc_max_E = args.maxexc
nopack    = args.nopack

#####################
#  INPUT AND FILES  #
#####################

# Packing results (o-* files) from the calculations into yambopy-friendly .json files
if nopack: # True by default, False if -np used
    print 'Packing ...'
    pack_files_in_folder(folder,mask=var)
    pack_files_in_folder(folder,mask='reference')
    print 'Packing done.'
else:
    print 'Packing skipped.'

# importing data from .json files in <folder>
print 'Importing...'
data = YamboAnalyser(folder)

# extract data according to relevant var
invars = data.get_inputfiles_tag(var)

# Get only files related to the convergence study of the variable,
# ordered to have a smooth plot
keys=[]
sorted_invars = sorted(invars.items(), key=operator.itemgetter(1))

for i in range(0,len(sorted_invars)):
    key=sorted_invars[i][0]
    if key.startswith(var) or key=='reference.json':
        keys.append(key)
print 'Files detected: ',keys

# unit of the input value
unit = invars[keys[0]]['variables'][var][1]

######################

# Output-file filename
os.system('mkdir -p analyse_bse')
outname = './analyse_bse/'+folder+'_'+var

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
    # Get excitons values (runs ypp once)
    a.get_excitons(min_intensity=exc_int,max_energy=exc_max_E,Degen_Step=exc_degen)
    # Write .json file with spectra and eigenenergies
    a.write_json(filename=outname)

    ### Loading data from .json file
    f = open(outname+'.json')
    data = json.load(f)
    f.close()
    print 'JSON file prepared and loaded.'

    ### Plotting the absorption spectra
    # BSE spectra
    plt.plot(data['E/ev[1]'], data['EPS-Im[2]'],label=jobname,lw=2)
#   # Axes : lines for exciton energies (disabled, would make a mess)
#   for n,exciton in enumerate(data['excitons']):
#       plt.axvline(exciton['energy'])

    ### Creating array with exciton values (according to settings)
    l = [inp]
    for n,exciton in enumerate(data['excitons']):
        if n <= exc_n-1:
            l.append(exciton['energy'])

    excitons.append(l)

if args.notext:
    header = 'Columns : '+var+' (in '+unit+') and "bright" excitons eigenenergies in order.'
    print excitons
    np.savetxt(outname+'.dat',excitons,header=header)
    #np.savetxt(outname,excitons,header=header,fmt='%1f')
    print outname+'.dat'
else:
    print '-nt flag : no text produced.'

if args.nodraw:
    plt.xlabel('$\omega$ (eV)')
    plt.gca().yaxis.set_major_locator(plt.NullLocator())
    plt.legend()
    #plt.draw()
    #plt.show()
    plt.savefig(outname+'.png', bbox_inches='tight')
    print outname+'.png'
else:
    print '-nd flag : no plot produced.'

print 'Done.'
