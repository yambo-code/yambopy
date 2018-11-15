# Copyright (C) 2015 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy

'''
Plot the excitonic weights on the brillouin zone using ypp
Plot the absorption spectrum with the excitonic energies
'''
from __future__ import print_function

from builtins import str
from builtins import range
import matplotlib
matplotlib.use('Agg') # prevents crashes if no X server present
from yambopy import *
from qepy import *
import json
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt, ceil
import argparse

parser = argparse.ArgumentParser(description='Plot excitonic weights in BZ or excitonic eigenenergies on optical spectrum using ypp.')
parser.add_argument('folder', help='Folder with SAVE and run.')
parser.add_argument('jobname', help='Job from which excitons are of interest.')
parser.add_argument('-ie','--intexc' , help='Minimum intensity for excitons to be considered bright', default=0.1)
parser.add_argument('-de','--degenexc', help='Energy threshold under which different peaks are merged (eV)', default=0.01)
parser.add_argument('-me','--maxexc', help='Energy threshold after which excitons are not read anymore (eV)', default=4.0)
parser.add_argument('-c','--cut', help='Set the limits of the BZ plot.', default=0.2)
parser.add_argument('-ny','--noypp', help='Skips the ypp calls and the .json file building.', action='store_false')
parser.add_argument('-nw','--noweight', help='Skips the plot of the weights in the BZ.', action='store_false')
parser.add_argument('-ns','--nospectrum', help='Skips the plot of the spectrum with the excitons eigenenergies.', action='store_false')
args = parser.parse_args()

jobname = args.jobname
path    = args.folder 
cut     = args.cut
exc_int   = args.intexc
exc_degen = args.degenexc
exc_max_E = args.maxexc

def get_var(dictionary,variables):
    """
    To have compatibility with different version of yambo
    We provide a list of different possible tags
    """
    for var in variables:
        if var in dictionary:
            return dictionary[var]
    raise ValueError( 'Could not find the variables %s in the output file'%str(variables) )

if args.noypp:
    print('Preparing JSON file. Calling ypp ...')
    ### Creating the 'absorptionspectra.json' file
    y = YamboOut(folder=path,save_folder=path)
    # Args : name of job, SAVE folder path, folder where job was run path
    a = YamboBSEAbsorptionSpectra(jobname,path=path)
    # Get excitons values
    a.get_excitons(min_intensity=exc_int,max_energy=exc_max_E,Degen_Step=exc_degen)
    # Get excitonic WFs
    a.get_wavefunctions(Degen_Step=exc_degen,repx=list(range(-1,2)),repy=list(range(-1,2)),repz=list(range(1)))
    # Write .json file
    a.write_json(filename=path+'_'+jobname)
else:
    print('YPP call disabled.')

### Loading data from .json file
f = open(path+'_'+jobname+'.json')
data = json.load(f)
f.close()
print('JSON file loaded.')

### Plotting the absorption spectra
print('Absorption spectra ...')
# BSE spectra
plt.plot(get_var(data,['E/ev','E/ev[1]']), get_var(data,['EPS-Im[2]' ]),label='BSE',lw=4)
# IP spectra
plt.plot(get_var(data,['E/ev','E/ev[1]']), get_var(data,['EPSo-Im[4]']),label='IP',lw=4)
# Axes : lines for exciton energies, labels
nexcitons = len(data['excitons'])
print('Number of excitons: ', nexcitons) 
for n,exciton in enumerate(data['excitons']):
    plt.axvline(exciton['energy'])
plt.xlabel('$\omega$ (eV)')
plt.ylabel('Intensity (arb. units)')
plt.gca().yaxis.set_major_locator(plt.NullLocator())
plt.legend()
#plt.draw()
#plt.show()
plt.savefig(path+'_'+jobname+'_abs.png', bbox_inches='tight')
print(path+'_'+jobname+'_abs.png')

### Lattice : necessary to determine bounds of exciton weights plot
print('Reciprocal lattice ...')
lat = np.array(data['lattice'])
rlat = rec_lat(lat)
x,y,z = np.array(rlat )
xmin,ymin,_ = -(x+y)*cut
xmax,ymax,_ = +(x+y)*cut


### Plot excitons
print('Excitons weights ...')
nx = int(ceil(sqrt(nexcitons)))
ny = int(ceil(nexcitons*1.0/nx))
print("cols:",nx)
print("rows:",ny)
cmap = plt.get_cmap("gist_heat_r")

fig = plt.figure(figsize=(nx*4,ny*4))
sorted_excitons = sorted(data['excitons'],key=lambda x: x['energy'])

for n,exciton in enumerate(sorted_excitons):
    #get data
    w   = np.array(exciton['weights'])
    qpt = np.array(exciton['qpts'])

    #plot
    ax = plt.subplot(ny,nx,n+1)
    ax.scatter(qpt[:,0], qpt[:,1], s=20, c=w, marker='H', cmap=cmap, lw=0, label="e: %lf"%exciton['energy'])

    # axis
    plt.xlim([-cut,cut])
    plt.ylim([-cut,cut])
    plt.gca().yaxis.set_major_locator(plt.NullLocator())
    plt.gca().xaxis.set_major_locator(plt.NullLocator())
    ax.set_aspect('equal')

#plt.draw()
#plt.show()
plt.savefig(path+'_'+jobname+'_ex.png', bbox_inches='tight')
print('path'+'_'+jobname+'_ex.png')
print('Done.')
