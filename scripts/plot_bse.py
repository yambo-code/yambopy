# Version : nov. 18th, 11:30
""" python plot_BSE.py <JobName> <PathToJobFolder>
The Job Folder normally contains the SAVE folder as well as <JobName> folder. Path should be at least '.'.
Yambo must be able to run, because ypp is called multiple times.
"""

import matplotlib
matplotlib.use('Agg') # prevents crashes if no X server present
from yambopy import *
from qepy import *
import json
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt
import sys

jobname = sys.argv[1]
path = sys.argv[2]

print 'Preparing JSON file. Calling ypp ...'
### Creating the 'absorptionspectra.json' file
y = YamboOut(folder=path,save_folder=path)
# Args : name of job, SAVE folder path, folder where job was run path
a = YamboBSEAbsorptionSpectra(jobname,path=path)
# Get excitons values
a.get_excitons(min_intensity=0.0005,max_energy=6,Degen_Step=0.01)
# Get excitonic WFs
a.get_wavefunctions(Degen_Step=0.01,repx=range(-1,2),repy=range(-1,2),repz=range(1))
# Write .json file
a.write_json(filename=path+'_'+jobname)

### Loading data from .json file
f = open(path+'_'+jobname+'.json')
data = json.load(f)
f.close()
print 'JSON file prepared and loaded.'

### Plotting the absorption spectra
print 'Absorption spectra ...'
# BSE spectra
plt.plot(data['E/ev[1]'], data['EPS-Im[2]'],label='BSE',lw=4)
# IP spectra
plt.plot(data['E/ev[1]'], data['EPSo-Im[4]'],label='IP',lw=4)
# Axes : lines for exciton energies, labels
for n,exciton in enumerate(data['excitons']):
    plt.axvline(exciton['energy'])
plt.xlabel('$\omega$ (eV)')
plt.gca().yaxis.set_major_locator(plt.NullLocator())
plt.legend()
#plt.draw()
#plt.show()
plt.savefig(path+'_'+jobname+'_abs.png', bbox_inches='tight')

### Lattice : necessary to determine bounds of exciton weights plot
print 'Reciprocal lattice ...'
lat = np.array(data['lattice'])
rlat = rec_lat(lat)
x,y,z = np.array(rlat )
cut=.65
xmin,ymin,_ = -(x+y)*cut
xmax,ymax,_ = +(x+y)*cut


### Plot excitons
print 'Excitons weights ...'
fig = plt.figure(figsize=(16,8))
# Issue if dvipng is not installed
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif',serif="Computer Modern Roman",size=20)

sorted_excitons = sorted(data['excitons'],key=lambda x: x['energy'])

#colormap
cmap = plt.get_cmap("viridis")

for n,exciton in enumerate(sorted_excitons[:8]):
    ax = plt.subplot(2,4,n+1)
    ax.set_aspect('equal', 'datalim')
    w = np.array(exciton['weights'])
    qpt = np.array(exciton['qpts'])
    ax.scatter(qpt[:,0], qpt[:,1], marker='H', s=2.5, color=[cmap(sqrt(c)) for c in w], label="e: %lf"%exciton['energy'])
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    plt.gca().yaxis.set_major_locator(plt.NullLocator())
    plt.gca().xaxis.set_major_locator(plt.NullLocator())
    ax.legend(prop={'size':12})

#plt.draw()
#plt.show()
plt.savefig(path+'_'+jobname+'_ex.png', bbox_inches='tight')

print 'Done.'
