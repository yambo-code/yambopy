# Copyright (c) 2016, Henrique Miranda
# All rights reserved.
#
#
# This file is part of the yambopy project
#
from yambopy import *
import json
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt

#colormap
cmap = plt.get_cmap("gist_heat_r")

#read file
f = open('absorptionspectra.json')
data = json.load(f)
f.close()

print "nexitons", len(data['excitons'])

#plot the absorption spectra
plt.plot(data['E/ev'], data['EPS-Im'],label='BSE',lw=4)
plt.plot(data['E/ev'], data['EPSo-Im'],label='IP',lw=4)
for n,exciton in enumerate(data['excitons']):
    plt.axvline(exciton['energy'])
plt.xlabel('$\omega$ (eV)')
plt.gca().yaxis.set_major_locator(plt.NullLocator())
plt.legend()
plt.draw()

lat = np.array(data['lattice'])
rlat = rec_lat(lat)
print "reciprocal lattice:"
for x in rlat:
    print ("%12.8lf "*3)%tuple(x)
x,y,z = np.array(rlat )
cut=.65
xmin,ymin,_ = -(x+y)*cut 
xmax,ymax,_ = +(x+y)*cut

#plot excitons
fig = plt.figure(figsize=(16,8))
plt.rc('text', usetex=True)
plt.rc('font', family='serif',serif="Computer Modern Roman",size=20)

sorted_excitons = sorted(data['excitons'],key=lambda x: x['energy'])

for n,exciton in enumerate(sorted_excitons[:8]):
    ax = plt.subplot(2,4,n+1)
    ax.set_aspect('equal', 'datalim')
    w = np.array(exciton['weights'])
    qpt = np.array(exciton['qpts'])
    ax.scatter(qpt[:,0], qpt[:,1], marker='H', s=10, color=[cmap(sqrt(c)) for c in w], label="e: %lf"%exciton['energy'])
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    plt.gca().yaxis.set_major_locator(plt.NullLocator())
    plt.gca().xaxis.set_major_locator(plt.NullLocator())
    ax.legend(prop={'size':12})

plt.draw()
plt.show()
