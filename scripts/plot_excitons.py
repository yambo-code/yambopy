# Copyright (C) 2015 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
# Plot the excitonic weights on the brillouin zone using ypp
from yambopy import *
import json
import numpy as np
from math import sqrt, ceil
import matplotlib.pyplot as plt

ang2bohr = 1.889725989
cut = 0.2 #set the plot limits

def get_var(dictionary,variables):
    """
    To have compatibility with different version of yambo
    We provide a list of different possible tags
    """
    for var in variables:
        if var in dictionary:
            return dictionary[var]
    raise ValueError( 'Could not find the variables %s in the output file'%str(variables) )
#
# read file
#
f = open('absorptionspectra.json')
data = json.load(f)
f.close()

#
# plot the absorption spectra
#
nexcitons = len(data['excitons'])
print "nexitons", nexcitons
plt.plot(get_var(data,['E/ev','E/ev[1]']), get_var(data,['EPS-Im[2]' ]),label='BSE',lw=4)
plt.plot(get_var(data,['E/ev','E/ev[1]']), get_var(data,['EPSo-Im[4]']),label='EPS',lw=4)
for n,exciton in enumerate(data['excitons']):
    plt.axvline(exciton['energy'])
plt.xlabel('$\\omega$ (eV)')
plt.ylabel('Intensity arb. units')
plt.draw()

#
# plot excitons
#

#dimensions
nx = int(ceil(sqrt(nexcitons)))
ny = int(ceil(nexcitons*1.0/nx))
print "cols:",nx
print "rows:",ny
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

#plt.subplot_tool()
plt.draw()
plt.show()
