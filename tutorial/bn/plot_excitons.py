from yambopy import *
import json
import numpy as np
from math import sqrt
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
print "nexitons", len(data['excitons'])
plt.plot(get_var(data,['E/ev','E/ev[1]']), get_var(data,['EPS-Im[2]' ]),label='BSE',lw=4)
plt.plot(get_var(data,['E/ev','E/ev[1]']), get_var(data,['EPSo-Im[4]']),label='EPS',lw=4)
for n,exciton in enumerate(data['excitons']):
    plt.axvline(exciton['energy'])
plt.draw()


#
# plot excitons
#
cmap = plt.get_cmap("gist_heat_r")
fig = plt.figure(figsize=(12,4))

sorted_excitons = sorted(data['excitons'],key=lambda x: x['energy'])

for n,exciton in enumerate(sorted_excitons):
    #get data
    w   = np.array(exciton['weights'])
    qpt = np.array(exciton['qpts'])

    #plot
    ax = plt.subplot(2,5,n+1)
    ax.scatter(qpt[:,0], qpt[:,1], s=20, c=w, cmap=cmap, lw=0, label="e: %lf"%exciton['energy'])

    # axis
    plt.xlim([-cut,cut])
    plt.ylim([-cut,cut])
    plt.gca().yaxis.set_major_locator(plt.NullLocator())
    plt.gca().xaxis.set_major_locator(plt.NullLocator())
    ax.set_aspect('equal', 'datalim')

#plt.subplot_tool()
plt.draw()
plt.show()
