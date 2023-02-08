from __future__ import print_function, division
import sys
import argparse
from qepy import *
from schedulerpy import *
from math import sqrt

# Matplotlib options
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(5,7))
ax  = fig.add_axes( [ 0.12, 0.10, 0.70, 0.80 ])

# k-points map
npoints = 50
path_kpoints = Path([ [[0.0, 0.0, 0.0],'$\Gamma$'],
                      [[0.5, 0.0, 0.0],'M'],
                      [[1./3,1./3,0.0],'K'],
                      [[0.0, 0.0, 0.0],'$\Gamma$']], [int(npoints*2),int(npoints),int(sqrt(5)*npoints)])

atom_1 = list(range(16))
atom_2 = list(range(16,32))

# Class Projwfc
# Class to run projwfc.x and create 
# atomic_proj.xml (comment if already done)
'''
proj = ProjwfcIn(prefix='bn')
proj.run(folder='bands')
'''

# Class ProjwfcXML
# Atom-projected band structure
band = ProjwfcXML(prefix='bn',path='bands',qe_version='7.0')
band.plot_eigen(ax,path_kpoints=path_kpoints,cmap='viridis',selected_orbitals=atom_1,selected_orbitals_2=atom_2)

# Plot colormap
#
import matplotlib as mpl
cmap=plt.get_cmap('viridis')
bx  = fig.add_axes( [ 0.88, 0.10, 0.05, 0.80 ])
norm = mpl.colors.Normalize(vmin=0.,vmax=1.)
cb1 = mpl.colorbar.ColorbarBase(bx, cmap=cmap, norm=norm,orientation='vertical',ticks=[0,1])
cb1.set_ticklabels(['B', 'N'])

plt.show()
