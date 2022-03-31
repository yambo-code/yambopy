from __future__ import print_function, division
import sys
import argparse
from qepy import *
from schedulerpy import *
from math import sqrt

# k-points map
npoints = 50
path_kpoints = Path([ [[0.0, 0.0, 0.0],'$\Gamma$'],
                      [[0.5, 0.0, 0.0],'M'],
                      [[1./3,1./3,0.0],'K'],
                      [[0.0, 0.0, 0.0],'$\Gamma$']], [int(npoints*2),int(npoints),int(sqrt(5)*npoints)])

atom_1 = [range(16)]
atom_2 = [range(16,32)]

# Xml database reading
xml = PwXML(prefix='bn',path='bands')
xml.plot_eigen(path_kpoints)

#proj = ProjwfcIn(prefix='pw') #,folder='bands/t0')
#proj.run(folder='bands/t0')

# Atom-projected band structure

#import matplotlib.pyplot as plt
#fig = plt.figure(figsize=(5,7))
#ax  = fig.add_axes( [ 0.10, 0.10, 0.70, 0.80 ])

#band = ProjwfcXML(prefix='pw',path='bands/t0',qe_version='6.7')
#band.plot_eigen(ax,path_kpoints=path_kpoints,selected_orbitals=atom_1,selected_orbitals_2=atom_2)
#plt.show()

