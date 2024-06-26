import sys
import argparse
from qepy import *
from schedulerpy import *

# Matplotlib options
import matplotlib.pyplot as plt

# k-points map
npoints = 50
path_kpoints = Path([ [[0.0, 0.0, 0.0 ],'G'],
                      [[0.0, 0.0, 1.0 ],'H'],
                      [[1./2,0.0,1./2.],'N'],
                      [[0.0, 0.0, 0.0 ],'G'],
                      [[1./2, 1./2, 1./2 ],'P'],
                      [[1./2,0.0,1./2. ],'N']], [npoints,npoints,npoints,npoints,npoints])

# Class Projwfc
# Class to run projwfc.x and create 
# atomic_proj.xml (comment if already done)
'''
proj = ProjwfcIn(prefix='pw')
proj.run(folder='bands/t0')
'''

# Class ProjwfcXML
# Atom-projected band structure. Size
band = ProjwfcXML(prefix='pw',path='bands/t0')
# print info
print(band)

# Manual selection of the lists of states by inspecting projwfc output
#s = [8]
#p = [0,1,2]
#d = [3,4,5,6,7]

# Automatic selection of the states
s = band.get_states_helper(orbital_query=['s'])
p = band.get_states_helper(orbital_query=['p'])
d = band.get_states_helper(orbital_query=['d'])

fig = plt.figure(figsize=(5,7))
ax  = fig.add_axes( [ 0.12, 0.10, 0.70, 0.80 ])

band.plot_eigen(ax,path_kpoints=path_kpoints,selected_orbitals=s,color='red',color_2='blue')
band.plot_eigen(ax,path_kpoints=path_kpoints,selected_orbitals=p,color='green',color_2='orange')
band.plot_eigen(ax,path_kpoints=path_kpoints,selected_orbitals=d,color='pink',color_2='black')

plt.show()
