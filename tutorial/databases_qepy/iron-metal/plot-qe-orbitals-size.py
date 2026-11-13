from qepy import ProjwfcXML,ProjwfcIn
from yambopy import BrillouinZone
import matplotlib.pyplot as plt

# k-points map
npoints = 50
path = 'GHNGPN'
extra_points = {'G': [0.0, 0.0, 0.0], 'H': [0.0, 0.0, 1.0], 'N': [1/2, 0.0, 1/2], 'P': [1/2, 1/2, 1/2]}
intervals = [npoints,npoints,npoints,npoints,npoints]

bz = BrillouinZone(
    path_string=path,
    extra_points=extra_points,
    intervals=intervals
)

# Class Projwfc
# Class to run projwfc.x and create 
# atomic_proj.xml (comment if already done)
'''
proj = ProjwfcIn(prefix='pw')
proj.run(folder='bands')
'''

# Class ProjwfcXML
# Atom-projected band structure. Size
band = ProjwfcXML(prefix='pw',path='bands')
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

#band.plot_eigen(ax,bz,selected_orbitals=s,color='pink',color_2='black')
#band.plot_eigen(ax,bz,selected_orbitals=p,color='green',color_2='orange')
band.plot_eigen(ax,path_kpoints=bz.get_indices(),selected_orbitals=d,color='red',color_2='blue')

plt.show()
