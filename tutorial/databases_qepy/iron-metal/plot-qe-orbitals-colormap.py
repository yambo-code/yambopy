from qepy import  ProjwfcXML, ProjwfcIn
from yambopy import BrillouinZone
import matplotlib.pyplot as plt

# k-points map
npoints = 50
intervals=[npoints,npoints,npoints,npoints,npoints]
extra_points = { 'H': [0.0, 0.0, 1.0], 'N': [1/2, 0.0, 1/2], 'P': [1/2, 1/2, 1/2] }
path = 'GHNGPN'

bz = BrillouinZone(
    path_string=path,
    extra_points=extra_points,
    intervals=intervals
)

atom_s = [8]
atom_p = [0,1,2]
atom_d = [3,4,5,6,7]

# Class Projwfc
# Class to run projwfc.x and create 
# atomic_proj.xml (comment if already done)
'''
proj = ProjwfcIn(prefix='pw')
proj.run(folder='bands')
'''

# Class ProjwfcXML
# Atom-projected band structure. Colormap
band = ProjwfcXML(prefix='pw',path='bands')

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

band.plot_eigen(ax,path_kpoints=bz.get_indices(),cmap='viridis',cmap2='rainbow',selected_orbitals=p,selected_orbitals_2=d)

# Plot colormap
#
import matplotlib as mpl
cmap =plt.get_cmap('viridis')
cmap2=plt.get_cmap('rainbow')
bx  = fig.add_axes( [ 0.84, 0.10, 0.03, 0.80 ])
cx  = fig.add_axes( [ 0.93, 0.10, 0.03, 0.80 ])
norm = mpl.colors.Normalize(vmin=0.,vmax=1.)
cb1 = mpl.colorbar.ColorbarBase(bx, cmap=cmap, norm=norm,orientation='vertical',ticks=[0,1])
cb1.set_ticklabels(['d', 'p'])
cb2 = mpl.colorbar.ColorbarBase(cx, cmap=cmap2, norm=norm,orientation='vertical',ticks=[0,1])
cb2.set_ticklabels(['d', 'p'])

plt.show()
