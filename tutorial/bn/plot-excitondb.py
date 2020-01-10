#from __future__ import print_function, division
from qepy import *
from yambopy import *
import matplotlib.pyplot as plt

npoints = 20

fig = plt.figure(figsize=(4,6))
ax  = fig.add_axes( [ 0.15, 0.15, 0.80, 0.80 ])

# Define path in reduced coordinates using Class Path

path = Path([ [[  0.0,  0.0,  0.0],'$\Gamma$'],
              [[  0.5,  0.0,  0.0],'M'],
              [[1./3.,1./3.,  0.0],'K'],
              [[  0.0,  0.0,  0.0],'$\Gamma$']], [int(npoints*2),int(npoints),int(sqrt(5)*npoints)] )

# Load databases

# SAVE database
save = YamboSaveDB.from_db_file(folder='bse_calc/SAVE')

# Lattice information
lat  = YamboLatticeDB.from_db_file(filename='bse_calc/SAVE/ns.db1')

# Exciton database read from db file
yexc = YamboExcitonDB.from_db_file(lat,filename='ndb.BS_diago_Q01',folder='bse_calc/yambo')

print("Ground state energy: %lf" % yexc.eigenvalues[0].real )
print("Intensity: %lf" % (yexc.get_intensities()[0]+yexc.get_intensities()[1]) )
print("1st-excited state energy: %lf" % yexc.eigenvalues[2].real )
print("Intensity: %lf" % (yexc.get_intensities()[2]+yexc.get_intensities()[3]) )

# List of states to be merged
states = [3,4]

# 1. Plot exciton weights in band structure NOT interpolated

exc_bands = yexc.get_exciton_bs(save,path,states,size=1.0)
exc_bands.plot_ax(ax,color_bands='grey',c_weights='red')

plt.show()

# 2. Plot exciton weights in band structure INTERPOLATED

fig = plt.figure(figsize=(4,6))
ax  = fig.add_axes( [ 0.15, 0.15, 0.80, 0.80 ])

exc_bands_inter = yexc.interpolate(save,path,states,lpratio=5,f=None,size=0.5,verbose=True)

exc_bands_inter.plot_ax(ax,color_bands='grey',c_weights='red',alpha_weights=0.5,c_label='$X_1$')

plt.show()

# 3. Plot exciton weights in a 2D map of the BZ

from matplotlib.patches import Polygon

fig = plt.figure(figsize=(4,4))
ax  = fig.add_axes( [ 0.15, 0.15, 0.80, 0.80 ])
lattice = lat.rlat

x1 = 1./3*lattice[0][:2]+1./3*lattice[1][:2]
x2 =-1./3*lattice[0][:2]+2./3*lattice[1][:2]
x3 =-2./3*lattice[0][:2]+1./3*lattice[1][:2]
x4 = -x1
x5 = -x2
x6 = -x3
hexagon = [x1,x2,x3,x4,x5,x6]

yexc.plot_exciton_2D_ax(ax,states,mode='hexagon',limfactor=0.8,scale=160)
ax.add_patch(Polygon(hexagon,closed=True,fill=False,color='w',lw=1.0))

plt.show()

