#from __future__ import print_function, division
from qepy import *
from yambopy import *
import matplotlib.pyplot as plt

npoints = 20

fig = plt.figure(figsize=(4,6))
ax  = fig.add_axes( [ 0.15, 0.15, 0.80, 0.80 ])

# Define path in reduced coordinates using Class Path

nk = 50
path = Path([ [[0.5,0.5,0.5],'$L$'         ],
              [[0.0,0.0,0.0],'$\Gamma$'    ],
              [[0.0,0.5,0.0],'$X$'         ],
              [[3./8.,3./4.,3./8],'$K$' ],
              [[0.0,0.0,0.0],'$\Gamma$'    ]], [nk,int(nk*2./sqrt(3.0)),int(nk/sqrt(6.)),int(nk*3./sqrt(3.0)/2.)] )

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
states = [1,2]

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
