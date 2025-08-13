#from __future__ import print_function, division
from qepy import *
from yambopy import *
import matplotlib.pyplot as plt
import os

npoints = 20
bse_dir = 'bse-exc'
#bse_dir = 'bse-magn'
fig = plt.figure(figsize=(4,6))
ax  = fig.add_axes( [ 0.15, 0.15, 0.80, 0.80 ])

# Define path in reduced coordinates using Class Path

path = Path([ [[  0.0,  0.0,  0.0],r'$\Gamma$'],
              [[  0.5,  0.0,  0.0],'M'],
              [[1./3.,1./3.,  0.0],'K'],
              [[  0.0,  0.0,  0.0],r'$\Gamma$']], [int(npoints*2),int(npoints),int(sqrt(5)*npoints)] )

# Load databases

# Electron energies database
elec = YamboElectronsDB.from_db_file(folder='./SAVE')

# Lattice information
lat  = YamboLatticeDB.from_db_file(filename='./SAVE/ns.db1')

# Exciton database read from db file
if os.path.isfile(bse_dir+'/ndb.BS_diago_Q01'):
    yexc = YamboExcitonDB.from_db_file(lat,filename='ndb.BS_diago_Q01',folder=bse_dir)
if os.path.isfile(bse_dir+'/ndb.BS_diago_Q1'):
    yexc = YamboExcitonDB.from_db_file(lat,filename='ndb.BS_diago_Q1',folder=bse_dir)

print("Ground state energy: %lf" % yexc.eigenvalues[0].real )
print("Intensity: %lf" % (yexc.get_intensities()[0].real+yexc.get_intensities()[1].real) )
print("1st-excited state energy: %lf" % yexc.eigenvalues[2].real )
print("Intensity: %lf" % (yexc.get_intensities()[2].real+yexc.get_intensities()[3].real) )

# List of states to be merged
states = [1,2]

# 1. Plot exciton weights in band structure NOT interpolated
exc_bands_up, exc_bands_dn = yexc.get_exciton_bs_spin_pol(elec,path,states,size_up=1.0, size_dw = 1.0)
exc_bands_up.plot_ax(ax,c_bands='red',c_weights='orange')
exc_bands_dn.plot_ax(ax,c_bands='blue',c_weights='orange')
plt.savefig('plot1.png')
plt.show()

# 2. Plot exciton weights in band structure INTERPOLATED
#PART UNDER DEVELOPMENT
fig = plt.figure(figsize=(4,6))
ax  = fig.add_axes( [ 0.15, 0.15, 0.80, 0.80 ])

# In case of problems with the interpolation, try to increase lpratio
exc_bands_inter_up, exc_bands_inter_dn  = yexc.interpolate_spin_pol(elec,path,states,lpratio=10,f=None,size=0.5,verbose=True)
exc_bands_inter_up.plot_ax(ax,c_bands='red',c_weights='orange',alpha_weights=0.5,c_label='$X_1$')
#plt.savefig('plot2.png')
plt.show()


# 3. Plot exciton weights in a 2D map of the BZ

fig = plt.figure(figsize=(4,4))
ax  = fig.add_axes( [ 0.15, 0.15, 0.80, 0.80 ])

#yexc.plot_exciton_2D_ax(ax,states,mode='hexagon',limfactor=0.8,scale= 600)
yexc.plot_exciton_2D_ax(ax,states,limfactor=0.8,scale= 600)
#plt.savefig('plot3.png')
plt.show()
