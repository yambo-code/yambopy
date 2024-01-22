#
# Author: Alejandro Molina-Sanchez
#
# Example of YamboQPDB Class 
#
from qepy import *
from yambopy import *
import matplotlib.pyplot as plt


# Define path in reduced coordinates using Class Path
npoints = 10
path = Path([ [[  0.0,  0.0,  0.0],'$\Gamma$'],
              [[  0.5,  0.0,  0.0],'M'],
              [[1./3.,1./3.,  0.0],'K'],
              [[  0.0,  0.0,  0.0],'$\Gamma$']], [int(npoints*2),int(npoints),int(sqrt(5)*npoints)] )

# Read Lattice information from SAVE
lat  = YamboSaveDB.from_db_file(folder='SAVE',filename='ns.db1')
# Read QP database
ydb  = YamboQPDB.from_db(filename='ndb.QP',folder='qp-gw')


# 1. Find scissor operator for valence and conduction bands

fig = plt.figure(figsize=(6,4))
ax  = fig.add_axes( [ 0.20, 0.20, 0.70, 0.70 ])
ax.set_xlabel('$E_{KS}$')
ax.set_ylabel('$E_{GW}$')

n_top_vb = 4
ydb.plot_scissor_ax(ax,n_top_vb)

plt.show()

# 2. Plot of KS and QP eigenvalues NOT interpolated along the path
n_top_vb = 3
ks_bs_0, qp_bs_0 = ydb.get_bs_path(lat,path)

fig = plt.figure(figsize=(4,5))
ax = fig.add_axes( [ 0.20, 0.20, 0.70, 0.70 ])

ks_bs_0.plot_ax(ax,legend=True,c_bands='r',label='KS')
qp_bs_0.plot_ax(ax,legend=True,c_bands='b',label='QP-GW')

plt.show()

# 3. Interpolation of KS and QP eigenvalues

ks_bs, qp_bs = ydb.interpolate(lat,path,what='QP+KS',lpratio=20)

fig = plt.figure(figsize=(4,5))
ax = fig.add_axes( [ 0.20, 0.20, 0.70, 0.70 ])

ks_bs.plot_ax(ax,legend=True,c_bands='r',label='KS')
qp_bs.plot_ax(ax,legend=True,c_bands='b',label='QP-GW')

plt.show()

# 4. Comparison of not-interpolaed and  interpolated eigenvalues

fig = plt.figure(figsize=(4,5))
ax = fig.add_axes( [ 0.20, 0.20, 0.70, 0.70 ])

ks_bs_0.plot_ax(ax,legend=True,c_bands='r',label='KS')
qp_bs_0.plot_ax(ax,legend=True,c_bands='b',label='QP-GW')
ks_bs.plot_ax(ax,legend=True,c_bands='g',label='KS int.')
qp_bs.plot_ax(ax,legend=True,c_bands='k',label='QP-GW int.')

plt.show()
