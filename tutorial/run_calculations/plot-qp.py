#
# Author: Alejandro Molina-Sanchez
#
# Example of YamboQPDB Class 
#
from qepy import *
from yambopy import *
import matplotlib.pyplot as plt


# Read Lattice information from SAVE
lat  = YamboLatticeDB.from_db_file(filename='gw/SAVE/ns.db1')

npoints = 10
# Brillouin-zone path on the lattice of the databases (hexagonal BN)
bz = BrillouinZone(cell=lat.lat, path_string='GMKG',
                   intervals=[int(npoints*2),int(npoints),int(np.sqrt(5)*npoints)])
# Read QP database
ydb  = YamboQPDB.from_db(filename='ndb.QP',folder='gw/yambo')


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
ks_bs_0, qp_bs_0 = ydb.get_bs_path(lat,bz)

fig = plt.figure(figsize=(4,5))
ax = fig.add_axes( [ 0.20, 0.20, 0.70, 0.70 ])

ks_bs_0.set_fermi(n_top_vb)
qp_bs_0.set_fermi(n_top_vb)

ks_bs_0.plot_ax(ax,legend=True,color_bands='r',label='KS')
qp_bs_0.plot_ax(ax,legend=True,color_bands='b',label='QP-GW')

plt.show()

# 3. Interpolation of KS and QP eigenvalues

ks_bs, qp_bs = ydb.interpolate(lat,bz,what='QP+KS',lpratio=20)

fig = plt.figure(figsize=(4,5))
ax = fig.add_axes( [ 0.20, 0.20, 0.70, 0.70 ])

ks_bs.set_fermi(n_top_vb)
qp_bs.set_fermi(n_top_vb)

ks_bs.plot_ax(ax,legend=True,color_bands='r',label='KS')
qp_bs.plot_ax(ax,legend=True,color_bands='b',label='QP-GW')

plt.show()

# 4. Comparison of not-interpolaed and  interpolated eigenvalues

fig = plt.figure(figsize=(4,5))
ax = fig.add_axes( [ 0.20, 0.20, 0.70, 0.70 ])

ks_bs_0.plot_ax(ax,legend=True,color_bands='r',label='KS')
qp_bs_0.plot_ax(ax,legend=True,color_bands='b',label='QP-GW')
ks_bs.plot_ax(ax,legend=True,color_bands='g',label='KS')
qp_bs.plot_ax(ax,legend=True,color_bands='k',label='QP-GW')

plt.show()
