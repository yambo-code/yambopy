from yambopy import YamboLatticeDB, YamboQPDB # Load yambo netcdf databases
from yambopy import BrillouinZone # Define path in k-space
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt

# Read Lattice information from SAVE
## Note: we do not expand the kpts because QP database is in the IBZ
lat = YamboLatticeDB.from_db_file(filename='SAVE/ns.db1',Expand=False)

# Define path in reduced coordinates using BrillouinZone class
path_string = 'GMKG'

bz = BrillouinZone(
    cell=lat.lat,
    path_string=path_string,
    npoints=100,
    optimise=True,  # or options, e.g. {'method': 'brute', 'parallel_segments': False}
    sym_red=lat.sym_red
)

# Read QP database
ydb  = YamboQPDB.from_db(filename='ndb.QP',folder='qp-gw')
n_top_vb = 3 # Top valence band index starting from 0
fermie = np.max(ydb.eigenvalues_qp[:,n_top_vb]) # Energy shift to top valence

# 1. Find scissor operator for valence and conduction bands

fig = plt.figure(figsize=(6,4))
ax  = fig.add_axes( [ 0.20, 0.20, 0.70, 0.70 ])
ax.set_xlabel('$E_{KS}$')
ax.set_ylabel('$E_{GW}$')

ydb.plot_scissor_ax(ax,n_top_vb+1)

plt.show()

# 2. Plot of KS and QP eigenvalues NOT interpolated along the path

ks_bs_0, qp_bs_0 = ydb.get_bs_path(lat,bz)

fig = plt.figure(figsize=(4,5))
ax = fig.add_axes( [ 0.20, 0.20, 0.70, 0.70 ])

ks_bs_0.plot_ax(ax,legend=True,c_bands='r',label='KS')
qp_bs_0.plot_ax(ax,legend=True,c_bands='b',fermie=fermie,label='QP-GW')

plt.show()

# 3. Interpolation of KS and QP eigenvalues

ks_bs, qp_bs = ydb.interpolate(lat,bz,what='QP+KS',lpratio=20)

fig = plt.figure(figsize=(4,5))
ax = fig.add_axes( [ 0.20, 0.20, 0.70, 0.70 ])

ks_bs.plot_ax(ax,legend=True,c_bands='r',label='KS')
qp_bs.plot_ax(ax,legend=True,c_bands='b',fermie=fermie,label='QP-GW')

plt.show()

# 4. Comparison of not-interpolated and  interpolated eigenvalues

fig = plt.figure(figsize=(4,5))
ax = fig.add_axes( [ 0.20, 0.20, 0.70, 0.70 ])

ks_bs_0.plot_ax(ax,legend=True,c_bands='r',label='KS')
qp_bs_0.plot_ax(ax,legend=True,c_bands='b',fermie=fermie,label='QP-GW')
ks_bs.plot_ax(ax,legend=True,c_bands='g',label='KS int.')
qp_bs.plot_ax(ax,legend=True,c_bands='k',fermie=fermie,label='QP-GW int.')

plt.show()
