#
# Author: Alejandro Molina-Sanchez
#
# Example of YamboQPDB Class 
#
from qepy import *
from yambopy import *
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(4,6))
#ax  = fig.add_axes( [ 0.15, 0.15, 0.80, 0.80 ])

# Define path in reduced coordinates using Class Path
npoints = 10
path = Path([ [[  0.0,  0.0,  0.0],'$\Gamma$'],
              [[  0.5,  0.0,  0.0],'M'],
              [[1./3.,1./3.,  0.0],'K'],
              [[  0.0,  0.0,  0.0],'$\Gamma$']], [int(npoints*2),int(npoints),int(sqrt(5)*npoints)] )

# Read Lattice information from SAVE
lat  = YamboSaveDB.from_db_file(folder='gw_flow/t0/SAVE',filename='ns.db1')
# Read QP database
y    = YamboQPDB.from_db(filename='ndb.QP',folder='gw_flow/t0/run')

#print(y.eigenvalues_dft)

#ks_bs, qp_bs = y.get_bs()
ks_bs, qp_bs = y.interpolate(lat,path,what='KS',lpratio=20)

ax = fig.add_axes( [ 0.10, 0.15, 0.40, 0.80 ])

ks_bs.plot_ax(ax)

plt.show()
