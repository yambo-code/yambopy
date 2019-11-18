#from __future__ import print_function, division
from qepy import *
from yambopy import *
import matplotlib.pyplot as plt

npoints = 5

fig = plt.figure(figsize=(4,4))
ax  = fig.add_axes( [ 0.25, 0.10, 0.55, 0.45 ])

#define path in reduced coordinates using Class Path

path = Path([ [[  0.0,  0.0,  0.0],'$\Gamma$'],
              [[  0.5,  0.0,  0.0],'M'],
              [[1./3.,1./3.,  0.0],'K'],
              [[  0.0,  0.0,  0.0],'$\Gamma$']], [int(npoints*2),int(npoints),int(sqrt(5)*npoints)] )

# Load databases

# SAVE database
save = YamboSaveDB.from_db_file(folder='bse_flow/t0/SAVE')

# Lattice information
lat  = YamboLatticeDB.from_db_file(filename='bse_flow/t0/SAVE/ns.db1')

# Exciton database
yexc = YamboExcitonDB.from_db_file(lat,filename='ndb.BS_diago_Q01',folder='bse_flow/t0/run')

# Plot exciton weights in band structure 


exc_bands = yexc.get_exciton_bs(save,path,[1,2],size=1.0)
exc_bands.plot_ax(ax,color_bands='grey',c_weights='red')

plt.show()

exc_bands_inter = yexc.interpolate(save,path,[1,2],lpratio=5,f=None,size=0.05,verbose=True)
exc_bands_inter.plot_ax(ax,color_bands='grey',c_weights='red',alpha_weights=1.0,c_label='$X_1$')

plt.show()

'''
print("case of reading the eigenvalues from a qpDB DB")
#need to have the ndb.QP file in the same folder where this script is
yqp = YamboQPDB.from_db()
ax = plt.gca()
yexc.plot_exciton_bs(ax, yqp, path, (1,2,), space='bands')
#plt.savefig('exciton_bs_qp.pdf')
plt.show()

if 0:
    print("plot exciton in the brillouin zone")
    kpoints, amplitude, phase = yexc.get_amplitudes_phases((1,0,))
    for n,k in enumerate(kpoints):
        x,y,z = k
        plt.text(x,y,n)
    plt.scatter(kpoints[:,0],kpoints[:,1],c=amplitude,s=65,marker='H')
    ax = plt.axes()
    ax.set_aspect('equal', 'datalim')
    plr.savefig('exciton_bz.pdf')
    plt.show()

if 1:
    print("case of reading the eigenvalues from a saveDB")
    ax = plt.gca()
    yexc.plot_exciton_bs(ax, ysave, path, (1,2,), args_plot={'c':'g'},space='bands')
    plt.savefig('exciton_bs.pdf')
    plt.show()

if 1:
    print("case of reading the eigenvalues from a qpDB DB")
    #nee to have the ndb.QP file in the same folder where this script is
    yqp = YamboQPDB()
    ax = plt.gca()
    yexc.plot_exciton_bs(ax, yqp, path, (1,2,), space='bands')
    plt.savefig('exciton_bs_qp.pdf')
    plt.show()

#old method
if 0:
    yw = YamboExcitonWeight('o-yambo.exc_weights_at_1')
    ax = plt.gca()
    yw.plot_exciton_bs(ax,path,(1,2),space='bands')
    plt.show()
'''
