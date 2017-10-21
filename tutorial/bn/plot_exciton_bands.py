from __future__ import print_function, division
from yambopy import *

#define path in reduced coordinates
path = [ [0.0, 0.0, 0.0],
         [0.5, 0.0, 0.0],
         [1./3,1./3,0.0],
         [0.0, 0.0, 0.0]]

#load databases
ysave = YamboSaveDB()
ylat  = YamboLatticeDB()
yexc = YamboExcitonDB(ylat,path='yambo')

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

