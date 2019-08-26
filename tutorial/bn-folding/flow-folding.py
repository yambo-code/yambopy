# Copyright (C) 2019 Henrique Pereira Coutada Miranda and Alejandro Molina-Sanchez
# All rights reserved.
#
# This file is part of yambopy
#
# Tutorial File of Unfolding Class 
# 
# 1. Calculations of Bands of the Primitive (PC) and Super Cell (SC) 
# 2. Unfolding of the SC onto the PC 

from numpy import sqrt
import argparse
import os
import shutil
from yambopy.data.structures import BN 
from qepy.lattice import Path
from yambopy.io.factories import PwNscfTasks, PwBandsTasks, PwRelaxTasks
from yambopy.flow import YambopyFlow, PwTask, E2yTask, YamboTask
from schedulerpy import Scheduler
from yambopy import yambopyenv
from qepy import PwXML
from qepy import Unfolding
from numpy import array, dot, sqrt, cross
import matplotlib.pyplot as plt

nscf_bands_pc = 6
nscf_bands_sc = 24
kpoints_pc = [6,6,1]
kpoints_sc = [3,3,1]
ecut = 20
npoints = 10

lattice_sc = dict(ibrav=4,celldm1=9.4,celldm3=1.27659)
atypes  = dict(B=[ 10.811,"B.pbe-mt_fhi.UPF"],
               N=[14.0067,"N.pbe-mt_fhi.UPF"])

atoms_sc  = [['N' ,  [0.0000000000  , 0.0000000000  ,-0.0000000000 ]],
             ['N' ,  [0.5000000000  , 0.0000000000  , 0.0000000000 ]],
             ['N' ,  [0.0000000000  , 0.5000000000  , 0.0000000000 ]],
             ['N' ,  [0.5000000000  , 0.5000000000  , 0.0000000000 ]],
             ['B' ,  [0.1666666667  , 0.3333333333  , 0.0000000000 ]],
             ['B' ,  [0.6666666667  , 0.3333333333  , 0.0000000000 ]],
             ['B' ,  [0.1666666667  , 0.8333333333  , 0.0000000000 ]],
             ['B' ,  [0.6666666667  , 0.8333333333  , 0.0000000000 ]]]

#occ = dict(occupations='smearing', smearing='mp', degauss=0.01)
BN_sc = dict(lattice=lattice_sc,atypes=atypes,atoms=atoms_sc)

path_kpoints_pc = Path([ [[0.0,  0.0, 0.0],'G'],
                         [[0.5,  0.0, 0.0],'M'],
                         [[1./3, 1./3,0.0],'K'],
                         [[0.0,  0.0, 0.0],'G']], [int(npoints*2),int(npoints),int(sqrt(5)*npoints)])

path_kpoints_sc = Path([ [[0.0,  0.0, 0.0],'G'],
                         [[1.0,  0.0, 0.0],'M'],
                         [[2./3, 2./3,0.0],'K'],
                         [[0.0,  0.0, 0.0],'G']], [int(npoints*2),int(npoints),int(sqrt(5)*npoints)])


def bands():

    pw_scf_pc,pw_bands_pc = PwBandsTasks(BN,kpoints_pc,ecut,nscf_bands_pc,path_kpoints_pc,spin="spinor",pseudo_dir="../../../bn/pseudos")
    bands_flow = YambopyFlow.from_tasks('bands_pc',[pw_scf_pc,pw_bands_pc])
    bands_flow.create(agressive=True)
    bands_flow.run()

    pw_scf_sc,pw_bands_sc = PwBandsTasks(BN_sc,kpoints_sc,ecut,nscf_bands_sc,path_kpoints_sc,spin="spinor",pseudo_dir="../../../bn/pseudos")
    bands_flow = YambopyFlow.from_tasks('bands_sc',[pw_scf_sc,pw_bands_sc])
    bands_flow.create(agressive=True)
    bands_flow.run()

## Unfolding Part
def plot():

    prefix_pc = 'pw'
    prefix_sc = 'pw'

    pc = PwXML(prefix=prefix_pc,path='bands_pc/t0')
    sc = PwXML(prefix=prefix_sc,path='bands_sc/t0')

    fold = Unfolding(prefix_pc=prefix_pc,path_pc='bands_pc/t0',prefix_sc=prefix_sc,path_sc='bands_sc/t0',spin="spinor")
    
    ax = plt.subplot(1,1,1)
    fold.plot_eigen_ax(ax,path=path_kpoints_pc)
    plt.show()

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Choose Yambopy Task.')
    parser.add_argument('-b' ,'--bands', action="store_true", help='Scf and Bands calculation Task')
    parser.add_argument('-p' ,'--plot',  action="store_true", help='Scf, Nscf and p2y calculation Task')
    args = parser.parse_args()

if args.bands:
   bands()
if args.plot:
   plot()

