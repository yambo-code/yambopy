# Copyright (C) 2018 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
# Tutorial File of Yambopy Tasks
# 
# 1. PwRelaxTask is organized in:
#    (a) Atomic relaxation
#    (b) Lattice relaxation
#    (c) SCF calculation
# 2. PwBandsTaks is organized in:
#    (a) SCF   calculation
#    (a) BANDS calculation
# 3. PwNscfTaks is organized in:
#    (a) SCF   calculation
#    (a) NSCF  calculation
#    (a) P2Y   calculation
from numpy import sqrt
import argparse
import os
import shutil
from yambopy.data.structures import BN 
from qepy.lattice import Path
from qepy.matdyn import Matdyn
from yambopy.io.factories import PwNscfTasks, PwBandsTasks, PwRelaxTasks
from yambopy.flow import YambopyFlow, PwTask, E2yTask, YamboTask
from schedulerpy import Scheduler
from yambopy import yambopyenv

nscf_bands = 60
kpoints = [6,6,1]
nscf_kpoints = [12,12,1]
ecut = 20
npoints = 10
pseudo_dir = '../pseudos'
path_kpoints = Path([ [[0.0, 0.0, 0.0],'G'],
                      [[0.5, 0.0, 0.0],'M'],
                      [[1./3,1./3,0.0],'K'],
                      [[0.0, 0.0, 0.0],'G']], [int(npoints*2),int(npoints),int(sqrt(5)*npoints)])

def relax():
    qe_relax_atoms_task, qe_relax_cell_task, qe_scf_task = PwRelaxTasks(BN,kpoints,ecut,cell_dofree='2Dxy',pseudo_dir=pseudo_dir)

    relax_flow = YambopyFlow.from_tasks('relax_flow',[qe_relax_atoms_task,qe_relax_cell_task,qe_scf_task])
    relax_flow.create(agressive=True)
    relax_flow.run()

def bands():
    pw_scf,pw_bands = PwBandsTasks(BN,kpoints,ecut,nscf_bands,path_kpoints,spin="spinor",pseudo_dir=pseudo_dir)
    bands_flow = YambopyFlow.from_tasks('bands_flow',[pw_scf,pw_bands])
    bands_flow.create(agressive=True)
    bands_flow.run()

def nscf():
    pw_scf,pw_nscf,p2y_task = PwNscfTasks(BN,kpoints,ecut,nscf_bands,nscf_kpoints,spin="spinor",pseudo_dir=pseudo_dir)
    nscf_flow = YambopyFlow.from_tasks('nscf_flow',[pw_scf,pw_nscf,p2y_task])
    nscf_flow.create(agressive=True)
    nscf_flow.run()
    print(nscf_flow)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Choose Yambopy Task.')
    parser.add_argument('-r' ,'--relax',       action="store_true", help='Structural relaxation Task')
    parser.add_argument('-b' ,'--bands',       action="store_true", help='Scf and Bands calculation Task')
    parser.add_argument('-n' ,'--nscf',        action="store_true", help='Scf, Nscf and p2y calculation Task')
    args = parser.parse_args()

if args.relax:
   relax()
if args.bands:
   bands()
if args.nscf:
   nscf()

