# Copyright (C) 2019 Alejandro Molina Sanchez - Henrique PC Miranda 
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

import argparse
import os
import shutil
from yambopy.data.structures import Si
from qepy import PwXML
from qepy.lattice import Path
from qepy.matdyn import Matdyn
from yambopy.io.factories import PwNscfTasks, PwBandsTasks, PwRelaxTasks
from yambopy.flow import YambopyFlow, PwTask, E2yTask, YamboTask
from schedulerpy import Scheduler
from yambopy import yambopyenv

nscf_bands   = 10
kpoints      = [2,2,2]
nscf_kpoints = [4,4,4]
ecut = 30
path_kpoints = Path([ [[1.0,1.0,1.0],'$\Gamma$'],
                      [[0.0,0.5,0.5],'$X$'],
                      [[0.0,0.0,0.0],'$\Gamma$'],
                      [[0.5,0.0,0.0],'$L$']], [20,20,20])
pseudo_dir = '../../pseudos'

def relax():
    qe_relax_atoms_task, qe_relax_cell_task, qe_scf_task = PwRelaxTasks(Si,kpoints,ecut,cell_dofree='all',pseudo_dir=pseudo_dir)

    relax_flow = YambopyFlow.from_tasks('relax_flow',[qe_relax_atoms_task,qe_relax_cell_task,qe_scf_task])
    relax_flow.create(agressive=True)
    relax_flow.run()

def bands():
    pw_scf,pw_bands = PwBandsTasks(Si,kpoints,ecut,nscf_bands,path_kpoints,pseudo_dir=pseudo_dir)
    bands_flow = YambopyFlow.from_tasks('bands_flow',[pw_scf,pw_bands])
    bands_flow.create(agressive=True)
    bands_flow.run()

def nscf():
    pw_scf,pw_nscf,p2y_task = PwNscfTasks(Si,kpoints,ecut,nscf_bands,nscf_kpoints,pseudo_dir=pseudo_dir)
    nscf_flow = YambopyFlow.from_tasks('nscf_flow',[pw_scf,pw_nscf,p2y_task])
    nscf_flow.create(agressive=True)
    nscf_flow.run()
    print(nscf_flow)

def plot_bands(show=True):
    xml = PwXML(prefix='pw',path='bands_flow/t0')
    xml.plot_eigen(path=path_kpoints,show=show)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Choose Yambopy Task.')
    parser.add_argument('-r' ,'--relax',       action="store_true", help='Structural relaxation Task')
    parser.add_argument('-b' ,'--bands',       action="store_true", help='Scf and Bands calculation Task')
    parser.add_argument('-n' ,'--nscf',        action="store_true", help='Scf, Nscf and p2y calculation Task')
    parser.add_argument('-p' ,'--plot',        action="store_true", help='Plot bands')
    args = parser.parse_args()

if args.relax:
   relax()
if args.bands:
   bands()
if args.nscf:
   nscf()
if args.plot:
   plot_bands()

