from numpy import sqrt,pi, array
import argparse
import os
import shutil
from qepy.lattice import Path
from qepy.matdyn import Matdyn
from yambopy.io.factories import PwNscfTasks, PwBandsTasks, PwRelaxTasks
from yambopy.flow import YambopyFlow, PwTask, E2yTask, YamboTask
from schedulerpy import Scheduler
from yambopy import yambopyenv

#sch = Scheduler.factory(scheduler="slurm",ntasks=8,walltime="10:00:00")
sch = Scheduler.factory(scheduler="bash")
#sch.add_module("qe/6.1")

pseudo_dir = '.'

from_bohr_to_ang = 0.529177249

kpoints, nscf_kpoints = [8,8,8], [4,4,4]
shiftk = [1,1,1]
nscf_bands, path_bands = 10, 10 
ecut = 45

npoints = 50
path_kpoints = Path([ [[0.0, 0.0, 0.0 ],'G'],
                      [[0.0, 0.0, 1.0 ],'H'],
                      [[1./2,0.0,1./2.],'N'],
                      [[0.0, 0.0, 0.0 ],'G'],
                      [[1./2, 1./2, 1./2 ],'P'],
                      [[1./2,0.0,1./2. ],'N']], [npoints,npoints,npoints,npoints,npoints])

lattice = dict(ibrav=3,celldm1=5.42)

atypes  = dict(Fe=[55.845  ,"Fe.rel-pbe-n-nc.UPF"])

atoms = [['Fe' , [  0.0000000000 , 0.0000000000 , 0.0 ]]]

occ  = dict(occupations='smearing', degauss=0.05)
elec = dict(mixing_mode='local-TF',mixing_beta=0.3)
ion  = dict(ion_dynamics='damp')
ce   = dict(cell_dynamics='sd')

material = dict(lattice=lattice,atypes=atypes,atoms=atoms,occupations=occ,electrons=elec,ions=ion,cell=ce)

# Relax Class
qe_atoms_task, qe_cell_task, qe_scf_task = PwRelaxTasks(material,kpoints,ecut,scf_conv_thr=1.e-8,cell_dofree='all',starting_magnetization=[2],spin="polarized",pseudo_dir=pseudo_dir)

# Bands Class
qe_bands_scf, qe_bands_bands = PwBandsTasks(material,kpoints,ecut,path_bands,path_kpoints,scf_conv_thr=1.e-8,starting_magnetization=[2],spin="polarized",pseudo_dir=pseudo_dir)
qe_bands_scf.pwinput.set_kpoints(kpoints,shiftk)

# Nscf Class
qe_nscf_scf, qe_nscf_nscf, p2y_task = PwNscfTasks(material,kpoints,ecut,nscf_bands,nscf_kpoints,scf_conv_thr=1.e-8,starting_magnetization=[2],spin="polarized",pseudo_dir=pseudo_dir)

def relax():
    # Scheduler
    qe_atoms_task.scheduler_setup(sch)
    qe_cell_task.scheduler_setup(sch)
    qe_scf_task.scheduler_setup(sch)

    relax_flow = YambopyFlow.from_tasks('relax',[qe_atoms_task,qe_cell_task,qe_scf_task])
    relax_flow.create(agressive=True)
    relax_flow.run()

def bands():

    # Scheduler
    qe_bands_scf.scheduler_setup(sch)
    qe_bands_bands.scheduler_setup(sch)

    bands_flow = YambopyFlow.from_tasks('bands',[qe_bands_scf,qe_bands_bands])
    bands_flow.create(agressive=True)
    bands_flow.run()

def nscf():

    # Scheduler
    qe_nscf_scf.scheduler_setup(sch)
    qe_nscf_nscf.scheduler_setup(sch)
    p2y_task.scheduler_setup(sch)

    nscf_flow = YambopyFlow.from_tasks('nscf3',[qe_nscf_scf,qe_nscf_nscf,p2y_task])
    nscf_flow.create(agressive=True)
    nscf_flow.run()

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
