from numpy import sqrt,pi, array
import argparse
import os
import shutil
from qepy.lattice import Path
from qepy.matdyn import Matdyn
from yambopy.io.factories import PwNscfTasks, PwBandsTasks, PwRelaxTasks
from yambopy.flow import YambopyFlow, PwTask, E2yTask, YamboTask
from schedulerpy import Scheduler
from yambopy.env import yambopyenv

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

# Bands Class
qe_bands_scf, qe_bands_bands = PwBandsTasks(material,kpoints,ecut,path_bands,path_kpoints,scf_conv_thr=1.e-8,starting_magnetization=[2],spin="polarized",pseudo_dir=pseudo_dir)
qe_bands_scf.pwinput.set_kpoints(kpoints,shiftk)

# Scheduler
qe_bands_scf.scheduler_setup(sch)
qe_bands_bands.scheduler_setup(sch)

# Run Flow
bands_flow = YambopyFlow.from_tasks('bands',[qe_bands_scf,qe_bands_bands])
bands_flow.create(agressive=True)
bands_flow.run()
