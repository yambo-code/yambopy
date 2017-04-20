from __future__ import print_function
##############################################################################
#
# Author: Alejandro Molina-Sanchez
# Run real-time simulations with Yambo
#
# Warning: Real-time simulations requires several data folders for running
# properly. Before using these scripts compulsively, it is recommended
# to understand the different run levels.
#
# Instructions:
# The dictionary 'job' stores the variable of the calculation
# calculation : 'collisions', 'pump', 'dissipation'
# folder-col  : collision data
# folder-run  : results (only work if collisions have been previously calculated)
# DG          : True or False if we use the double-grid
# nodes       : cpu-dependent variable
# cores       :          "
#
##############################################################################
#from __future__ import print_function
from builtins import str
from yambopy     import *
from schedulerpy import *
import argparse

# Select the run-level : 'collision', 'pump', 'dissipation'
parser = argparse.ArgumentParser(description='Real-time simulation')
parser.add_argument('-c' ,'--collisions' ,action="store_true")
parser.add_argument('-p' ,'--pump'       ,action="store_true")
parser.add_argument('-d' ,'--dissipation',action="store_true")
args = parser.parse_args()

yambo_rt = 'yambo_rt'
folder   = 'rt-6x6'

# Generation of the input file

#os.system('cd %s; mkdir -p inputs' % folder)

if args.collisions:
  print('Collisions')
  run = YamboIn('%s -r -e -v hsex'%yambo_rt,folder=folder)
elif args.pump:
  print('Time-dependent with electric field')
  run = YamboIn('%s -q p'%yambo_rt,folder=folder)
elif args.dissipation:
  print('Time-dependent with electric field and electron-phonon scattering')
  run = YamboIn('%s -s p -q p'%yambo_rt,folder=folder)
else:
  print('Invalid calculation type')
  exit()

# Proportionality constant of the DOS(T)
def dos_temperature(temperature):
  ac, bc, cc = 0.506760, 0.000069, 0.000002
  av, bv, cv = 0.437618, 0.000070, 0.000002
  dos_c = ac + bc*temperature + cc*temperature**2
  dos_v = av + bv*temperature + cv*temperature**2
  return [dos_c,dos_v]

job = dict()
job['radiative']    = False
job['folder-col']   = 'col-hxc' # name of collisions folder (path folder/folder-col)
# the input file and run-folder names are defined at the end
job['folder-gkkp']  = 'gkkp'
job['DG']           = (False,'dg-60x60') # Double-grid (True|False) and DG folder
job['nodes']        = 1
job['cores']        = 12
job['ppn']          = job['cores']
job['threads']      = job['cores']*job['nodes']
job['name']         = 'MoS2-dynamics'
job['yambo_module'] = 'yambo/master-intel'
job['temperature']  = 0

# System Common variables
run['DBsIOoff']      = "GF P J"
#run['DBsIOoff']      = "GF CARRIERs"
run['RT_CPU']        = "%s.1.1.1"  % (job['nodes']*job['cores']) # [PARALLEL] CPUs for each role
run['RT_ROLEs']      = "k.b.q.qp"    # [PARALLEL] CPUs roles (k,b,q,qp)
run['HXC_Potential'] = "HARTREE+SEX" # [SC] SC HXC Potential
#run['RT_Threads']    = 1

# Collision variables
if args.collisions:
  run['FFTGvecs']      = [20  ,'Ha']
  run['EXXRLvcs']      = [1000,'mHa']  # Equivalent to BSENGBlk in the BSE run-level
  run['HARRLvcs']      = [20  ,'Ha']   # Equivalent to BSENGexx in the BSE run-level
  run['CORRLvcs']      = [1000,'mHa']  #
  run['COLLBands']     = [25,28]
  run['NGsBlkXs']      = [1000,'mHa']
  run['BndsRnXs']      = [1, 60]
  run['RandQpts']      = 1000000         # Only in layers
  run['RandGvec']      = [ 1 , 'RL' ]
  run['CUTGeo']        = 'box z'         # Only in layers
  run['CUTBox']        = [0,0,38]
  run['X_all_q_CPU']   = "1.%d.1.1" % (job['nodes']*job['cores'])     # [PARALLEL] CPUs for each role
  run['X_all_q_ROLEs'] = "q.k.c.v"       # [PARALLEL] CPUs roles (q,k,c,v)
  #run.arguments.append('ALLGHAR')
  run.write('%s/collisions.in'%folder)

# Common time-dependent variable
if args.pump or args.dissipation:
  run['RTBands']    = [25,28]
  #run['GfnQP_Wv']   = [0.0,dos_temperature(job['temperature'])[1],0.00]  # Only for dissipation
  #run['GfnQP_Wc']   = [0.0,dos_temperature(job['temperature'])[0],0.00]
  run['GfnQP_Wv']   = [0.025,0.0,0.0]
  run['GfnQP_Wc']   = [0.025,0.0,0.0]
  run['GfnQP_E']    = [1.04588, 1.00, 1.00]            # [EXTQP BSK BSS] E parameters  (c/v) eV|adim|adim
  # Time-propagation
  run['RTstep']     = [5.0,'as'] # Real-time step
  run['NETime']     = [600,'fs'] # Simulation duration
  run['Integrator'] = "RK2 RWA"
  run['IOtime']     = [ [ 50.0, 00.10, 1.00], 'fs' ]   # [RT] Time between to consecutive I/O (J,P,CARRIERs - GF - OUTPUT)
  #run['IOtime']     = [ [ 0.10, 0.10, 1.00], 'fs' ]   # [RT] Time between to consecutive I/O (J,P,CARRIERs - GF - OUTPUT)
  # Pump Pulse
  run['Field1_Int']      = [ 1E3 , 'kWLm2']    # Intensity pulse
  run['Field1_Dir']      = [0.0,1.0,0.0]       # Polarization pulse
  run['Field1_Dir_circ'] = [1.0,0.0,0.0]       # Polarization pulse (second axis for circular)
  run['Field1_pol']      = "linear"            # Polarization type (linear or circular)
  run['Field1_kind'] = "QSSIN"                 # [RT Field1] Kind(SIN|RES|ANTIRES|GAUSS|DELTA|QSSIN)
  run['Field1_Damp'] = [ 70.0,'fs']            # Damping (width of pulse)
  run['Field1_Freq'] = [[2.0,2.0],'eV']
  if job['radiative']:
    run.arguments.append('el_photon_scatt')

# Pumping with finite pulse and electron-phonon dissipation
if args.dissipation:
# Interpolation
  run['BoseTemp']        = [ job['temperature'], 'K']
  run['LifeExtrapSteps'] = [ [0.05,0.05], 'fs' ] # Extrapolation time for phonon-assisted dissipation
  run['ElPhModes']       = [ 1, 9]
  run.arguments.append('LifeExtrapolation')   # Commented:   Lifetimes are constant

# Run -- Creating

print('Collisions        ',job['folder-col'])
print('Number of nodes   ',job['nodes'])
print('Number of cores   ',job['cores'])

oarsub = oarsub(nodes=job['nodes'],core=job['cores'],dependent=0,name=job['name'],walltime="24:00:00")
oarsub.add_command('module load %s'%job['yambo_module'])
oarsub.add_command('export OMP_NUM_THREADS=1')

# Collisions
if args.collisions:
  oarsub.add_command('cd %s; mpirun -hostfile \$OAR_NODEFILE %s -F collisions.in -J %s -C %s'%(folder,yambo_rt,job['folder-col'],job['folder-col']))
  oarsub.write('%s/collisions.ll' % folder )
  oarsub.run()
  print('running yambo-collision')

# Time-dependent without dissipation
if args.pump:

  # name of run - change here if you want different variables
  if run['Field1_kind'] == 'DELTA':
    job['folder-run'] = 'DELTA-%.0e' % ( run['Field1_Int'][0])
  else:
    job['folder-run'] = '%s' %(run['Field1_kind'])
    job['folder-run'] += '-%.0e-%sfs-%seV-%sK' % ( run['Field1_Int'][0], run['Field1_Damp'][0], run['Field1_Freq'][0][0],job['temperature'])

  # writing input file
  os.system('mkdir -p %s/%s'%(folder,job['folder-run']))
  run.write('%s/%s/pulse.in'% (folder,job['folder-run']))

  # submission script
  oarsub.add_command('cd %s/%s; mpirun -hostfile \$OAR_NODEFILE %s -F pulse.in -J \'pulse,..//%s\' -C pulse -I \'../\''%(folder,job['folder-run'],yambo_rt,job['folder-col']) )
  oarsub.write('%s/%s/pulse.ll'%(folder,job['folder-run']))
  oarsub.run()
  print('running Dynamics without dissipation in folder: ' + str(job['folder-run']))

# Time-dependent with dissipation
if args.dissipation:

  # name of run - change here if you want different variables
  if run['Field1_kind'] == 'DELTA':
    job['folder-run'] = 'DELTA-D-%.0e' % ( run['Field1_Int'][0] )
  else:
    job['folder-run'] = '%s' %(run['Field1_kind'])
    job['folder-run'] += '-D-%seV-%sK-%sfs' % ( run['Field1_Freq'][0][0],job['temperature'],run['LifeExtrapSteps'][0][0])
    #job['folder-run'] += '-D-%.0e-%sfs-%seV-%sK' % ( run['Field1_Int'][0], run['Field1_Damp'][0], run['Field1_Freq'][0][0],job['temperature'])
  if job['DG'][0]:
    job['folder-run'] += '-DG'

  # writing input file
  os.system('mkdir -p %s/%s'%(folder,job['folder-run']))
  run.write('%s/%s/pulse.in'% (folder,job['folder-run']))

  # submission script
  if job['DG'][0]:
    oarsub.add_command('cd %s/%s; mpirun -hostfile \$OAR_NODEFILE %s -F pulse.in -J \'pulse,..//%s,..//%s,..//%s\' -C pulse -I \'../\''%(folder,job['folder-run'],yambo_rt,job['folder-col'],job['folder-gkkp'],job['DG'][1]) )
    print('Double Grid enabled')
  else:
    oarsub.add_command('cd %s/%s; mpirun -hostfile \$OAR_NODEFILE %s -F pulse.in -J \'pulse,..//%s,..//%s\' -C pulse -I \'../\''%(folder,job['folder-run'],yambo_rt,job['folder-col'],job['folder-gkkp']) )
  oarsub.write('%s/%s/pulse.ll'%(folder,job['folder-run']))
  oarsub.run()
  print('running Dynamics with dissipation in folder: ' + str(job['folder-run']))
