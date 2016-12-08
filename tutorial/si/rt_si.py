##############################################################################
#
# Author: Alejandro Molina-Sanchez
# Run real-time simulations with yambo
# 
# Warning: Real-time simulations requires several data folders for running
# properly. Before using this scripts compulsively is recommended
# to understand the different run levels.
#
# Instructions: 
# The dictionary 'job' is a personal choice to store useful instructions. This
# is the serial version but one can add options for running in parallel and 
# any other thing. Feel free to play with it.
# calculation : 'collision', 'negf', 'dissipation'
# folder-col  : collision data
# folder-run  : results (only work if collisions have been previously calculated)
# DG          : True or False if we use the double-grid (not yet implemented)
#
# Calculations are done inside the folder rt (feel free to rename it)
#
##############################################################################
#from __future__ import print_function
from yambopy import *
from qepy import *
import argparse
import os

# Select the run-level : 'collision', 'pump', 'dissipation'
parser = argparse.ArgumentParser(description='Example of real-time simulation')
parser.add_argument('-c' ,'--collisions')
parser.add_argument('-p' ,'--pump')
parser.add_argument('-d' ,'--dissipation')
args = parser.parse_args()

p2y      = 'p2y'
yambo    = 'yambo'
yambo_rt = 'yambo_rt'
ypp_rt   = 'ypp_rt'
ypp_ph   = 'ypp_ph'
folder   = 'rt'

job = dict()
job['folder-run']   = ''                # Optional additional job identifier
job['folder-col']   = 'col-hxc'              # Collisions folder
job['folder-gkkp']  = 'GKKP'                 # gkkp folder
job['DG']           = (False,'dg-4x4x4')     # Double-grid folder
job['temperature']  = 0.0                    # Temperature phonon bath

# check if the database is present
if not os.path.isdir('database'):
    os.mkdir('database')

#check if the nscf cycle is present
if os.path.isdir('nscf/si.save'):
    print('nscf calculation found!')
else:
    print('nscf calculation not found!')
    exit() 

#check if the SAVE folder is present
if not os.path.isdir('database/SAVE'):
    print('preparing yambo database')
    os.system('cd nscf/si.save; %s ;%s ; mv SAVE ../../database' % (p2y,yambo))

#check if the rt folder is present
if os.path.isdir('%s/SAVE'%folder):
  print('Symmetries for carrier dynamics ready') 
if not os.path.isdir('%s/SAVE'%folder):
  breaking_symmetries([1,0,0],folder=folder)

if args.collisions:
  print 'Collisions'  
  run = YamboIn('%s -r -e -v hsex'%yambo_rt,folder=folder)
elif args.pump:
  print 'Time-dependent with electric field'
  run = YamboIn('%s -q p'%yambo_rt,folder=folder)
elif args.dissipation:
  print 'Time-dependent with electric field and electron-phonon scattering'
  run = YamboIn('%s -s p -q p'%yambo_rt,folder=folder)
else:
  print 'Invalid calculation type'
  exit()

# Collision variables
if args.collisions:
  run['FFTGvecs']  = [5,'Ha']
  run['HARRLvcs']  = [5,'Ha']     # Hartree term: Equivalent to BSENGexx in the BSE run-level
  run['EXXRLvcs']  = [500,'mHa']  # Forck term:   Equivalent to BSENGBlk in the BSE run-level
  run['CORRLvcs']  = [500,'mHa']  # Correlation term: Not appearing in BSE. 
  run['NGsBlkXs']  = [500,'mHa']
  run['BndsRnXs' ] = [1,30]       # Static screening
  run['COLLBands'] = [2,7]        # Electron-Hole states
  run['HXC_Potential']  = 'HARTREE+SEX' 
  run.write('%s/03_COLLISION'%folder)

# Common time-dependent variable
if args.pump or args.dissipation:
  run['RTBands']    = [2,7]        # Electron-Hole states
  run['GfnQP_Wv']   = [0.05,0.00,0.00]    # Constant damping valence
  run['GfnQP_Wc']   = [0.05,0.00,0.00]    # Constant damping conduction
  run['GfnQP_E']    = [0.00, 1.00, 1.00]  # [EXTQP BSK BSS] E parameters  (c/v) eV|adim|adim
  run['HXC_Potential']  = 'HARTREE+SEX' 
  # Time-propagation 
  run['RTstep']     = [   5.0,'as']
  run['NETime']     = [   1.0,'ps']
  run['Integrator'] = "RK2 RWA"
  run['IOtime']     = [ [10.000, 10.000, 1.000], 'fs' ] 
  # Pump Pulse
  run['Field1_Int']       = [ 1E4 , 'kWLm2']    # Intensity pulse
  run['Field1_Dir']       = [1.0,0.0,0.0]  # Polarization pulse
  run['Field1_Dir_circ']  = [0.0,1.0,0.0]  # Polarization pulse
  run['Field1_pol']       = "linear"       # Polarization type (linear or circular) 
  run['Field1_kind']      = "QSSIN"        # [RT Field1] Kind(SIN|RES|ANTIRES|GAUSS|DELTA|QSSIN)
  run['Field1_Damp']      = [ 50,'fs']
  run['Field1_Freq']      = [[2.3,2.3],'eV']

# Pumping with finite pulse and electron-phonon dissipation
if args.pump or args.dissipation:
# Interpolation 
  run['LifeInterpKIND']  = 'FLAT'
  run['LifeInterpSteps'] = [ [1.0,1.0], 'fs' ] 

# Submission in serial

# Collisions

if args.collisions:
  print('running yambo-collision')
  os.system('cd %s; %s -F 03_COLLISION -J %s'%(folder,yambo_rt,job['folder-col']))

# Dynamics without dissipation and without/with Double Grid 

if args.pump:
  print('running pumping with finite pulse')
  run.write('%s/04_PUMP'%folder)
  if run['Field1_kind'] == 'DELTA':
    jobname = '%s%s-%.0e' % (job['folder-run'], run['Field1_kind'], run['Field1_Int'][0] )
  else:
    jobname = '%s%s-%.0e-%sfs-%seV-%sK' % ( job['folder-run'], run['Field1_kind'],run['Field1_Int'][0], run['Field1_Damp'][0], run['Field1_Freq'][0][0],job['temperature'] )
  print('running NEGF in folder: %s' % job['folder-run'])
  if job['DG'][0]:
    print('with Double Grid from folder %s'%job['DG'][1])
  else:
    print 'cd %s ; %s -F 04_PUMP -J \'%s,%s\' -C %s'%(folder,yambo_rt,jobname,job['folder-col'],jobname)
    os.system ('cd %s; %s -F 04_PUMP -J \'%s,%s\' -C %s'%(folder,yambo_rt,jobname,job['folder-col'],jobname))

# Time-dependent with a pulse and dissipation and without/with Double Grid 

if args.dissipation:
  run.write('%s/05_DISS'%folder)
  print('running pumping with finite pulse and with electron-phonon scattering')
  print('this run level needs the GKKP folder to run')
  job['folder-run'] += 'dneq'
  if run['Field1_kind'] == 'DELTA':
    jobname = '%s%s-%.0e' % (job['folder-run'], run['Field1_kind'], run['Field1_Int'][0] )
  else:
    jobname = '%s%s-%.0e-%sfs-%seV-%sK' % ( job['folder-run'], run['Field1_kind'],run['Field1_Int'][0], run['Field1_Damp'][0], run['Field1_Freq'][0][0],job['temperature'] )
  if job['DG'][0]:
    print('with Double Grid from folder %s'%job['DG'][1])
    print('%s -F 06_DISS -J \'%s,%s,%s\''%(yambo_rt,jobname,job['folder-col'],job['DG'][1]))
  else:
    print 'cd %s; %s -F 05_DISS -J \'%s,%s,%s\' -C %s'%(folder,yambo_rt,jobname,job['folder-col'],job['folder-gkkp'],jobname)
    os.system( 'cd %s; %s -F 05_DISS -J \'%s,%s,%s\' -C %s'%(folder,yambo_rt,jobname,job['folder-col'],job['folder-gkkp'],jobname) )
