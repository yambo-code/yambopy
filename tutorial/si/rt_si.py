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
# calculation : 'collision', 'tdsex', 'negf', 'dissipation'
# folder-col  : collision data
# folder-run  : results (only work if collisions have been previously calculated)
# DG          : True or False if we use the double-grid
#
# Calculations are done inside the folder rt (feel free to rename it)
#
##############################################################################
#from __future__ import print_function
from yambopy import *
from qepy import *
import os

yambo    = 'yambo'
yambo_rt = 'yambo_rt'
ypp_rt   = 'ypp_rt'
ypp_ph   = 'ypp_ph'
folder   = 'rt'

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
    os.system('cd nscf/si.save; p2y')
    os.system('cd nscf/si.save; yambo')
    os.system('mv nscf/si.save/SAVE database')

#check if the rt folder is present
if os.path.isdir('%s/SAVE'%folder):
  print('symmetries for carrier dynamics ready') 
if not os.path.isdir('%s/SAVE'%folder):
  breaking_symmetries([1,0,0],folder=folder)

# you can select the calculation among : 'collision', 'tdsex', 'pump', 'dissipation'

job = dict()
job['calculation']  = 'collision'
job['folder-run']   = 'Diss'
job['folder-col']   = 'COLLISION'
job['folder-gkkp']  = 'GKKP'
job['DG']           = (False,'DG-60x60')

if job['calculation']=='collision':
  print 'Collision'  
  run = YamboIn('%s -r -e -v c -V all'%yambo_rt,folder=folder)
elif job['calculation']=='tdsex':
  run = YamboIn('%s -q p -v c -V all'%yambo_rt,folder=folder)
elif job['calculation']=='pump':
  run = YamboIn('%s -q p -v c -V all'%yambo_rt,folder=folder)
elif job['calculation']=='dissipation':
  run = YamboIn('%s -s p -q p -v c -V all'%yambo_rt,folder=folder)
else:
  print 'Invalid calculation type'
  exit()

# System Common variables
run['FFTGvecs']  = [5,'Ha']
run['EXXRLvcs']  = [5,'Ha']

# Collision variables
if job['calculation']=='collision':
  run['NGsBlkXs']  = [ 100,'mHa']
  run['BndsRnXs' ] = [1,30]
  run['COLBands']  = [2,7]
  run.write('%s/03_COLLISION'%folder)

# Common time-dependent variable
if job['calculation']=='tdsex' or 'pump' or 'dissipation':
  run['RTBands']    = [2,7]
  run['GfnQP_Wv']   = [0.01,0.00,0.00]
  run['GfnQP_Wc']   = [0.01,0.00,0.00]
  run['GfnQPdb']    = 'none' 
  run['Potential']  = 'COHSEX' 
  # Time-propagation 
  run['RTstep']     = [0.100,'fs']
  run['NETime']     = [ 300,'fs']
  run['Integrator'] = "RK2 RWA"
  run['IOtime']     = [ [2.000, 2.000, 2.000], 'fs' ] 
  # Pump Pulse
  run['Field1_Int']       = [ 1E6 , 'kWLm2']    # Intensity pulse
  run['Field1_Dir']       = [1.0,0.0,0.0]  # Polarization pulse
  run['Field1_Dir_circ']  = [0.0,1.0,0.0]  # Polarization pulse
  run['Field1_pol']       = "linear"            # Polarization type (linear or circular) 

# Time-dependent COHSEX -- DELTA PULSE
if job['calculation']=='tdsex':
  run['Field1_kind'] = "DELTA"
  run.write('%s/04_TDSEX'%folder)

# Pumping with finite pulse
if job['calculation']=='pump' or 'dissipation':
  run['Field1_kind'] = "QSSIN"
  run['Field1_Damp'] = [  50,'fs']
  run['Field1_Freq'] = [[2.5,0.0],'eV']

if job['calculation']=='pump':
  run.write('%s/05_NEGF'%folder)

# Pumping with finite pulse and electron-phonon dissipation
if job['calculation']=='dissipation':
# Interpolation 
  run['LifeInterpKIND']  = 'FLAT'
  run['LifeInterpSteps'] = [ [4.0,1.0], 'fs' ] 
  run.write('%s/06_DISS'%folder)

# Run

# Collisions

if job['calculation'] == 'collision':
  print('running yambo-collision')
  os.system('cd %s; %s -F 03_COLLISION -J %s'%(folder,yambo_rt,job['folder-col']))

# Time-dependent without/with Double Grid 

if job['calculation'] == 'tdsex':
  job['folder-run'] += 'tds-int-%s' % ( str(int(run['Field1_Int'][0])) )
  print('running TD-COHSEX in folder: ' + str(job['folder-run']))
  if not os.path.isdir('%s/'%folder+job['folder-col']):
    print 'Collisions not found'
    exit()
  if job['DG'][0]:
    print('with Double Grid from folder %s'%job['DG'][1])
    os.system('cd %s; %s -F 04_TDSEX -J \'%s,%s,%s\' -C %s'%(folder,yambo_rt,job['folder-run'],job['folder-col'],job['DG'][1],job['folder-run']))
  else:
    os.system ('cd %s; %s -F 04_TDSEX -J \'%s,%s\' -C %s'%(folder,yambo_rt,job['folder-run'],job['folder-col'],job['folder-run']))

# Time-dependent with a pulse and without/with Double Grid 

if job['calculation'] == 'pump':
  print('running pumping with finite pulse')
  job['folder-run'] += 'negf-int-%s-damp-%sfs-freq-%seV' % ( str(int(run['Field1_Int'][0])), str(run['Field1_Damp'][0]),run['Field1_Freq'][0][0] )
  print('running NEGF in folder: ' + str(job['folder-run']))
  if not os.path.isdir('%s/'%folder+job['folder-col']):
    print 'Collisions not found'
    exit()
  if job['DG'][0]:
    print('with Double Grid from folder %s'%job['DG'][1])
  else:
    os.system ('cd %s; %s -F 05_NEGF -J \'%s,%s\' -C %s'%(folder,yambo_rt,job['folder-run'],job['folder-col'],job['folder-run']))

# Time-dependent with a pulse and dissipation and without/with Double Grid 

if job['calculation'] == 'dissipation':
  print('running pumping with finite pulse and with electron-phonon scattering')
  print('this run level needs the GKKP folder to run')
  if not os.path.isdir('%s/'%folder+job['folder-col']):
    print 'Collisions not found'
    exit()
  if os.path.isdir('%s/GKKP'%folder):
    print('Gkkp files exists')
    if job['DG'][0]:
      print('with Double Grid from folder %s'%job['DG'][1])
      print('%s -F 06_DISS -J \'%s,%s,%s\''%(yambo_rt,job['folder-run'],job['folder-col'],job['DG'][1]))
    else:
      os.system ('cd %s; %s -F 06_DISS -J \'%s,%s,%s\' -C %s'%(folder,yambo_rt,job['folder-run'],job['folder-col'],job['folder-gkkp'],job['folder-run']))
  else:
    print('No gkkp files. Calculation stop')
    print('You may run gkkp_si.py')
    exit()
