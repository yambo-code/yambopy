#
# Author: Alejandro Molina-Sanchez
# Run real-time simulations with yambo
# We put this script in the same folder than the SAVE directory 
#
#from __future__ import print_function
from yambopy.inputfile import *
from pwpy.inputfile import *
from pwpy.outputxml import *

yambo    = 'yambo'
yambo_rt = 'yambo_rt'
ypp_rt   = 'ypp_rt'

'''
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
    log_p2y   = shell.run(['p2y'],  cwd="nscf/si.save")
    print(log_p2y.output,file=open("p2y.log","w"))
    log_yambo = shell.run(['yambo'],cwd="nscf/si.save")
    print(log_yambo.output,file=open("yambo.log","w"))
    shell.run('mv nscf/si.save/SAVE database'.split())
'''

#check if the FixSymm folder is present
#if not os.path.isdir('FixSymm'):
#  print('breaking symmetries')
#  shell.run('mv FixSymm/SAVE SAVE'.split())
   
def realtime(name_folder):
  a, b, c, d = YamboIn('%s -r -e -v c -V all'%yambo_rt,folder='FixSymm'),YamboIn('%s -q p -v c -V all'%yambo_rt,folder='FixSymm'),YamboIn('%s -q p -v c -V all'%yambo_rt,folder='FixSymm'),YamboIn('%s -s p -q p -v c -V all'%yambo_rt,folder='FixSymm')
  os.system('rm FixSymm/yambo.in')
  return a, b, c, d

job = dict()

job['calculation'] = 'pump'
job['folder-col']  = 'COLLISION'
job['folder-run']  = 'NEGF-1e5-80fs'
job['DG']          = (False,'DG-60x60')
job['nodes']       = 1
job['cores']       = 1

col, tdsex, pump, diss = realtime('.')

# System Common variables
col['FFTGvecs']  = [2,'Ha']
col['EXXRLvcs']  = [2,'Ha']
col['SCBands']   = [[4,5],'']

tdsex['FFTGvecs'], tdsex['EXXRLvcs'], tdsex['SCBands'] = col['FFTGvecs'], col['EXXRLvcs'], col['SCBands']
pump['FFTGvecs'] , pump['EXXRLvcs'] ,  pump['SCBands'] = col['FFTGvecs'], col['EXXRLvcs'], col['SCBands']
diss['FFTGvecs'] , diss['EXXRLvcs'] ,  diss['SCBands'] = col['FFTGvecs'], col['EXXRLvcs'], col['SCBands']

# Collision
col['NGsBlkXp']  = [ 100,'mHa']
col['BndsRnXs' ] = [[1,20],'']
#col['CUTGeo']    = 'box z'
#col['CUTBox']    = [[0,0,38], '']
col.write('FixSymm/03_COLLISION')

# Time-dependent COHSEX -- DELTA PULSE
# System
tdsex['GfnQP_Wv']   = [ [0.04,0.00,0.00],'' ]
tdsex['GfnQP_Wc']   = [ [0.04,0.00,0.00],'' ]
tdsex['GfnQPdb']    = 'none' 
# Time-propagation 
tdsex['RTstep']     = [0.01,'fs']
tdsex['NETime']     = [  80,'fs']
tdsex['Integrator'] = "RK2 RWA"
tdsex['IOtime']     = [ [0.05, 0.05, 0.01], 'fs' ] 
# Pump Pulse
tdsex['Probe_Int']  = [ 1E5 , 'kWLm2']
tdsex['Probe_kind'] = "DELTA"
tdsex['Probe_Damp'] = [0,'fs']
tdsex['Probe_Freq'] = [[0.0,0.0],'eV']
tdsex.write('FixSymm/04_TDSEX')

# Pumping with finite pulse
# System
pump['GfnQP_Wv']   = [ [0.04,0.00,0.00],'' ]
pump['GfnQP_Wc']   = [ [0.04,0.00,0.00],'' ]
pump['GfnQPdb']    = 'none' 
# Time-propagation 
pump['RTstep']     = [0.01,'fs']
pump['NETime']     = [   1,'ps']
pump['Integrator'] = "RK2 RWA"
pump['IOtime']     = [ [0.05, 0.05, 0.01], 'fs' ] 
# Pump Pulse
pump['Probe_Int']  = [ 1E5 , 'kWLm2']
pump['Probe_kind'] = "QSSIN"
pump['Probe_Damp'] = [100,'fs']
pump['Probe_Freq'] = [[2.5,0.0],'eV']
pump.write('FixSymm/05_NEGF')

# Pumping with finite pulse
# System
diss['GfnQP_Wv']   = [ [0.04,0.00,0.00],'' ]
diss['GfnQP_Wc']   = [ [0.04,0.00,0.00],'' ]
diss['GfnQPdb']    = 'none' 
# Time-propagation 
diss['RTstep']     = [0.01,'fs']
diss['NETime']     = [   1,'ps']
diss['Integrator'] = "RK2 RWA"
diss['IOtime']     = [ [0.05, 0.05, 0.01], 'fs' ] 
diss['LifeInterpKIND']  = 'FLAT'
diss['LifeInterpSteps'] = [ [10.0,2.5], 'fs' ] 
# Pump Pulse
diss['Probe_Int']  = [ 1E5 , 'kWLm2']
diss['Probe_kind'] = "QSSIN"
diss['Probe_Damp'] = [100,'fs']
diss['Probe_Freq'] = [[2.5,0.0],'eV']
diss.write('FixSymm/06_DISS')

# Run 

# Collisions
if job['calculation'] == 'collision':
  print('running yambo-collision')
  #log_yambo = shell.run(('%s -F 03_COLLISION -J '%yambo,'COLLISION').split(),cwd='gw')
  print ('cd FixSymm; %s -F 03_COLLISION -J %s'%(yambo_rt,job['folder-col']))
  os.system('cd FixSymm; %s -F 03_COLLISION -J %s'%(yambo_rt,job['folder-col']))
# Time-dependent without/with Double Grid 
if job['calculation'] == 'tdsex':
  print('running TD-COHSEX')
  if job['DG'][0]:
    print('with Double Grid from folder %s'%job['DG'][1])
    print('cd FixSymm; %s -F 04_TDSEX -J \'%s,%s,%s\''%(yambo_rt,job['folder-run'],job['folder-col'],job['DG'][1]))
    os.system('cd FixSymm; %s -F 04_TDSEX -J \'%s,%s,%s\''%(yambo_rt,job['folder-run'],job['folder-col'],job['DG'][1]))
  else:
    print ('cd FixSymm; %s -F 04_TDSEX -J \'%s,%s\''%(yambo_rt,job['folder-run'],job['folder-col']))
    os.system ('cd FixSymm; %s -F 04_TDSEX -J \'%s,%s\''%(yambo_rt,job['folder-run'],job['folder-col']))
if job['calculation'] == 'pump':
  print('running pumping with finite pulse')
  #print('Probe_Damp %s %s'%str(pump['Probe_Damp']) )
  if job['DG'][0]:
    print('with Double Grid from folder %s'%job['DG'][1])
    print('%s -F 05_NEGF -J \'%s,%s,%s\''%(yambo_rt,job['folder-run'],job['folder-col'],job['DG'][1]))
  else:
    print ('%s -F 05_NEGF -J \'%s,%s\''%(yambo_rt,job['folder-run'],job['folder-col']))
    os.system ('cd FixSymm; %s -F 05_NEGF -J \'%s,%s\''%(yambo_rt,job['folder-run'],job['folder-col']))
if job['calculation'] == 'dissipation':
  print('running pumping with finite pulse')
  #print('Probe_Damp %s %s'%str(pump['Probe_Damp']) )
  if job['DG'][0]:
    print('with Double Grid from folder %s'%job['DG'][1])
    print('%s -F 06_DISS -J \'%s,%s,%s\''%(yambo_rt,job['folder-run'],job['folder-col'],job['DG'][1]))
  else:
    print ('%s -F 06_DISS -J \'%s,%s\''%(yambo_rt,job['folder-run'],job['folder-col']))
