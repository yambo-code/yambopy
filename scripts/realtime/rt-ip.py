#########################################################
#  YAM(BO)PY(THON) Library
#
#  Generation of Yambo input files using python
#
#  Authors: A Molina-Sanchez, HPC Miranda
#
#  January 2016
#########################################################
# Calculation of spectra in the Independant Particle approximation
# Author: Alejandro Molina-Sanchez
#########################################################
# Collection of real-time scripts:
# cohsex-rt.py | bse-rt.py | merging-rt.py | ip-rt.py
#########################################################

#########################################################
from __future__ import print_function
from yambopy     import *
from schedulerpy import *

############## SETTINGS ##############

yambo_module = 'yambo/master-intel'
yambo_rt     = 'yambo_rt'
ypp_rt       = 'ypp_rt'

#source       = 'rt-D-1.94eV-0K-2.0fs-DG'
source       = 'QSSIN-1e+03-70.0fs-1.94eV-0K' # no dissipation
folder_kerr  = 'kerr-24x24'
folder_rt    = 'rt-24x24'

ip_nodes =  1
ip_cores =  12

time_probe = range(0,610,150)
#time_probe=(0,)
print(time_probe)


dir_pump   = '../%s/%s/' % (folder_rt,source)
link_pump  = '../%s/%s/' % (folder_rt,source)
dir_inputs = 'inputs'
os.system('cd %s; mkdir -p %s'%(folder_kerr,dir_inputs))


ip = YamboIn('%s -o b'%yambo_rt,folder=folder_kerr) # NEQ COHSEX
# Common variables
ip['DBsIOoff'] = 'DIP'
ip['FFTGvecs'] = [ 20   , 'Ha'  ]
ip['BSKmod']   = 'IP'
ip['BSENGexx'] = [ 20   , 'Ha'  ]
ip['BSENGBlk'] = [ 1000 , 'mHa' ]
ip['NGsBlkXs'] = [ 1000 , 'mHa' ]
ip['BndsRnXs'] = [1 ,60]
ip['BSEBands'] = [25,28]
#paralelization
ip['X_all_q_ROLEs'] = 'q.k.c.v'
ip['X_all_q_CPU']   = '1.12.1.1'
ip['BS_ROLEs']      = 'k.eh.t'
ip['BS_CPU']        = '12.1.1'
ip['BSKmod']        = "IP"
ip['BEnRange'] =  [ [0.0  , 4.0] , 'eV'  ]    # Energy range spectra
ip['BLongDir'] =  [ 1.0, 0.0, 0.0 ]           # [BSS] [cc] Electric Field
ip['BDmRange'] =  [ [0.100 , 0.100] , 'eV'  ]   # Width
ip['BEnSteps'] =    1000                      # Energy steps
#ip.arguments.append('EvalKerr')
#ip.arguments.append('ForceEqTrans')


# Submission of the jobs for each time

print('Running IP calculations...')
print('RT source calculation: %s \n' % source)

for time in time_probe:
    print('Time of carriers database %d' % time)
    ip['XfnRTdb'] = 'f @ %d fs < ../%s/%s/ndb.RT_carriers' % ( time , folder_rt, link_pump )
    ip['KfnRTdb'] = 'f @ %d fs < ../%s/%s/ndb.RT_carriers' % ( time , folder_rt, link_pump )
    ip['RfnRTdb'] = 'f @ %d fs < ../%s/%s/ndb.RT_carriers' % ( time , folder_rt, link_pump )
    nameip        = 'IP-%s-t%d'              % ( source, time )
    print(nameip)
    ip.write('%s/%s/%s.in' %(folder_kerr, dir_inputs, nameip))
    yambo = oarsub(core=ip_cores,dependent=0,name='ip-rt',walltime="02:00:00")
    yambo.add_command('module load %s'%yambo_module)
    yambo.add_command('export OMP_NUM_THREADS=1')
    yambo.add_command('cd %s; mpirun -machinefile \$OAR_NODEFILE %s -F %s/%s.in -J %s -C %s'%(folder_kerr, yambo_rt, dir_inputs, nameip, nameip, nameip))
    yambo.run()
    yambo.clean()
