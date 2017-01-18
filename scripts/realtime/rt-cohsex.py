#########################################################
#  YAM(BO)PY(THON) Library
#
#  Generation of Yambo input files using python
#
#  Authors: A Molina-Sanchez, HPC Miranda
#
#  January 2016
#########################################################
# Calculation of COHSEX corrections for different times
# in real-time calculations
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

folder_rt    = 'rt-24x24'
#source       = 'rt-D-1.94eV-0K-0.1fs-DG'
source       = '0.05fs'
#source       = 'QSSIN-1e+03-70.0fs-1.94eV-0K' # no dissipation
folder_kerr  = 'kerr-24x24'
# Folder folder_cohsex in other -rt scripts is for the ndb.QP @ t=0.

CSRTmode = 'XG'  #  X: Screening, G: GFs

nodes =  4
cores =  12

time_probe = range(0,610,150)
#time_probe=(0,)
print(time_probe)

#######################################


dir_pump   = '../%s/%s/' % (folder_rt,source)
link_pump  = '../%s/%s/' % (folder_rt,source)
dir_inputs = 'inputs'
os.system('cd %s; mkdir -p %s'%(folder_kerr,dir_inputs))


cs = YamboIn('%s -r -p c -g n'%yambo_rt,folder=folder_kerr) # NEQ COHSEX
#coulomb cutoff
cs['DBsIOoff'] = 'DIP'
cs['RandQpts'] =  1000000
cs['RandGvec'] =  [ 1, 'RL' ]
cs['CUTGeo']   =  "box z"
cs['CUTBox']   =  [ 0, 0, 38]
# Common variables
cs['FFTGvecs'] = [ 20   , 'Ha'  ]
cs['EXXRLvcs'] = [ 20   , 'Ha'  ]
cs['NGsBlkXs'] = [ 1000 , 'mHa' ]
cs['BndsRnXs'] = [1,60]
cs['QPkrange'] = [1,576,25,28]
#paralelization
cs['X_all_q_ROLEs'] = 'q.k.c.v'
cs['X_all_q_CPU']   = '1.24.2.1'
cs['SE_ROLEs']      = 'q.qp.b'
cs['SE_CPU']        = '1.48.1'
# Flags
cs.arguments.append('ExtendOut')


# Submission of the jobs for each time
print('Running COHSEX calculation...')
print('RT source calculation: %s \n' % source)

for time in time_probe:
    print('Time of carriers database %d' % time)
    if CSRTmode == 'X':
      cs['XfnRTdb'] = 'f @ %d fs < ./%s/ndb.RT_carriers' % ( time , link_pump )
    elif CSRTmode == 'G':
      cs['GfnRTdb'] = 'f @ %d fs < ./%s/ndb.RT_carriers' % ( time , link_pump )
    elif CSRTmode == 'XG':
      cs['XfnRTdb'] = 'f @ %d fs < ./%s/ndb.RT_carriers' % ( time , link_pump )
      cs['GfnRTdb'] = 'f @ %d fs < ./%s/ndb.RT_carriers' % ( time , link_pump )
    else:
      print('Error in the RT run level')
      exit()
    namecs        = 'C-%s-%s-t%d'              % ( CSRTmode, source, time )
    print(namecs)
    cs.write('%s/%s/%s.in' %(folder_kerr, dir_inputs, namecs))
    yambo = oarsub(nodes=nodes,core=cores,dependent=4078486,name='cohsex',walltime="10:00:00")
    yambo.add_command('module load %s'%yambo_module)
    yambo.add_command('export OMP_NUM_THREADS=1')
    #yambo.add_command('mpirun -npernode 12 -x OMP_NUM_THREADS -x PATH -x LD_LIBRARY_PATH -hostfile \$OAR_NODEFILE %s -F %s/%s.in -J %s -C %s'%(yambo_rt,dir_inputs,namecs,namecs,namecs))
    #yambo.add_command('mpirun -machinefile \$OAR_NODEFILE %s -F %s/%s.in -J %s -C %s'%(yambo_rt,dir_inputs,namecs,namecs,namecs))
    yambo.add_command('mpirun -hostfile \$OAR_NODEFILE %s -F %s/%s.in -J %s -C %s'%(yambo_rt,dir_inputs,namecs,namecs,namecs))
    yambo.write('%s/%s.ll' % (folder_kerr, namecs))
    os.system('cd %s; sh %s.ll'% (folder_kerr, namecs))
    yambo.clean()
