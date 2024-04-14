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
from builtins import range
from yambopy     import *
from schedulerpy import *


############## SETTINGS ##############

yambo_module = 'yambo/intel-4.1'
yambo_rt     = 'yambo_rt'

folder_rt    = 'rt-6x6'
source       = 'QSSIN-1e+03-70.0fs-2.0eV-0K'

CSRTmode = 'XG'  #  X: Screening, G: GFs

nodes =  1
cores =  12

time_probe = list(range(0,610,150))
#time_probe=(0,)
print(time_probe)


#######################################

cs = YamboIn('%s -r -p c -g n'%yambo_rt,folder=folder_rt) # NEQ COHSEX
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
cs['X_all_q_CPU']   = '1.%d.2.1'%(nodes*cores/2)
cs['SE_ROLEs']      = 'q.qp.b'
cs['SE_CPU']        = '1.%d.1'%(nodes*cores)
# Flags
cs.arguments.append('ExtendOut')


# Submission of the jobs for each time
print('Running COHSEX calculation...')
print('RT source calculation: %s \n' % source)

for time in time_probe:
    print('Time of carriers database %d' % time)
    if CSRTmode == 'X':
      cs['XfnRTdb'] = 'f @ %d fs < ./pulse/ndb.RT_carriers' % ( time )
    elif CSRTmode == 'G':
      cs['GfnRTdb'] = 'f @ %d fs < ./pulse/ndb.RT_carriers' % ( time )
    elif CSRTmode == 'XG':
      cs['XfnRTdb'] = 'f @ %d fs < ./pulse/ndb.RT_carriers' % ( time )
      cs['GfnRTdb'] = 'f @ %d fs < ./pulse/ndb.RT_carriers' % ( time )
    else:
      print('Error in the RT run level')
      exit()
    namecs        = 'C-%s-t%d'              % ( CSRTmode, time )
    print(namecs)
    cs.write('%s/%s/%s.in' %(folder_rt, source, namecs))
    yambo = oarsub(nodes=nodes,core=cores,dependent=0,name='cohsex',walltime="10:00:00")
    yambo.add_command('module load %s'%yambo_module)
    yambo.add_command('export OMP_NUM_THREADS=1')
    yambo.add_command('cd %s/%s ; mpirun -hostfile \$OAR_NODEFILE %s -F %s.in -J %s -C %s -I \'../\''%(folder_rt,source,yambo_rt,namecs,namecs,namecs))
    yambo.write('%s/%s/%s.ll' % (folder_rt, source, namecs))
    yambo.run()
