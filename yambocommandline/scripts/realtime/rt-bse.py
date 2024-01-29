#########################################################
#  YAM(BO)PY(THON) Library
#
#  Generation of Yambo input files using python
#
#  Authors: A Molina-Sanchez, HPC Miranda
#
#  January 2016
#########################################################
# Bethe-Salpeter spectra calculations at given times
# of a Real-Time simulation
# Author: Alejandro Molina-Sanchez
#########################################################

#########################################################
from __future__ import print_function
from builtins import range
from yambopy     import *
from schedulerpy import *

############## SETTINGS ##############

yambo_module = 'yambo/intel-4.1'
yambo_rt     = 'yambo_rt'

source       = 'QSSIN-1e+03-70.0fs-2.0eV-0K'
folder_rt    = 'rt-24x24'
folder_gw    = 'gw-24x24' # Where the ndb.QP is located (see l. 93)

BSRTmode = 'XRK' #  X: Screening, R: Residuals, K: Kernel
CSRTmode = 'XG'  #  X: Screening, G: GFs
QPdata   = 'N'   #  E: Equilibrium QPs, N: Non-Equilibrium QPs, L: Scissor operator

bs_nodes =  1
bs_cores =  12

time_probe = list(range(0,610,150))
#time_probe=(0,)
print(time_probe)

##########################################################
# Scissors Operator
# GW for 6x6  12x12   18x18   24x24   30x30   42x42  51x51
gw = [2.0, 1.41334, 1.16886, 1.04588, 0.98218, 0.9214, 0.8995]
##########################################################


bs = YamboIn('%s  -r -b -o b -k sex -y d -Q'%yambo_rt,folder=folder_rt)  # BS
#bs = YamboIn('%s -o g -Q'%yambo_rt,vim=False)  # BS
bs['DBsIOoff'] = ''
bs['BSKmod']   = "SEX"                # [BSE] IP/Hartree/HF/ALDA/SEX/BSfxc
bs['X_all_q_ROLEs'] = 'q.k.c.v'
bs['X_all_q_CPU']   = '1.%d.2.1'%(bs_nodes*bs_cores/2)  # bs_nodes bs_cores
bs['BS_ROLEs']      = 'k.eh.t'
bs['BS_CPU']        = '%d.1.1'%(bs_nodes*bs_cores)
bs['RandQpts'] =  1000000
bs['RandGvec'] =  [ 1, 'RL' ]
bs['CUTGeo']   =  "box z"
bs['CUTBox']   =  [ 0.00, 0.00, 38.0]
# Common variables
#bs['KfnQP_E']  = [gw[3], 1.00, 1.00]       # [EXTQP BSK BSS] E parameters  (c/v) eV|adim|adim
bs['KfnQP_E']  = [0.0, 1.00, 1.00]       # [EXTQP BSK BSS] E parameters  (c/v) eV|adim|adim
bs['FFTGvecs'] = [ 20   , 'Ha'  ]
bs['BSENGexx'] = [ 20   , 'Ha'  ]
bs['BSENGBlk'] = [ 1000 , 'mHa' ]
bs['NGsBlkXs'] = [ 1000 , 'mHa' ]
bs['BndsRnXs'] = [1 ,60]
bs['BSEBands'] = [25,28]

bs['Gauge']    =  "length"
bs['BSEmod']   =  "causal"
bs['BSSmod']   =  "d"
bs['BEnRange'] =  [ [0.0  , 4.0] , 'eV'  ]    # Energy range spectra
bs['BLongDir'] =  [ 1.0, 0.0, 0.0 ]           # [BSS] [cc] Electric Field
bs['BDmRange'] =  [ [0.100 , 0.100] , 'eV'  ]   # Width
bs['BEnSteps'] =    1000                      # Energy steps
bs.arguments.append('WRbsWF')
#bs.arguments.append('ForceEqTrans')

# Submission of the jobs for each time

print('Running BSE calculations...')
print('RT source calculation: %s \n' % source)

for time in time_probe:
    print('Time of carriers database %d' % time)
    if 'X' in BSRTmode:
      bs['XfnRTdb'] = 'f @ %d fs < ./pulse/ndb.RT_carriers' % ( time )
    if 'R' in BSRTmode:
      bs['RfnRTdb'] = 'f @ %d fs < ./pulse/ndb.RT_carriers' % ( time )
    if 'K' in BSRTmode:
      bs['KfnRTdb'] = 'f @ %d fs < ./pulse/ndb.RT_carriers' % ( time )
    if 'E' in QPdata:
      bs['KfnQPdb'] = 'E < ../../%s/ndb.QP' % folder_gw         # GW database
      namebs      = 'B-%s-%s-t%d' % ( BSRTmode, QPdata, time )
    if 'N' in QPdata:
      namebs      = 'B-%s-%s-t%d' % ( BSRTmode, CSRTmode, time )
      name_merged = 'M-%s-t%d'    % ( CSRTmode, time )
      bs['KfnQPdb'] = 'E < ./%s/ndb.QP' % name_merged  # GW + NEQ_COHSEX - EQ_COHSEX
    if 'L' in QPdata:
      namebs      = 'B-L-t%d' % ( time )
      print('Use LDA+Scissor')
    bs.write('%s/%s/%s.in' %(folder_rt, source, namebs))
    yambo = oarsub(nodes=bs_nodes,core=bs_cores,dependent=0,name='bse-rt',walltime="10:00:00")
    yambo.add_command('module load %s'%yambo_module)
    yambo.add_command('export OMP_NUM_THREADS=1')
    yambo.add_command('cd %s/%s ; mpirun -hostfile \$OAR_NODEFILE %s -F %s.in -J %s -C %s -I \'../\''%(folder_rt,source,yambo_rt,namebs,namebs,namebs))
    yambo.write('%s/%s/%s.ll' % (folder_rt, source, namebs))
    yambo.run()

