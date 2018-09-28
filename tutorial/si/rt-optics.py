from __future__ import print_function
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
# Author: Alejandro Molina-Sanchez
# Run real-time simulations with yambo
#########################################################

#########################################################
#from __future__ import print_function
from yambopy.inputfile import *
from pwpy.inputfile import *
from pwpy.outputxml import *
from oarfile import *

yambo_rt     = 'yambo_rt'
ypp_rt       = 'ypp_rt'
yambo_module = 'devel-5420'

Set      = 'M'   #  C: Cohsex, M: Merging, B: BSE
BSRTmode = 'XRK' #  X: Screening, R: Residuals, K: Kernel
CSRTmode = 'XG'  #  X: Screening, G: GFs
QPdata   = 'N'   #  E: Equilibrium QPs, N: Non-Equilibrium QPs

cs_nodes =  4 
cs_cores = 12 
bs_nodes =  4 
bs_cores = 12 

time_probe = [0,500,1000] 
source     = 'Diss'
exc_pump   = 1.71
damp_pump  = 70

dir_pump   = '../RT/FixSymm/%s-1e+05-%dfs-%.2feV' % ( source, damp_pump, exc_pump)
link_pump  = 'Freq%.2feV-Damp%dfs-%s'            %  ( exc_pump, damp_pump, source)

# Check RT simulations exists

print('Checking link...')
if not os.path.isdir(link_pump):
  os.system('ln -s ' + dir_pump + ' ' + link_pump)

# A. Cohsex Step

# Select the Yambo Run-level
dir_inputs = 'inputs'
os.system('mkdir -p %s'%dir_inputs)
cs = YamboIn('%s -r -p c -g n'%yambo_rt,vim=False) # NEQ COHSEX
#coulomb cutoff
cs['DBsIOoff'] = 'DIP'
cs['RandQpts'] =  1000000
cs['RandGvec'] =  [ 1, 'RL' ]
cs['CUTGeo']   =  "box z"
cs['CUTBox']   =  [ 0, 0, 38.0]
# Common variables
cs['FFTGvecs'] = [ 15   , 'Ha'  ]
cs['EXXRLvcs'] = [ 15   , 'Ha'  ]
cs['NGsBlkXs'] = [ 1200 , 'mHa' ]
cs['BndsRnXs'] = [1,70]
cs['QPkrange'] = [1,91,25,28]
#paralelization
cs['X_all_q_ROLEs'] = 'q.k.c.v'
cs['X_all_q_CPU']   = '1.48.1.1'
cs['SE_ROLEs']      = 'q.qp.b'
cs['SE_CPU']        = '1.48.1'
# Flags
cs.arguments.append('ExtendOut')

# B. Merging DBs Steps

db = YamboIn('ypp_rt -q m',filename='ypp.in')
db['Z_input'] = 1.0
db['Actions_and_names'] = [['\"C\"', '\"./COHSEX-T0/ndb.QP\"', '\n\"N\"', '', '\n\"E\"', '\"./GW/ndb.QP\"'],'']

# C. Bethe-Salpeter Step
bs = YamboIn('%s  -r -b -o b -k sex -y d -Q'%yambo_rt,vim=False)  # BS
bs['DBsIOoff'] = 'DIP'
bs['RandQpts'] =  1000000
bs['RandGvec'] =  [ 1, 'RL' ]
bs['CUTGeo']   =  "box z"
bs['CUTBox']   =  [ 0.00, 0.00, 38.0]
# Common variables
bs['FFTGvecs'] = [ 15   , 'Ha'  ]
bs['BSENGexx'] = [ 15   , 'Ha'  ]
bs['BSENGBlk'] = [ 1200 , 'mHa' ]
bs['NGsBlkXs'] = [ 1200 , 'mHa' ]
bs['BndsRnXs'] = [1 ,70]
bs['BSEBands'] = [25,28]

bs['BSEmod']   =  "resonant"
bs['BSSmod']   =  "d"
bs['BEnRange'] =  [ [0.0  , 5.0] , 'eV'  ]    # Energy range spectra
bs['BDmRange'] =  [ [0.05 , 0.05] , 'eV'  ]   # Width
bs['BEnSteps'] =    1000                      # Energy steps
bs.arguments.append('WRbsWF')
bs['X_all_q_ROLEs'] = 'q.k.c.v'
bs['X_all_q_CPU']   = '1.48.1.1'
bs['BS_ROLEs']      = 'k.eh.t'
bs['BS_CPU']        = '48.1.1'

if 'C' in Set:
# Carriers database NEQ COHSEX
  for time in time_probe:
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
    namecs        = 'cohsex-%s-%s%.2feV-damp%dfs-t%dfs'    % ( CSRTmode, source, exc_pump, damp_pump, time )
    cs.write('%s/%s.in' %(dir_inputs,namecs))
    yambo = oarsub(nodes=cs_nodes,core=cs_cores,dependent=0,name='yambo_bse',walltime="10:00:00")
    yambo.add_command('module load yambo/%s'%yambo_module)
    yambo.add_command('mpirun -hostfile \$OAR_NODEFILE %s -F %s/%s.in -J %s -C %s'%(yambo_rt,dir_inputs,namecs,namecs,namecs))
    yambo.run()
    yambo.clean()

# Merging database

if 'M' in Set:
  for time in time_probe:
    namecs      = 'cohsex-%s-%s%.2feV-damp%dfs-t%dfs' % ( CSRTmode, source, exc_pump, damp_pump, time )
    merged_file = 'merged-%s-%s%.2feV-damp%dfs-t%dfs' % ( CSRTmode, source, exc_pump, damp_pump, time )
    db['Actions_and_names'][0][3] = '\"./%s/ndb.QP\"'%namecs
    db.write('%s/%s-merging.in'% (dir_inputs,namecs))
    os.system('%s -F %s/%s-merging.in -J %s'% (ypp_rt,dir_inputs,namecs,merged_file))

# Carriers database BS 

if 'B' in Set:
  for time in time_probe:
    if 'X' in BSRTmode:
      bs['XfnRTdb'] = 'f @ %d fs < ./%s/ndb.RT_carriers' % ( time , link_pump )
    if 'R' in BSRTmode:
      bs['RfnRTdb'] = 'f @ %d fs < ./%s/ndb.RT_carriers' % ( time , link_pump )
    if 'K' in BSRTmode:
      bs['KfnRTdb'] = 'f @ %d fs < ./%s/ndb.RT_carriers' % ( time , link_pump )
    if 'E' in QPdata:
      bs['KfnQPdb'] = 'E < ./GW/ndb.QP'                # GW
      namebs      = 'bse-%s-%s-%s%.2feV-damp%dfs-t%dfs' % ( BSRTmode, QPdata, source, exc_pump, damp_pump, time )
    if 'N' in QPdata:
      namebs      = 'bse-%s-%s-%s%.2feV-damp%dfs-t%dfs' % ( BSRTmode, CSRTmode, source, exc_pump, damp_pump, time )
      name_merged = 'merged-%s-%s%.2feV-damp%dfs-t%dfs' % ( CSRTmode, source, exc_pump, damp_pump, time )
      bs['KfnQPdb'] = 'E < ./%s/ndb.QP' % name_merged  # GW + NEQ_COHSEX - EQ_COHSEX
    bs.write('%s/%s.in' %(dir_inputs,namebs))
    yambo = oarsub(nodes=bs_nodes,core=bs_cores,dependent=0,name='yambo_bse',walltime="10:00:00")
    yambo.add_command('module load yambo/%s'%yambo_module)
    yambo.add_command('mpirun -hostfile \$OAR_NODEFILE %s -F %s/%s.in -J %s -C %s'%(yambo_rt,dir_inputs,namebs,namebs,namebs))
    yambo.run()

# Submission
#create job files
print('running calculation...')
print('Pump pulse features:')
print('Type of calculations %s' % Set)
print('Type of dynamics: %s ' % source)
print('Frequency: %.2f    eV' % exc_pump)
print('Damping:   %d      fs' % damp_pump)
#print 'Time:      %d      fs' % time_probe
print('           ')
