#########################################################
#  YAM(BO)PY(THON) Library
#
#  Generation of Yambo input files using python
#
#  Authors: A Molina-Sanchez, HPC Miranda
#
#  January 2016
#########################################################
# Merges the QP corrections from equilibrium (GW), carriers
# at t=0 and at t=tau.
# /!\ Doesn't schedule jobs, calls yambo/ypp directly.
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

ypp_rt       = 'ypp_rt'

source       = 'rt-D-1.94eV-0K-0.1fs-DG'
#source       = 'QSSIN-1e+03-70.0fs-1.94eV-0K' # no dissipation
folder_kerr  = 'kerr-24x24'
folder_rt    = 'rt-24x24'
folder_gw    = 'gw-24x24' # folder with ndb.QP from GW calculation (check line 60)
folder_cohsex= 'cohsex-24x24' # contains ndb.QP @ t = 0

CSRTmode = 'XG'  #  X: Screening, G: GFs


time_probe = range(0,610,150)
#time_probe=(0,)
print(time_probe)

dir_pump   = '../%s/%s/' % (folder_rt,source)
link_pump  = '../%s/%s/' % (folder_rt,source)
dir_inputs = 'inputs'
os.system('cd %s; mkdir -p %s'%(folder_kerr,dir_inputs))


# Merging DBs Steps
db = YamboIn('%s -q m' % ypp_rt,filename='ypp.in',folder=folder_kerr)
db['Z_input'] = 1.0
db['Actions_and_names'] = [['\"C\"', '\"../%s/ndb.QP\"' % folder_cohsex, '\n\"N\"', '', '\n\"E\"', '\"../%s/FixSymm/24x24_run/ndb.QP\"' % folder_gw],'']
#db['Actions_and_names'] = [['\"C\"', '\"../%s/ndb.QP\"' % folder_cohsex, '\n\"N\"', ''],'']


# Submission of the jobs for each time

print('Starting merges.')
print('RT source calculation: %s \n' % source)

for time in time_probe:
    print('Time of carriers database %d' % time)
    namecs      = 'C-%s-%s-t%d' % ( CSRTmode, source, time )
    merged_file = 'M-%s-%s-t%d' % ( CSRTmode, source, time )
    db['Actions_and_names'][0][3] = '\"./%s/ndb.QP\"'%namecs # Fetches ndb.QP from COHSEX at time tau
    db.write('%s/%s/%s-merging.in'% (folder_kerr, dir_inputs, namecs))
    os.system('cd %s; %s -F %s/%s-merging.in -J %s'% (folder_kerr,ypp_rt,dir_inputs,namecs,merged_file))
