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
from builtins import range
from yambopy     import *
from schedulerpy import *

############## SETTINGS ##############

ypp_rt       = 'ypp_rt'

source       = 'QSSIN-1e+03-70.0fs-2.0eV-0K'
folder_rt    = 'rt-6x6'
# folder with ndb.QP from GW calculation (check line 46)
# sym must have been broken
folder_gw    = 'gw-6x6'

CSRTmode = 'XG'  #  X: Screening, G: GFs


time_probe = list(range(0,610,150))
#time_probe=(0,)
print(time_probe)


# Merging DBs Steps
db = YamboIn('%s -q m' % ypp_rt,filename='ypp.in',folder=folder_rt)
db['Z_input'] = 1.0
db['Actions_and_names'] = [['\"C\"', '', '\n\"N\"', '', '\n\"E\"', '\"../../%s/ndb.QP\"' % folder_gw],'']
#db['Actions_and_names'] = [['\"C\"', '', '\n\"N\"', ''],''] # Not including GW


# Submission of the jobs for each time

print('Starting merges.')
print('RT source calculation: %s \n' % source)

for time in time_probe:
    print('Time of carriers database %d' % time)
    namecs0     = 'C-%s-t0'  % ( CSRTmode )
    namecs      = 'C-%s-t%d' % ( CSRTmode, time )
    merged_file = 'M-%s-t%d' % ( CSRTmode, time )
    db['Actions_and_names'][0][1] = '\"./%s/ndb.QP\"'%namecs0 # Fetches ndb.QP from COHSEX at time 0
    db['Actions_and_names'][0][3] = '\"./%s/ndb.QP\"'%namecs # Fetches ndb.QP from COHSEX at time tau
    db.write('%s/%s/%s.in'% (folder_rt, source, merged_file))
    os.system('cd %s/%s; %s -F %s.in -J %s -I ..'% (folder_rt,source,ypp_rt,merged_file,merged_file))
