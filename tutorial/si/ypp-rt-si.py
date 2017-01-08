##############################################################################
#
# Author: Alejandro Molina-Sanchez
# Run real-time simulations with yambo
# 
# Warning: Real-time simulations requires several data folders for running
# properly. Before using this scripts compulsively is recommended
# to understand the different run levels.
#
# This script calls ypp_rt for post-processing
#
##############################################################################
from yambopy     import *
from schedulerpy import *
import sys
import argparse

folder  = 'rt'
jobname = 'dneqQSSIN-1e+01-75fs-2.23eV-0.0K'

ypp = YamboIn('ypp_rt -n o b -V all',folder=folder,filename='ypp.in')

ypp['TimeRange']   = [ [ 0, 500] , 'fs' ]
ypp['TimeStep']    = [ 50, 'fs' ]
ypp['cooIn']       = "rlu"
ypp['BANDS_steps'] = 10
ypp['QPkrange']    = [1,6,2,7]
ypp['BKpts'] = [[0.5,0.5,0.5],[0.0,0.0,0.0],[0.0,0.5,0.0],[0.25,0.5,-0.25],[0.0,0.0,0.0]]

ypp.write('%s/ypp-obands.in' % folder)
os.system('cd %s; ypp_rt -F ypp-obands.in -J %s' % (folder,jobname) )

ypp = YamboIn('ypp_rt -n o t -V all',folder=folder,filename='ypp.in')
ypp['QPkrange']    = [1,6,2,7]
ypp['TimeRange']   = [ [ 0, 500] , 'fs' ]
ypp['TimeStep']    = [ 50, 'fs' ]
ypp.write('%s/ypp-otime.in' % folder)
os.system('cd %s; ypp_rt -F ypp-otime.in -J %s' % (folder,jobname) )
