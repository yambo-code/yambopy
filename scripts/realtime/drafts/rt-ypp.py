from __future__ import print_function
from yambopy     import *
from schedulerpy import *
import sys
import glob
import argparse

"""
Automation of ypp_rt calls.
   - Spectra : BSE spectra of real-time simulations (NOT on time t)
   - Occupation : Calculates occupation at different times
"""

parser = argparse.ArgumentParser(description='')
parser.add_argument('-f' ,'--folder'    , help='Folder with real-time simulations')
parser.add_argument('-j' ,'--job'       , help='Name of job (ex: DELTA-1E+03)')
parser.add_argument('-s' ,'--spectra'   , action="store_true", help='Do a BSE calculation')
parser.add_argument('-o' ,'--occupation', action="store_true", help='Compute occupation at different times')
args = parser.parse_args()

folder = args.folder
job = args.job

### Occupation path
path = [[0.0,0.0,0.0],[0.5,0.5,0.0]]


###
# This snippet of code can be used to automatize the "spectra" calculation over several jobs, using a prefix (here QSSIN in the if startswith)

##calculation = sys.argv[2]
#calculations = []
##calculations = glob.glob('QSSIN*')
#ls = os.listdir("./%s/"%folder)
#for i,calc in enumerate(ls):
#  if calc.startswith('QSSIN')==True:
#    calculations.append(calc)
#
#print calculations
###

if args.spectra:
  #for calculation in calculations:
  #  print calculation
    run = YamboIn('ypp_rt -t X -V all',folder=folder,filename='ypp.in')
    run['EnRngeRt'] = [ [0,5], 'eV']
    run['ETStpsRt'] = 1000
    run.arguments.append('SkipJP_IO')
    run.write('%s/ypp-spectra.in' % folder)
    os.system('cd %s;ypp_rt -F ypp-spectra.in -J %s' % (folder,job))

if args.occupation:
    run = YamboIn('ypp_rt -n o b -V all',folder=folder,filename='ypp.in')

    run['TimeRange']   = [ [ 0, 50] , 'fs' ]
    run['TimeStep']    = [ 50, 'fs' ]
    run['cooIn']       = "rlu"
    run['BANDS_steps'] = 10
    run['QPkrange']    = [1,576,25,28]
    run.arguments.append('NNInterp')
    run['BKpts'] = [path,'']
    run.write('%s/ypp-obands.in' % folder)
    os.system('cd %s; ypp_rt -F ypp-obands.in -J %s' % (folder,job) )

