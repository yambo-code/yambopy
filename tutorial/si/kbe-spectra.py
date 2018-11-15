from __future__ import print_function
##############################################################################
#
# Author: Alejandro Molina-Sanchez
# Run real-time simulations with yambo
# 
# Warning: Real-time simulations requires several data folders for running
# properly. Before using this scripts compulsively is recommended
# to understand the different run levels.
#
# This script plots the KBE (delta-pulse) and BSE spectra 
#
##############################################################################
from yambopy     import *
from schedulerpy import *
import sys
import argparse
from numpy import loadtxt

print('-f  Folder containing the rt simulation')
print('-b  Folder containing the bse simulation')
print('-j  Jobname of the rt simulation')
print('-s  Jobname of bse simulation')

parser = argparse.ArgumentParser(description='Map of a double-grid')
parser.add_argument('-f' ,'--folder'   ,help='Folder containing the rt simulation')
parser.add_argument('-b' ,'--bsefolder',help='Folder containing the bse simulation')
parser.add_argument('-j' ,'--jobname'  ,help='Jobname of the rt simulation')
parser.add_argument('-s' ,'--bsename'  ,help='Jobname of bse simulation')

args = parser.parse_args()

folder    = args.folder
jobname   = args.jobname
bsefolder = args.bsefolder
bsename   = args.bsename

run = YamboIn('ypp_rt -t X -V all',folder=folder,filename='ypp.in')
run['EnRngeRt'] = [ [0,10], 'eV']
run['ETStpsRt'] = 1000
run.arguments.append('SkipJP_IO')
run.write('%s/ypp.in' % folder)
os.system('cd %s; ypp_rt -F ypp.in -J %s' % (folder,jobname))

kbe = loadtxt('%s/o-%s.YPP-eps_along_E'  % (folder,jobname))
bse = loadtxt('%s/o-%s.eps_q1_diago_bse' % (bsefolder,bsename))

plt.plot(kbe[:,0],kbe[:,1],label='KBE')
plt.plot(bse[:,0],bse[:,1],label='BSE')

plt.legend()
plt.show()
