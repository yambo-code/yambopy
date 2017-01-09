##############################################################################
#
# Author: Alejandro Molina-Sanchez
# Run real-time simulations with yambo
# 
# Warning: Real-time simulations requires several data folders for running
# properly. Before using this scripts compulsively is recommended
# to understand the different run levels.
#
# This script map a fine grid to a coarse grid
#
##############################################################################
#from __future__ import print_function
from sys import argv
from yambopy     import *
import argparse

print ('This script map a fine grid to a coarse grid.')
print ('It requires three arguments')
print ('1: -i  folder with the fine grid')
print ('2: -o  folder of the RT simulation')
print ('3: -dg name for the folder hosting the double-grid')

parser = argparse.ArgumentParser(description='Map of a double-grid')
parser.add_argument('-i'  ,'--input'    , help='Folder containing the SAVE folder of the double grid')
parser.add_argument('-o'  ,'--output'   , help='Output folder (contains the rt simulation SAVE)')
parser.add_argument('-dg' ,'--folder_dg',help='Folder containing the mapped double grid')
args = parser.parse_args()

print(args.input)
print(args.output)
print(args.folder_dg)
folder_in  = args.input 
folder_out = args.output
folder_dg  = args.folder_dg

sym = YamboIn('ypp_rt -m',folder=folder_out,filename='ypp.in')
sym['DbGd_DB1_paths']= [ ["'../%s'" % folder_in], '' ]
sym.arguments.append('noBZExpand')
sym.write('%s/map-dg.in' % (folder_out))
os.system('cd %s; ypp_rt -F map-dg.in' % (folder_out))
os.system('cd %s; mkdir -p %s; mv SAVE/ndb.Double_Grid %s/' % (folder_out, folder_dg, folder_dg)) 

