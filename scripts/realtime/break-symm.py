##############################################################################
#
# Author: Alejandro Molina-Sanchez
# Run real-time simulations with yambo
# 
# Warning: Real-time simulations requires several data folders for running
# properly. Before using this scripts compulsively is recommended
# to understand the different run levels.
#
# This script breaks the symmetries of a given nscf run and
# prepares the rt-folder for the RT simulations
#
##############################################################################
from __future__ import print_function
from yambopy     import *
import argparse

print('This script breaks the symmetries of a given nscf run and')
print('prepares the rt-folder for the RT simulations')
print('Any previous "database" folder will be erased')
print('arg1: folder with nscf data')
print('arg2: folder for RT simulation')
print('arg3: prefix')
print('arg4: symmetry')
print('Example: python break-symm.py -i nscf -o rt-100 -p si -s 100')

parser = argparse.ArgumentParser(description='Map of a double-grid')
parser.add_argument('-i' ,'--input'    , help='Folder with nscf data')
parser.add_argument('-o' ,'--output'   , help='Folder for RT simulation')
parser.add_argument('-p' ,'--prefix'   ,help='Prefix of nscf calculation')
parser.add_argument('-s' ,'--symmetry' ,help='Choice symmetry: 100 010 110')
args = parser.parse_args()

print('Folder of nscf data     ===>>>  ' ,args.input)
print('Folder of RT simulation ===>>>  ' ,args.output)
print('Prefix                  ===>>>  ' ,args.prefix)
print('Symmetry                ===>>>  ' ,args.symmetry)

nscf_folder = args.input
rt_folder   = args.output
prefix      = args.prefix
symm        = args.symmetry
  
# Check rt_folder
if os.path.isdir(rt_folder):
    print('Real-time folder already exists')
    raise SystemExit
# Generation of the database folder
if os.path.isdir('database'):
    os.system('rm -rf database')
os.system('cd %s/%s.save ; p2y -O ../../database' % (nscf_folder, prefix))
os.system('cd database; yambo')
  
# Breaking of symmetries

if symm   == '100':
    breaking_symmetries([1,0,0], [0,0,0], rt_folder)
elif symm == '010':
    breaking_symmetries([0,1,0], [0,0,0], rt_folder)
elif symm == '110':
    breaking_symmetries([1,0,0], [0,1,0], rt_folder)
