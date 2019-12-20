#
# Author: Henrique Pereira Coutada Miranda
# Run a GW calculation using yambo
#
from __future__ import print_function
from yambopy import *
from qepy import *
import subprocess

def run(filename):
    """ Function to be called by the optimize function """
    folder = filename.split('.')[0]
    print(filename, folder)
    os.system('cd bse_conv; yambo -F %s -J %s -C %s 2> %s.log'%(filename,folder,folder,folder))

if not os.path.isdir('database'):
    os.mkdir('database')

#check if the nscf data is present
if os.path.isdir('nscf/si.save'):
    print('nscf calculation found!')
else:
    print('nscf calculation not found!')
    exit()

#check if the SAVE folder is present
if not os.path.isdir('database/SAVE'):
    print('preparing yambo database')
    os.system('cd nscf/si.save; p2y')
    os.system('cd nscf/si.save; yambo')
    os.system('mv nscf/si.save/SAVE database')

#if bse folder is not present, create it
if not os.path.isdir('bse_conv'):
    os.mkdir('bse_conv')
    os.system('cp -r database/SAVE bse_conv')

#create the yambo input file
y = YamboIn.from_runlevel('yambo -b -o b -k sex -y h -V all',folder='bse_conv')

#list of variables to optimize and the values they might take
conv = { 'FFTGvecs': [[2,5,10,15,20],'Ry'],
         'NGsBlkXs': [[0,1,2,5], 'Ry'],
         'BndsRnXs': [[1,10],[1,20],[1,30]] }

y.optimize(conv,folder='bse_conv',run=run)
