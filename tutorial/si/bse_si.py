#
# Author: Henrique Pereira Coutada Miranda
# Run a GW calculation using yambo
#
from __future__ import print_function
from yambopy.inputfile import *
from pwpy.inputfile import *
from pwpy.outputxml import *

yambo = "yambo"

if not os.path.isdir('database'):
    os.mkdir('database')

#check if the nscf cycle is present
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

if not os.path.isdir('bse'):
    os.mkdir('bse')
    os.system('cp -r database/SAVE bse')

#create the yambo input file
y = YamboIn('yambo -b -o b -k sex -y d -V all',folder='bse')

y['FFTGvecs'] = [15,'Ry']
y['NGsBlkXs'] = [1,'Ry']
y['BndsRnXs'] = [[1,30],'']
y.arguments.append('WFbuffIO')
y.write('bse/yambo_run.in')

print('running yambo')
os.system('cd  bse; %s -F yambo_run.in -J yambo'%yambo)


