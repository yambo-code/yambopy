#
# Author: Henrique Pereira Coutada Miranda
# Run a GW+BSE calculation using Yambo
#
from __future__ import print_function
from yambopy import *
from qepy import *

yambo =  'yambo'

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

if not os.path.isdir('gw_bse'):
    os.mkdir('gw_bse')
    os.system('cp -r database/SAVE gw_bse')

#create the gw yambo input file
y = YamboIn('%s -d -p c -g n -V all'%yambo,folder='gw_bse')
QPKrange,_ = y['QPkrange']
y['QPkrange'] = [QPKrange[:2]+[6,10],'']
y['FFTGvecs'] = [15,'Ry']
y['NGsBlkXs'] = [1,'Ry']
y['BndsRnXs'] = [[1,30],'']
y.arguments.append('WFbuffIO')
y.write('gw_bse/yambo_run.in')

print('running gw')
os.system('cd gw_bse; %s -F yambo_run.in -J yambo'%yambo)

#creathe the bse input file
y = YamboIn('%s -b -o b -k sex -y d -V all'%yambo,folder='gw_bse')
y['FFTGvecs'] = [15,'Ry']
y['NGsBlkXs'] = [1,'Ry']
y['BndsRnXs'] = [[1,30],'']
y.write('gw_bse/yambo_run.in')

#run the bse calculation using the dielectric function from gw
print('running bse')
os.system('cd gw_bse; %s -F yambo_run.in -J yambo'%yambo)
