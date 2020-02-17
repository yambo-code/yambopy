#
# Author: Henrique Pereira Coutada Miranda
# Run a GW calculation using Yambo
#
from __future__ import print_function
from yambopy import *
from qepy import *
from schedulerpy import *

yambo =  'yambo'
scheduler = Scheduler.factory

#check if the nscf cycle is present
if os.path.isdir('nscf/bn.save'):
    print('nscf calculation found!')
else:
    print('nscf calculation not found!')
    exit()

#check if the SAVE folder is present
if not os.path.isdir('database/SAVE'):
    print('preparing yambo database')
    p2y_run = scheduler()
    p2y_run.add_command('mkdir -p database')    
    p2y_run.add_command('cd nscf/bn.save; p2y > p2y.log')
    p2y_run.add_command('yambo > yambo.log')
    p2y_run.add_command('mv SAVE ../../database')
    p2y_run.run()

if not os.path.islink('gw/SAVE'):
    s = scheduler()
    s.add_command('mkdir -p gw')
    s.add_command('cd gw; ln -s ../database/SAVE')
    s.run()

#create the yambo input file
y = YamboIn.from_runlevel('%s -p p -g n -V all'%yambo,folder='gw')
    
y['EXXRLvcs'] = [60,'Ry']       # Self-energy. Exchange
y['BndsRnXp'] = [1,40]          # Screening. Number of bands
y['NGsBlkXp'] = [3,'Ry']        # Cutoff Screening
y['GbndRnge'] = [1,30]          # Self-energy. Number of bands

#read values from QPkrange
values, units = y['QPkrange']
kpoint_start, kpoint_end, band_start, band_end = values
#set the values of QPkrange
y['QPkrange'] = [kpoint_start,kpoint_end,2,6]

y.write('gw/yambo_run.in')

print('running yambo')
yambo_run = scheduler()
yambo_run.add_command('cd gw; %s -F yambo_run.in -J yambo'%yambo)
yambo_run.run()
