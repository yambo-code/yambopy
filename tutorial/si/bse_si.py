#
# Run a BSE calculation using yambo
#
from __future__ import print_function
from yambopy import *
from qepy import *
from schedulerpy import *

# scheduler
scheduler = Scheduler.factory

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

if not os.path.isdir('bse_calc'):
    os.mkdir('bse_calc')
    os.system('cp -r database/SAVE bse_calc')

#create the yambo input file
y = YamboIn.from_runlevel('yambo -r -b -o b -k sex -y d -V all',folder='bse_calc')

y['FFTGvecs'] = [5,'Ha']
y['BSENGexx'] = [5,'Ha']
y['NGsBlkXs'] = [800,'mHa']
y['BSENGBlk'] = [800,'mHa']
y['BndsRnXs'] = [1,20]
y['BSEBands'] = [2, 7]
y['RandQpts'] = 1000000
y['BEnSteps'] = 1000
y.arguments.append('WRbsWF')
y.arguments.append('ALLGexx')
y.write('bse_calc/yambo_run.in')

print('running yambo')
shell=scheduler()
shell.add_command('cd  bse_calc; %s -F yambo_run.in -J yambo -C yambo'%yambo)
shell.run()
shell.clean()
print('done!')
