#
# Run a GW calculation using Yambo
#
from __future__ import print_function
from yambopy import *
from qepy import *
from schedulerpy import *

# scheduler
scheduler = Scheduler.factory

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

if not os.path.isdir('gw_calc'):
    os.mkdir('gw_calc')
    os.system('cp -r database/SAVE gw_calc')

#create the yambo input file
y = YamboIn.from_runlevel('%s -p p -g n'%yambo,folder='gw_calc')
QPKrange,_ = y['QPkrange']
y['QPkrange'] = [QPKrange[:2]+[2,7],'']
y['FFTGvecs'] = [2000,'RL']
y['NGsBlkXp'] = [10,'RL']
y['BndsRnXp'] = [1,20]
y['GbndRnge'] = [1,20]
y.write('gw_calc/yambo_run.in')

print('running yambo')
shell=scheduler()
shell.add_command('cd gw_calc; %s -F yambo_run.in -J yambo -C yambo'%yambo)
shell.run()
shell.clean()
print('done!')

print('plotting all QP eigenvalues...')

print('packing in json files...')
pack_files_in_folder('gw_calc')
print('done!')

print('reading json files in folder gw_calc ...')
ya = YamboAnalyser('gw_calc')
ya.plot_gw(plot=True)
print('done!')
