#
# Author: Henrique Pereira Coutada Miranda
# Run a GW calculation using Yambo
#
from __future__ import print_function
from yambopy.inputfile import *
from pwpy.inputfile import *
from pwpy.outputxml import *
import spur

shell = spur.LocalShell()
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
    log_p2y   = shell.run(['p2y'],  cwd="nscf/si.save")
    print(log_p2y.output,file=open("p2y.log","w"))
    log_yambo = shell.run(['yambo'],cwd="nscf/si.save")
    print(log_yambo.output,file=open("yambo.log","w"))
    shell.run('mv nscf/si.save/SAVE database'.split())

if not os.path.isdir('gw'):
    os.mkdir('gw')
    shell.run('cp -r database/SAVE gw'.split())

#create the yambo input file
y = YamboIn('%s -d -p c -g n -V all'%yambo,folder='gw')
QPKrange,_ = y['QPkrange']
y['QPkrange'] = [QPKrange[:2]+[6,10],'']
y['FFTGvecs'] = [15,'Ry']
y['NGsBlkXp'] = [1,'Ry']
y['BndsRnXp'] = [[1,30],'']
y.arguments.append('WFbuffIO')
y.write('gw/yambo_run.in')

print('running yambo')
log_yambo = shell.run(('%s -F yambo_run.in -J yambo'%yambo).split(),cwd='gw')
print(log_yambo.output,file=open("yambo_gw.log","w"))
