#
# Author: Henrique Pereira Coutada Miranda
# Run a GW+BSE calculation using Yambo
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

if not os.path.isdir('gw_bse'):
    os.mkdir('gw_bse')
    shell.run('cp -r database/SAVE gw_bse'.split())

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
log_yambo = shell.run(('%s -F yambo_run.in -J yambo'%yambo).split(),cwd='gw_bse')
print(log_yambo.output,file=open("yambo_gw.log","w"))

#creathe the bse input file
shell.run('rm gw_bse/yambo.in'.split())
y = YamboIn('%s -b -o b -k sex -y d -V all'%yambo,folder='gw_bse')
y['QPkrange'] = [QPKrange[:2]+[6,10],'']
y['FFTGvecs'] = [15,'Ry']
y['NGsBlkXs'] = [1,'Ry']
y['BndsRnXs'] = [[1,30],'']
y.write('gw_bse/yambo_run.in')

#run the bse calculation using the dielectric function from gw
print('running bse')
log_yambo = shell.run(('%s -F yambo_run.in -J yambo'%yambo).split(),cwd='gw_bse')
print(log_yambo.output,file=open("yambo_bse.log","w"))
