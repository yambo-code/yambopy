#
# Author: Henrique Pereira Coutada Miranda
# Run a GW calculation using yambo
#
from __future__ import print_function
from yambopy.inputfile import *
from pwpy.inputfile import *
from pwpy.outputxml import *
import spur

yambo = "yambo"

if not os.path.isdir('database'):
    os.mkdir('database')

shell = spur.LocalShell()

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

if not os.path.isdir('bse'):
    os.mkdir('bse')
    shell.run('cp -r database/SAVE bse'.split())

#create the yambo input file
y = YamboIn('yambo -b -o b -k sex -y d -V all',folder='bse')

y['FFTGvecs'] = [15,'Ry']
y['NGsBlkXs'] = [1,'Ry']
y['BndsRnXs'] = [[1,30],'']
y.arguments.append('WFbuffIO')
y.write('bse/yambo_run.in')

print('running yambo')
log_yambo = shell.run(('%s -F yambo_run.in -J yambo'%yambo).split(),cwd='bse')
print(log_yambo.output,file=open("yambo_bse.log","w"))


