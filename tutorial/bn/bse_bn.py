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
if os.path.isdir('nscf/bn.save'):
    print('nscf calculation found!')
else:
    print('nscf calculation not found!')
    exit()

#check if the SAVE folder is present
if not os.path.isdir('database/SAVE'):
    print('preparing yambo database')
    log_p2y   = shell.run(['p2y'],  cwd="nscf/bn.save")
    print(log_p2y.output,file=open("p2y.log","w"))
    log_yambo = shell.run(['yambo'],cwd="nscf/bn.save")
    print(log_yambo.output,file=open("yambo.log","w"))
    shell.run('mv nscf/bn.save/SAVE database'.split())

#check if the SAVE folder is present
if not os.path.isdir('database_double/SAVE'):
    print('preparing yambo database')
    log_p2y   = shell.run(['p2y'],  cwd="nscf_double/bn.save")
    print(log_p2y.output,file=open("p2y.log","w"))
    log_yambo = shell.run(['yambo'],cwd="nscf_double/bn.save")
    print(log_yambo.output,file=open("yambo.log","w"))
    shell.run('mv nscf_double/bn.save/SAVE database_double'.split())

if not os.path.isdir('bse'):
    os.mkdir('bse')
    shell.run('cp -r database/SAVE bse'.split())

#initialize the double grid
f = open('bse/ypp.in','w')
f.write("""kpts_map
%DbGd_DB1_paths
"../database_double"
%""")
f.close()
shell.run('ypp',cwd='bse')

#create the yambo input file
y = YamboIn('yambo -b -o b -k sex -y d -V all',folder='bse')

y['FFTGvecs'] = [15,'Ry']
y['NGsBlkXs'] = [1,'Ry']
y['BndsRnXs'] = [[1,30],'']
y['BSEBands'] = [[3,6],'']
y.arguments.append('WFbuffIO')
y.arguments.append('WRbsWF')
y.write('bse/yambo_run.in')

print('running yambo')
log_yambo = shell.run(('%s -F yambo_run.in -J yambo'%yambo).split(),cwd='bse')
print(log_yambo.output,file=open("yambo_bse.log","w"))


