#
# Author: Henrique Pereira Coutada Miranda
# Run a BSE calculation using yambo
# one job per q-point for the dielectric function
# This script has a problem with the current version of yambo.
# the indexation of the dielectric function database fragments is made using the databases present
# and not the q point (is this true?)
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
    print('preparing yambo database2')
    log_p2y   = shell.run(['p2y'],  cwd="nscf_double/bn.save")
    print(log_p2y.output,file=open("p2y.log","w"))
    log_yambo = shell.run(['yambo'],cwd="nscf_double/bn.save")
    print(log_yambo.output,file=open("yambo.log","w"))
    shell.run('mv nscf_double/bn.save/SAVE database_double'.split())

if not os.path.isdir('bse_par_bug'):
    os.mkdir('bse_par_bug')
    shell.run('cp -r database/SAVE bse_par_bug'.split())

#initialize the double grid
if not os.path.isfile('bse_par_bug/SAVE/ndb.Double_Grid'):
    f = open('bse_par_bug/ypp.in','w')
    f.write("""kpts_map
    %DbGd_DB1_paths
    "../database_double"
    %""")
    f.close()
    shell.run('ypp',cwd='bse_par_bug')

#create the yambo input file
y = YamboIn('yambo -b -V all',folder='bse_par_bug')

y['FFTGvecs'] = [15,'Ry']
y['NGsBlkXs'] = [1,'Ry']
y['BndsRnXs'] = [[1,30],'']
y['BSEBands'] = [[3,6],'']
y.arguments.append('WRbsWF')
y.write('bse_par_bug/yambo_run.in')
_,nkpoints = y['QpntsRXs'][0]

#prepare the q-points input files
f = open('jobs.sh','w')
for nk in xrange(1,int(nkpoints)+1):
    y['QpntsRXs'] = [[nk,nk],'']
    y.write('bse_par_bug/yambo_q%d.in'%(nk))
    f.write('cd bse_par_bug; %s -F yambo_q%d.in -J yambo\n'%(yambo,nk))
f.close()

print('running separate yambo files')
log_yambo = shell.run('parallel :::: jobs.sh'.split())
print(log_yambo.output,file=open("yambo_par_bse.log","w"))

print('running final yambo')
y = YamboIn('yambo -b -o b -k sex -y d -V all',folder='bse_par_bug')
y['FFTGvecs'] = [15,'Ry']
y['NGsBlkXs'] = [1,'Ry']
y['BndsRnXs'] = [[1,30],'']
y['BSEBands'] = [[3,6],'']
y.arguments.append('WRbsWF')
y.write('bse_par_bug/yambo_run.in')
shell.run(('%s -F yambo_run.in -J yambo'%yambo).split(),cwd='bse_par_bug')
print('done!')
