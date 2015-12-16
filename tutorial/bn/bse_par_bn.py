#
# Author: Henrique Pereira Coutada Miranda
# Run a BSE calculation using yambo
# one job per q-point for the dielectric function
#
from __future__ import print_function
from yambopy.inputfile import *
from pwpy.inputfile import *
from pwpy.outputxml import *

yambo = "yambo"

if not os.path.isdir('database'):
    os.mkdir('database')

#check if the nscf cycle is present
if os.path.isdir('nscf/bn.save'):
    print('nscf calculation found!')
else:
    print('nscf calculation not found!')
    exit()

#check if the SAVE folder is present
if not os.path.isdir('database/SAVE'):
    print('preparing yambo database')
    os.system('cd nscf/bn.save; p2y')
    os.system('cd nscf/bn.save; yambo')
    os.system('mv nscf/bn.save/SAVE database')

#check if the SAVE folder is present
if not os.path.isdir('database_double/SAVE'):
    print('preparing yambo database2')
    os.system('cd nscf_double/bn.save; p2y')
    os.system('cd nscf_double/bn.save; yambo')
    os.system('mv nscf_double/bn.save/SAVE database_double')

if not os.path.isdir('bse_par'):
    os.mkdir('bse_par')
    os.system('cp -r database/SAVE bse_par')

#initialize the double grid
if not os.path.isfile('bse_par/SAVE/ndb.Double_Grid'):
    f = open('bse_par/ypp.in','w')
    f.write("""kpts_map
    %DbGd_DB1_paths
    "../database_double"
    %""")
    f.close()
    os.system('cd bse_par; ypp')

#create the yambo input file
y = YamboIn('yambo -b -o b -V all',folder='bse_par')

y['FFTGvecs'] = [15,'Ry']
y['NGsBlkXs'] = [1,'Ry']
y['BndsRnXs'] = [[1,30],'']
y['BSEBands'] = [[3,6],'']
y.arguments.append('WRbsWF')
y.write('bse_par/yambo_run.in')
_,nkpoints = y['QpntsRXs'][0]

#prepare the q-points input files
f = open('jobs.sh','w')
for nk in xrange(1,int(nkpoints)+1):
    y['QpntsRXs'] = [[nk,nk],'']
    y.write('bse_par/yambo_q%d.in'%(nk))
    if nk != 1:
        f.write('cd bse_par; %s -F yambo_q%d.in -J %d\n'%(yambo,nk,nk))
f.close()

#calculate first q-point and dipoles
os.system('cd bse_par; %s -F yambo_q1.in -J 1'%yambo)
#copy dipoles to save
os.system('cp bse_par/1/ndb.dipoles* bse_par/SAVE')

print('running separate yambo files')
os.system('parallel :::: jobs.sh')

#gather all the files
os.system('cp merge_eps.py bse_par')
os.system('cd bse_par; python merge_eps.py')

os.system('rm bse_par/yambo.in')
y = YamboIn('yambo -b -o b -k sex -y d -V all',folder='bse_par')
y['FFTGvecs'] = [15,'Ry']
y['NGsBlkXs'] = [1,'Ry']
y['BndsRnXs'] = [[1,30],'']
y['BSEBands'] = [[3,6],'']
y.arguments.append('WRbsWF')
y.write('bse_par/yambo_run.in')
os.system('cd bse_par; %s -F yambo_run.in -J yambo'%yambo)


