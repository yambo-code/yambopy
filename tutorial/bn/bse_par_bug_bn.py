#
# Author: Henrique Pereira Coutada Miranda
# Run a BSE calculation using yambo
# one job per q-point for the dielectric function
# This script has a problem with the current version of yambo.
# the indexation of the dielectric function database fragments is made using the databases present
# and not the q point (is this true?)
# We calculate the q point 5 but its stored in fragment 1, in the last step of the calculation
# all th q-points are calculated except 1st qpoint but the data of qpoint 5 is there isntead
#
from __future__ import print_function
from yambopy import *
from qepy import *
import os

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

if not os.path.isdir('bse_par_bug'):
    os.mkdir('bse_par_bug')
    os.system('cp -r database/SAVE bse_par_bug')

#initialize the double grid
if not os.path.isfile('bse_par_bug/SAVE/ndb.Double_Grid'):
    f = open('bse_par_bug/ypp.in','w')
    f.write("""kpts_map
    %DbGd_DB1_paths
    "../database_double"
    %""")
    f.close()
    os.system('cd bse_par_bug; ypp')

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

print('running yambo q-point 5 calculation')
os.system('cd bse_par_bug; %s -F yambo_q5.in -J yambo'%yambo)

print('running final yambo')
y = YamboIn('yambo -b -o b -k sex -y d -V all',folder='bse_par_bug')
y['FFTGvecs'] = [15,'Ry']
y['NGsBlkXs'] = [1,'Ry']
y['BndsRnXs'] = [[1,30],'']
y['BSEBands'] = [[3,6],'']
y.arguments.append('WRbsWF')
y.write('bse_par_bug/yambo_run.in')
os.system('cd bse_par_bug; %s -F yambo_run.in -J yambo'%yambo)
print('done!')
