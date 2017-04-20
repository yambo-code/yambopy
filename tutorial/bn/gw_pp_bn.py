#
# Author: Henrique Pereira Coutada Miranda
# Date: 18/10/2015
# Parallelize a GW calculation using Yambo
# In this example we take advantage of the fact that the q points of the dielectric function
# are totally independent calculations so its a good idea to calculate them in separate jobs 
# These example runs locally, if you want to make it work behind a queing system you need to
# modify it accordingly
#
from __future__ import print_function
from __future__ import division
from builtins import range
from past.utils import old_div
from yambopy import *
from qepy import *

yambo =  'yambo'
folder = 'gw_pp_par'

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

if not os.path.isdir(folder):
    os.mkdir(folder)
    os.system('cp -r database/SAVE %s'%folder)

#create the yambo input file
y = YamboIn('%s -p p -V all'%yambo,folder=folder)
_,nqpoints = y['QpntsRXp'][0] 


y['FFTGvecs'] = [15,'Ry']
y['NGsBlkXp'] = [1,'Ry']
y['BndsRnXp'] = [[1,30],'']

#prepare the q-points input files
f = open('jobs.sh','w')
for nq in range(1,int(nqpoints)+1):
    y['QpntsRXp'] = [[nq,nq],'']
    y.write('%s/yambo_q%d.in'%(folder,nq))
    if nq != 1:
        f.write('cd %s; %s -F yambo_q%d.in -J %d\n'%(folder,yambo,nq,nq))
f.close()

#calculate first q-point and dipoles
os.system('cd %s; %s -F yambo_q1.in -J 1'%(folder,yambo))
#copy dipoles to save
os.system('cp %s/1/ndb.dip* %s/SAVE'%(folder,folder))
print('running separate yambo files')
os.system('parallel :::: jobs.sh')

#gather all the files
os.system('mkdir %s/yambo'%folder)
os.system('cd %s; cp 1/ndb.pp* yambo'%folder) 
os.system('cd %s; cp */ndb.pp_fragment_* yambo'%folder) 

y = YamboIn('yambo -p p -g n -V all',folder=folder)
QPKrange,_ = y['QPkrange']
y['QPkrange'] = [QPKrange[:2]+[3,6],'']
y['FFTGvecs'] = [15,'Ry']
y['NGsBlkXp'] = [1,'Ry']
y['BndsRnXp'] = [[1,30],'']
y.write('%s/yambo_run.in'%folder)
os.system('cd %s; %s -F yambo_run.in -J yambo'%(folder,yambo))

y.write('%s/yambo_run.in'%folder)

print('running yambo')
os.system('cd %s; %s -F yambo_run.in -J yambo'%(folder,yambo))

#pack the files in .json files
pack_files_in_folder(folder)

#plot the results using yambm analyser
ya = YamboAnalyser()
print(ya)
print('plot all qpoints')
ya.plot_gw('qp')
print('plot along a path')
path = [[   0,   0,   0],
        [ 0.5,   0,   0],
        [old_div(1.,3), old_div(1,3),   0],
        [   0,   0,   0]]
ya.plot_gw_path('qp',path)

print('done!')
