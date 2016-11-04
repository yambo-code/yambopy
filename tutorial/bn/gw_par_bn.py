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
from yambopy import *
from qepy import *

yambo =  'yambo'

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

if not os.path.isdir('gw_par'):
    os.mkdir('gw_par')
    os.system('cp -r database/SAVE gw_par')

#create the yambo input file
y = YamboIn('%s -d -V all'%yambo,folder='gw_par')
_,nqpoints = y['QpntsRXd'][0] 


y['FFTGvecs'] = [15,'Ry']
y['NGsBlkXd'] = [1,'Ry']
y['BndsRnXd'] = [[1,30],'']
y['ETStpsXd'] = [50,'']
#Here we write a range that does not make sense (upper bownd smaller than lower bound)
#So that yambo has to calculate it as it would be done when running a calculation in the gw0 runlevel
#The following code has to be cahnged in src/pol_function/X_em1.F
#     if (self_detect_E_range) then
#       call X_eh_setup(-iq,X,Xen,Xk,minmax_ehe)
#       deallocate(X_poles)
#       Xw%er=minmax_ehe
#     endif
#+     ! Start modified by Henrique Miranda
#+     ! If the EnRngeXd variable range does not make sense (upper bownd smaller than lower bound)
#+     ! we detect the range as in the case self_detect_E_range = .True.
#+     if (Xw%er(2) .lt. Xw%er(1)) then
#+       call X_eh_setup(-iq,X,Xen,Xk,minmax_ehe)
#+       deallocate(X_poles)
#+       Xw%er=minmax_ehe
#+     endif
#+     ! end modified by Henrique Miranda

y['EnRngeXd'] = [[1,-1],'eV']

#prepare the q-points input files
f = open('jobs.sh','w')
for nq in xrange(1,int(nqpoints)+1):
    y['QpntsRXd'] = [[nq,nq],'']
    y.write('gw_par/yambo_q%d.in'%nq)
    if nq != 1:
        f.write('cd gw_par; %s -F yambo_q%d.in -J %d\n'%(yambo,nq,nq))
f.close()

#calculate first q-point and dipoles
os.system('cd gw_par; %s -F yambo_q1.in -J 1'%yambo)
#copy dipoles to save
os.system('cp gw_par/1/ndb.dip* gw_par/SAVE')
print('running separate yambo files')
os.system('parallel :::: jobs.sh')

#gather all the files
#os.system('cp merge_eps.py gw_par')
os.system('mkdir gw_par/yambo')
os.system('cd gw_par; cp 1/ndb.em1* yambo') 
os.system('cd gw_par; cp */ndb.em1?_fragment_* yambo') 

y = YamboIn('yambo -d -g n -V all',folder='gw_par')
QPKrange,_ = y['QPkrange']
y['QPkrange'] = [QPKrange[:2]+[3,6],'']
y['FFTGvecs'] = [15,'Ry']
y['NGsBlkXd'] = [1,'Ry']
y['BndsRnXd'] = [[1,30],'']
y['ETStpsXd'] = [50,'']
y.write('gw_par/yambo_run.in')
os.system('cd gw_par; %s -F yambo_run.in -J yambo'%yambo)

y.write('gw_par/yambo_run.in')

print('running yambo')
os.system('cd gw_par; %s -F yambo_run.in -J yambo'%yambo)

#pack the files in .json files
pack_files_in_folder('gw_par')

#plot the results using yambm analyser
ya = YamboAnalyser()
print(ya)
print('plot all qpoints')
ya.plot_gw('qp')
print('plot along a path')
path = [[   0,   0,   0],
        [ 0.5,   0,   0],
        [1./3, 1/3,   0],
        [   0,   0,   0]]
ya.plot_gw_path('qp',path)

print('done!')
