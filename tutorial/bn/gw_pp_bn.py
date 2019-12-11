#
# Author: Henrique Pereira Coutada Miranda
# Date: 18/10/2015
# Parallelize a GW calculation using Yambo
# In this example we take advantage of the fact that the q points of the dielectric function
# are totally independent calculations so its a good idea to calculate them in separate jobs 
# These example runs locally, if you want to make it work behind a queing system you need to
# modify it accordingly
#
#from __future__ import print_function, import division
from yambopy import *
from qepy import *
from schedulerpy import *
import matplotlib.pyplot as plt
import argparse

#parse options
parser = argparse.ArgumentParser(description='Test the yambopy script.')
parser.add_argument('-c', '--calc', action="store_true", help='calculate the manually parallelised GW')
parser.add_argument('-p', '--plot', action="store_true", help='plot the results')
args = parser.parse_args()

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

yambo =  'yambo'
folder = 'gw_pp_par'
scheduler = Scheduler.factory

#check if the nscf cycle is present
if os.path.isdir('nscf/bn.save'):
    print('nscf calculation found!')
else:
    print('nscf calculation not found!')
    exit() 

#check if the SAVE folder is present
if not os.path.isdir('database/SAVE'):
    print('preparing yambo database')
    p2y_run = scheduler()
    p2y_run.add_command('mkdir -p database')
    p2y_run.add_command('cd nscf/bn.save; p2y > p2y.log')
    p2y_run.add_command('yambo > yambo.log')
    p2y_run.add_command('mv SAVE ../../database/')
    p2y_run.run()

if not os.path.isdir('%s/SAVE'%folder):
    s = scheduler()
    s.add_command('mkdir -p %s'%folder)
    s.add_command('cd %s; cp -r ../database/SAVE .'%folder)
    s.run()

if args.calc:
    #create the yambo input file for dipoles and plasmon pole
    y = YamboIn.from_runlevel('%s -p p -V all'%yambo,folder=folder)
    _,nqpoints = y['QpntsRXp'][0] 

    y['FFTGvecs'] = [15,'Ry']
    y['NGsBlkXp'] = [1,'Ry']
    y['BndsRnXp'] = [[1,30],'']

    #prepare the q-points input files
    for nq in range(1,int(nqpoints)+1):
        y['QpntsRXp'] = [[nq,nq],'']
        y.write('%s/yambo_q%d.in'%(folder,nq))
    #calculate first q-point and dipoles
    shell = scheduler()
    shell.add_command('cd %s; %s -F yambo_q1.in -J 1'%(folder,yambo))
    shell.add_command('mkdir -p yambo')
    #gather the dipoles in save
    shell.add_command('cp 1/ndb.dipoles* SAVE') 
    shell.add_command('cp 1/ndb.pp SAVE')
    #run all
    print('running separate yambo files')
    for nq in range(2,int(nqpoints)+1): shell.add_command('%s -F yambo_q%s.in -J %s'%(yambo,str(nq),str(nq)))
    shell.run()

    #create the yambo input file for gw
    y = YamboIn.from_runlevel('yambo -p p -g n -V all',folder=folder)
    QPKrange,_ = y['QPkrange']
    y['QPkrange'] = [QPKrange[:2]+[3,6],'']
    y['FFTGvecs'] = [15,'Ry']
    y['NGsBlkXp'] = [1,'Ry']
    y['BndsRnXp'] = [[1,30],'']
    y.write('%s/yambo_run.in'%folder)

    yambo_run = scheduler()
    yambo_run.add_command('cd %s'%folder)
    #gather all files
    for nq in range(1,int(nqpoints)+1): 
        yambo_run.add_command('cp {0}/ndb.pp_fragment* yambo'.format(str(nq))) 
        yambo_run.add_command('mv l-{0}* r-{1}* {2}'.format(str(nq),str(nq),str(nq)))
    print('running gw')
    yambo_run.add_command('%s -F yambo_run.in -J yambo'%yambo)
    yambo_run.run()

if args.plot:

    # Define path in reduced coordinates using Class Path
    npoints = 10
    path = Path([ [[  0.0,  0.0,  0.0],'$\Gamma$'],
                  [[  0.5,  0.0,  0.0],'M'],
                  [[1./3.,1./3.,  0.0],'K'],
                  [[  0.0,  0.0,  0.0],'$\Gamma$']], [int(npoints*2),int(npoints),int(sqrt(5)*npoints)] )

    # Read Lattice information from SAVE
    lat  = YamboSaveDB.from_db_file(folder='%s/SAVE'%folder,filename='ns.db1')
    # Read QP database
    y    = YamboQPDB.from_db(filename='ndb.QP',folder='%s/yambo'%folder)

    # 2. Plot of KS and QP eigenvalues NOT interpolated along the path
    ks_bs_0, qp_bs_0 = y.get_bs_path(lat,path)

    fig = plt.figure(figsize=(4,5))
    ax = fig.add_axes( [ 0.20, 0.20, 0.70, 0.70 ])

    ks_bs_0.plot_ax(ax,legend=True,color_bands='r',label='KS')
    qp_bs_0.plot_ax(ax,legend=True,color_bands='b',label='QP-GW')

    plt.show()
