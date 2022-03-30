#
# Author: Henrique Pereira Coutada Miranda
# Run a IP calculation using yambo
#
from __future__ import print_function
import sys
from yambopy import *
from qepy import *
import argparse
from schedulerpy import *
import matplotlib.pyplot as plt

#parse options
parser = argparse.ArgumentParser(description='Test the yambopy script.')
parser.add_argument('-dg','--doublegrid', action="store_true", help='Use double grid')
parser.add_argument('-c', '--calc', action="store_true", help='calculate the IP absorption')
parser.add_argument('-p', '--plot', action="store_true", help='plot the results')
args = parser.parse_args()

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

yambo  = 'yambo'
p2y    = 'p2y'
ypp    = 'ypp'
folder = 'ip'
prefix = 'bn'
scheduler = Scheduler.factory

#check if the SAVE folder is present
if not os.path.isdir('database/SAVE'):
    print('preparing yambo database')
    p2y_run = scheduler()
    p2y_run.add_command('mkdir -p database')
    p2y_run.add_command('cd nscf/bn.save; %s > %s.log'%(p2y,p2y))
    p2y_run.add_command('%s > %s.log'%(yambo,yambo))
    p2y_run.add_command('mv SAVE ../../database/')
    p2y_run.run()

if not os.path.islink('%s/SAVE'%folder):
    s = scheduler()
    s.add_command('mkdir -p %s'%folder)
    s.add_command('cd %s; ln -s ../database/SAVE .'%folder)
    if not args.doublegrid: s.add_command('cd .. ; rm -f database/SAVE/ndb.Double_Grid')
    s.run()

#initialize the double grid
if args.doublegrid and not os.path.isfile('database/SAVE/ndb.Double_Grid'):

    #check if the double grid nscf cycle is present
    if os.path.isdir('nscf_double/%s.save'%prefix):
        print('nscf_double calculation found!')
    else:
        print('nscf_double calculation not found!')
        exit()

    if not os.path.isdir('database_double/SAVE'):
        print('preparing yambo double database')
        shell = scheduler()
        shell.add_command('cd nscf_double/%s.save; %s; %s'%(prefix,p2y,yambo))
        shell.add_command('cd ../../')
        shell.add_command('mkdir -p database_double')
        shell.add_command('mv nscf_double/%s.save/SAVE database_double'%prefix)
        shell.run()

    #initialize the double grid
    print("creating double grid")

    yppin = YamboIn.from_runlevel('%s -m',filename='ypp.in',executable=ypp,folder='database')

    yppin['DbGd_DB1_paths'] = ["../database_double"]
    yppin.arguments.append('SkipCheck')

    yppin.write('database/ypp.in')

    shell = scheduler()
    shell.add_command('cd database; %s'%ypp)
    shell.add_command('cd ../%s ; rm -rf yambo o-*'%folder)
    #print(shell)
    shell.run()

if args.calc:
    #create the yambo input file
    y = YamboIn.from_runlevel('%s -o g -V all'%yambo,folder=folder)

    y['FFTGvecs'] = [30,'Ry']
    y['BndsRnXs'] = [1,30]
    y['QpntsRXd'] = [[1,1],'']
    y['ETStpsXd'] = 500

    y.write('%s/yambo_run.in'%folder)

    print('running yambo')
    yambo_run = scheduler()
    yambo_run.add_command('cd %s; %s -F yambo_run.in -J yambo'%(folder,yambo))
    if args.doublegrid: yambo_run.add_command('cd ..; rm database/SAVE/ndb.Double_Grid')
    yambo_run.run()

if args.plot:
    # Plot absorption spectrum
    data=np.genfromtxt('%s/o-yambo.eps_q1_ip'%folder,usecols=(0,1))
    fig = plt.figure(figsize=(4,5))
    ax = fig.add_axes( [ 0.20, 0.20, 0.70, 0.70 ])

    plt.plot(data[:,0],data[:,1],'-',c='b',label='IP Absorption')
    plt.legend()

    plt.show()
