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

yambo = "yambo"
folder = 'ip'
scheduler = Scheduler.factory

#check if the SAVE folder is present
if not os.path.isdir('database/SAVE'):
    print('preparing yambo database')
    p2y_run = scheduler()
    p2y_run.add_command('mkdir -p database')
    p2y_run.add_command('cd nscf/bn.save; p2y > p2y.log')
    p2y_run.add_command('yambo > yambo.log')
    p2y_run.add_command('mv SAVE ../../database/')
    p2y_run.run()

if not os.path.islink('%s/SAVE'%folder):
    s = scheduler()
    s.add_command('mkdir -p %s'%folder)
    s.add_command('cd %s; ln -s ../database/SAVE .'%folder)
    s.run()

#initialize the double grid
if args.doublegrid:
    print("creating double grid")
    f = open('%s/ypp.in'%folder,'w')
    f.write("""kpts_map
    %DbGd_DB1_paths
    "../database_double"
    %""")
    f.close()
    ypp_run = scheduler()
    ypp_run.add_command('cd %s; ypp'%folder)
    ypp_run.run()

if args.calc:
    #create the yambo input file
    y = YamboIn.from_runlevel('yambo -o g -V all',folder=folder)

    y['FFTGvecs'] = [30,'Ry']
    y['BndsRnXs'] = [1,30]
    y['QpntsRXd'] = [[1,1],'']
    y['ETStpsXd'] = 500

    y.write('%s/yambo_run.in'%folder)

    print('running yambo')
    yambo_run = scheduler()
    yambo_run.add_command('cd %s; %s -F yambo_run.in -J yambo'%(folder,yambo))
    yambo_run.run()

if args.plot:
    # Plot absorption spectrum
    data=np.genfromtxt('%s/o-yambo.eps_q1_ip'%folder,usecols=(0,1))
    fig = plt.figure(figsize=(4,5))
    ax = fig.add_axes( [ 0.20, 0.20, 0.70, 0.70 ])

    plt.plot(data[:,0],data[:,1],'-',c='b',label='IP Absorption')
    plt.legend()

    plt.show()
