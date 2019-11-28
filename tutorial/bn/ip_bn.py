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
    p2y_run.add_command('cd nscf/bn.save; yambo > yambo.log')
    p2y_run.add_command('mv nscf/bn.save/SAVE database')
    p2y_run.run()

if not os.path.isdir('%s/SAVE'%folder):
    s = scheduler()
    s.add_command(folder)
    s.add_command('cp -r database/SAVE %s'%folder)
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
    #pack in a json file
    y = YamboOut(folder)
    y.pack()
