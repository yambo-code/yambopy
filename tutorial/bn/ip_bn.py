#
# Author: Henrique Pereira Coutada Miranda
# Run a IP calculation using yambo
#
from __future__ import print_function
import sys
from yambopy import *
from qepy import *
import argparse

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

#check if the SAVE folder is present
if not os.path.isdir('database/SAVE'):
    print('preparing yambo database')
    os.system('mkdir -p database')
    os.system('cd nscf/bn.save; p2y > p2y.log')
    os.system('cd nscf/bn.save; yambo > yambo.log')
    os.system('mv nscf/bn.save/SAVE database')

if not os.path.isdir(folder):
    os.mkdir(folder)
    os.system('cp -r database/SAVE %s'%folder)

#initialize the double grid
if args.doublegrid:
    print("creating double grid")
    f = open('%s/ypp.in'%folder,'w')
    f.write("""kpts_map
    %DbGd_DB1_paths
    "../database_double"
    %""")
    f.close()
    os.system('cd %s; ypp'%folder)

if args.calc:
    #create the yambo input file
    y = YamboIn('yambo -o g -V all',folder=folder)

    y['FFTGvecs'] = [30,'Ry']
    y['BndsRnXs'] = [1,30]
    y['QpntsRXd'] = [[1,1],'']
    y['ETStpsXd'] = 500
    
    y.write('%s/yambo_run.in'%folder)

    print('running yambo')
    os.system('cd %s; %s -F yambo_run.in -J yambo'%(folder,yambo))

if args.plot:
    #pack in a json file
    y = YamboOut(folder)
    y.pack()
