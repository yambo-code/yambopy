#
# Author: Henrique Pereira Coutada Miranda
# Run a BSE calculation using yambo
# one job per q-point for the dielectric function
#
from __future__ import print_function
from yambopy import *
from pwpy import *
import argparse

#parse options
parser = argparse.ArgumentParser(description='Parallel BSE calculation')
parser.add_argument('-dg' ,'--doublegrid',  action="store_true", help='Run the calculation using a double grid')
parser.add_argument('-r' ,'--run',  action="store_true", help='Run the calculation')
parser.add_argument('-p' ,'--plot', action="store_true", help='Run the analysis')
args = parser.parse_args()

yambo = "yambo"
prefix = "mos2"

if not os.path.isdir('database'):
    os.mkdir('database')

#check if the nscf cycle is present
if os.path.isdir('nscf/%s.save'%prefix):
    print('nscf calculation found!')
else:
    print('nscf calculation not found!')
    exit()

#check if the SAVE folder is present
if not os.path.isdir('database/SAVE'):
    print('preparing yambo database')
    os.system('cd nscf/%s.save; p2y'%prefix)
    os.system('cd nscf/%s.save; yambo'%prefix)
    os.system('mv nscf/%s.save/SAVE database'%prefix)

#check if the SAVE folder is present
if not os.path.isdir('database_double/SAVE'):
    print('preparing yambo database2')
    os.system('cd nscf_double/%s.save; p2y'%prefix)
    os.system('cd nscf_double/%s.save; yambo'%prefix)
    os.system('mv nscf_double/%s.save/SAVE database_double'%prefix)

if not os.path.isdir('bse_par'):
    os.mkdir('bse_par')
    os.system('cp -r database/SAVE bse_par')

#initialize the double grid
if args.doublegrid:
    print("creating double grid")
    f = open('bse_par/ypp.in','w')
    f.write("""kpts_map
    %DbGd_DB1_paths
    "../database_double"
    %""")
    f.close()
    os.system('cd bse_par; ypp')

if args.run:
    #create the yambo input file
    y = YamboIn('yambo -r -b -o b -V all',folder='bse_par')

    y['FFTGvecs'] = [15,'Ry']
    y['NGsBlkXs'] = [1,'Ry']
    y['BndsRnXs'] = [[1,40],'']
    y.write('bse_par/yambo_run.in')

    #get the number of q-points
    _,nkpoints = y['QpntsRXs'][0]

    #prepare the q-points input files
    f = open('jobs.sh','w')
    for nk in xrange(1,int(nkpoints)+1):
        y['QpntsRXs'] = [[nk,nk],'']
        y.write('bse_par/yambo_q%d.in'%(nk))
        if nk != 1:
            f.write('cd bse_par; %s -F yambo_q%d.in -J yambo_q%d -C yambo_q%d\n'%(yambo,nk,nk,nk))
    f.close()

    #calculate first q-point and dipoles
    os.system('cd bse_par; %s -F yambo_q1.in -J yambo_q1 -C yambo_q1'%yambo)
    #copy dipoles to save
    os.system('cp bse_par/yambo_q1/ndb.dip* bse_par/SAVE')
    #run jobs using gnuparallel
    os.system('parallel :::: jobs.sh')

    #gather all the files
    os.system('cp merge_eps.py bse_par')
    os.system('cd bse_par; python merge_eps.py')

    y = YamboIn('yambo -r -b -o b -k sex -y d -V all',folder='bse_par')
    y['FFTGvecs'] = [15,'Ry']
    y['NGsBlkXs'] = [1,'Ry']
    y['BndsRnXs'] = [[1,40],'']
    y['BSEBands'] = [8,11]
    y['BEnSteps'] = 500
    y['BEnRange'] = [[0.0,6.0],'eV']
    y.arguments.append('WRbsWF')

    y.write('bse_par/yambo_run.in')
    os.system('cd bse_par; %s -F yambo_run.in -J yambo'%yambo)

if args.plot:
    #collect the data
    pack_files_in_folder('bse_par')

    #plot the results using yambo analyser
    y = YamboAnalyser()
    y.plot_bse('eps')
