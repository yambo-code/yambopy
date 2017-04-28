#
# Author: Henrique Pereira Coutada Miranda
# Run a BSE calculation using yambo
# one job per q-point for the dielectric function
#
from __future__ import print_function
from builtins import map
from builtins import range
from yambopy import *
from qepy import *
import multiprocessing
import argparse
import sys

yambo = "yambo"
folder = "bse_par"

def databases():
    #check if the nscf cycle is present
    if os.path.isdir('nscf/bn.save'):
        print('nscf calculation found!')
    else:
        print('nscf calculation not found!')
        exit()

    #check if the SAVE folder is present
    if not os.path.isdir('database/SAVE'):
        if not os.path.isdir('database'):
            os.mkdir('database')
        print('preparing yambo database')
        os.system('cd nscf/bn.save; p2y')
        os.system('cd nscf/bn.save; yambo')
        os.system('mv nscf/bn.save/SAVE database')

    #check if the SAVE folder is present
    if args.doublegrid:
        if not os.path.isdir('database_double/SAVE'):
            print('preparing yambo database')
            os.system('cd nscf_double/bn.save; p2y')
            os.system('cd nscf_double/bn.save; yambo')
            os.system('mv nscf_double/bn.save/SAVE database_double')

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

def run_job(job):
    print(job)
    os.system(job)

def run(nthreads=1,cut=False):
    databases()

    #create the yambo input file
    y = YamboIn('yambo -r -b -o b -V all',folder=folder)

    if cut:
        y['CUTGeo'] = 'box z'
        y['CUTBox'] = [0,0,10]

        y['RandQpts'] = 1000000
        y['RandGvec'] = [1,'Ry']

    y['FFTGvecs'] = [30,'Ry']
    y['NGsBlkXs'] = [1,'Ry']
    y['BndsRnXs'] = [[1,30],'']
    y.write('%s/yambo_run.in'%folder)

    #get the number of q-points
    startk,endk = list(map(int,y['QpntsRXs'][0]))

    #prepare the q-points input files
    jobs = []
    for nk in range(1,endk+1):
        y['QpntsRXs'] = [[nk,nk],'']
        y.write('%s/yambo_q%d.in'%(folder,nk))
        if nk != 1:
            jobs.append('cd %s; %s -F yambo_q%d.in -J yambo_q%d -C yambo_q%d 2> log%d'%(folder,yambo,nk,nk,nk,nk))

    #calculate first q-point and dipoles
    os.system('cd %s; %s -F yambo_q1.in -J yambo_q1 -C yambo_q1'%(folder,yambo))
    #copy dipoles to save
    os.system('cp %s/yambo_q1/ndb.dip* %s/SAVE'%(folder,folder))

    p = multiprocessing.Pool(nthreads)
    p.map(run_job, jobs)

    #gather all the files
    if not os.path.isdir('%s/yambo'%folder):
        os.mkdir('%s/yambo'%folder)
    os.system('cp %s/yambo_q1/ndb.em* %s/yambo'%(folder,folder))
    os.system('cp %s/*/ndb.em*_fragment* %s/yambo'%(folder,folder))

    y = YamboIn('yambo -r -b -o b -k sex -y d -V all',folder=folder)
    y['FFTGvecs'] = [30,'Ry']
    y['NGsBlkXs'] = [1,'Ry']
    y['BndsRnXs'] = [[1,30],'']
    y['BSEBands'] = [[3,6],'']
    y['BEnSteps'] = [500,'']
    y['BEnRange'] = [[0.0,10.0],'eV']
    y['KfnQP_E']  = [2.91355133,1.0,1.0] #some scissor shift
    y.arguments.append('WRbsWF')
    y.write('%s/yambo_run.in'%folder)

    print('running yambo')
    os.system('cd %s; %s -F yambo_run.in -J yambo'%(folder,yambo))

def plot():
    #collect the data
    pack_files_in_folder(folder)

    #plot the results using yambo analyser
    y = YamboAnalyser()
    y.plot_bse( ['eps','diago'] )

if __name__ == "__main__":
    #parse options
    parser = argparse.ArgumentParser(description='Test the yambopy script.')
    parser.add_argument('-dg','--doublegrid', action="store_true", help='Use double grid')
    parser.add_argument('-r' ,'--run',        action="store_true", help='Run the calculation')
    parser.add_argument('-c' ,'--cut',     action="store_true", help='Use coulomb cutoff')
    parser.add_argument('-p' ,'--plot',     action="store_true", help='Run the analysis')
    parser.add_argument('-t' ,'--nthreads',   default=2,           help='Run the analysis')
    args = parser.parse_args()
    nthreads = int(args.nthreads)

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    cut = args.cut 
    if args.run:
        run(nthreads,cut)
    if args.plot:
        plot()

