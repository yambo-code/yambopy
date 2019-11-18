#
# Author: Henrique Pereira Coutada Miranda
# Run a BSE calculation using yambo
# one job per q-point for the dielectric function
#
from __future__ import print_function
from builtins import map
from builtins import range
from multiprocessing import Pool
from yambopy import *
from qepy import *
import multiprocessing
import argparse
import sys

yambo = "yambo"
prefix = "mos2"
folder = "bse_par"

def databases():
    #check if the nscf cycle is present
    if os.path.isdir('nscf/%s.save'%prefix):
        print('nscf calculation found!')
    else:
        print('nscf calculation not found!')
        exit()

    #check if the SAVE folder is present
    if not os.path.isdir('database/SAVE'):
        if not os.path.isdir('database'):
            os.mkdir('database')
        print('preparing yambo database')
        os.system('cd nscf/%s.save; p2y'%prefix)
        os.system('cd nscf/%s.save; %s'%(prefix,yambo))
        os.system('mv nscf/%s.save/SAVE database'%prefix)

    if not os.path.isdir(folder):
        os.mkdir(folder)
        os.system('cp -r database/SAVE %s'%folder)

def run_job(job):
    print(job)
    os.system(job)

def run(nthreads=1):
    databases()

    #create the yambo input file
    y = YamboIn('%s -r -b -o b -V all'%yambo,folder=folder)

    y['FFTGvecs'] = [15,'Ry']
    y['NGsBlkXs'] = [1,'Ry']
    y['BndsRnXs'] = [1,40]
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

    y = YamboIn('%s -r -b -o b -k sex -y d -V all'%yambo,folder=folder)
    y['FFTGvecs'] = [15,'Ry']
    y['NGsBlkXs'] = [1,'Ry']
    y['BndsRnXs'] = [1,40]
    y['BSEBands'] = [8,11]
    y['BEnSteps'] = 500
    y['BEnRange'] = [[0.0,6.0],'eV']
    y.arguments.append('WRbsWF')

    y.write('%s/yambo_run.in'%folder)
    os.system('cd %s; %s -F yambo_run.in -J yambo'%(folder,yambo))

def plot():
    #collect the data
    pack_files_in_folder(folder)

    #plot the results using yambo analyser
    y = YamboAnalyser()
    y.plot_bse('eps')

if __name__ == "__main__":
    #parse options
    parser = argparse.ArgumentParser(description='Parallel BSE calculation')
    parser.add_argument('-r' ,'--run',      action="store_true", help='Run the calculation')
    parser.add_argument('-p' ,'--plot',     action="store_true", help='Run the analysis')
    parser.add_argument('-t' ,'--nthreads', default=2, help='Run the analysis')
    args = parser.parse_args()
    nthreads = int(args.nthreads)

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    if args.run:
        run(nthreads)
    if args.plot:
        plot()
