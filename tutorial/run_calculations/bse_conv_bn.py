#
# Author: Henrique Pereira Coutada Miranda
# Run a BSE calculation using yambo
#
from __future__ import print_function
from builtins import range
import sys
from yambopy import *
from yambocommandline import*
from qepy import *
from schedulerpy import *
import argparse
import shutil
import matplotlib.pyplot as plt

yambo  = 'yambo'
p2y    = 'p2y'
prefix = 'bn'
folder = 'bse_conv'

scheduler = Scheduler.factory

def create_save():
    #check if the nscf cycle is present
    if os.path.isdir('nscf/%s.save'%prefix):
        print('nscf calculation found!')
    else:
        print('nscf calculation not found!')
        exit()

    #check if the SAVE folder is present
    if not os.path.isdir('database'):
        print('preparing yambo database')
        shell = scheduler()
        shell.add_command('cd nscf/%s.save; %s; %s'%(prefix,p2y,yambo))
        shell.add_command('cd ../../')
        shell.add_command('mkdir -p database')
        shell.add_command('mv nscf/%s.save/SAVE database'%prefix)
        shell.run()

    #create the folder to run the calculation
    if not os.path.isdir(folder):
        shell = scheduler()
        shell.add_command('mkdir -p %s'%folder)
        shell.add_command('cp -r database/SAVE %s/'%folder)
        shell.run()

def bse_convergence(what='screening',threads=1,nohup=False):
    if nohup: nohup = 'nohup'
    else:     nohup = ''

    #create the yambo input file
    y = YamboIn.from_runlevel('-X s -o b -k sex -y d -V all',executable=yambo,folder=folder)

    #default variables
    y['BSEBands'] = [4,5]
    y['BEnSteps'] = 500
    y['BSEEhEny'] = [[1.0,10.0],'eV']
    y['BEnRange'] = [[2.0,12.0],'eV']
    y['BSENGexx'] = [10,'Ry']
    y['KfnQP_E']  = [2.91355133,1.0,1.0] #some scissor shift
    y['DBsIOoff'] = 'BS' #turn off writting BSE hamiltonian DB (better performance)

    print(what)

    if what == 'screening':
        #list of variables to optimize the screening screening
        conv = { 'FFTGvecs': [[10,10,15,20,30],'Ry'],
                 'NGsBlkXs': [[1,1,2,3,5,6], 'Ry'],
                 'BndsRnXs': [[1,10],[1,10],[1,20],[1,30],[1,40]] }
    else:
        # converged parameters for epsilon
        y['FFTGvecs'] = [30,'Ry']
        y['NGsBlkXs'] = [2,'Ry']
        y['BndsRnXs'] = [[1,40],'Ry']

        # choose a large number of Bands in BSE
        # BSEEhEny will choose which transitions to include
        y['BSEBands'] = [1,10]

        #list of variables to optimize the BSE
        conv = { 'BSEEhEny': [[[1,10],[1,10],[1,12],[1,14]],'eV'],
                 'BSENGBlk': [[0,0,1,2], 'Ry'],
                 'BSENGexx': [[10,10,15,20],'Ry']}

    def run(filename):
        """
        Function to be called by the optimize function
        """
        path = filename.split('.')[0]
        print(filename, path)
        shell = scheduler()
        shell.add_command('cd %s'%folder)
        shell.add_command('%s mpirun -np %d %s -F %s -J %s -C %s 2> %s.log'%(nohup,threads,yambo,filename,path,path,path))
        shell.add_command('touch %s/done'%path)
        if not os.path.isfile("%s/%s/done"%(folder,path)):
            shell.run()

    y.optimize(conv,folder='bse_conv',run=run,ref_run=False)

def analyse():
    #pack the files in .json files
    pack_files_in_folder(folder)

    paths = []
    #get folder names
    for dirpath,dirnames,filenames in os.walk(folder):
        #ignore the root folder
        if dirpath == folder:
            continue

        #check if there are some output files in the folder
        if ([ f for f in filenames if 'o-' in f ]):
            paths.append( dirpath.split('/')[-1] )

    for path in paths:
        print( path )
        #get the absorption spectra
        a = YamboBSEAbsorptionSpectra(path)

    print( "To plot the data run:" )
    print( "python bse_conv_bn.py -p -e" )
    print( "python bse_conv_bn.py -p -b" )

def plot(what):
    #plot the results using yambo analyser
    y = YamboAnalyser(folder)
    print(y)

    fig = plt.figure(figsize=(10,8))
    if what == "screening":
        ax = plt.subplot(3,1,1)
        y.plot_bse(['eps','FFTGvecs'],ax=ax)
        ax = plt.subplot(3,1,2)
        y.plot_bse(['eps','NGsBlkXs'],ax=ax)
        ax = plt.subplot(3,1,3)
        y.plot_bse(['eps','BndsRnXs'],ax=ax)
        plt.tight_layout()
        #plt.savefig('screening_conv_plot.png')
        plt.show()
    else:
        ax = plt.subplot(3,1,1)
        y.plot_bse(['eps','BSEEhEny'],ax=ax)
        ax = plt.subplot(3,1,2)
        y.plot_bse(['eps','BSENGBlk'],ax=ax)
        ax = plt.subplot(3,1,3)
        y.plot_bse(['eps','BSENGexx'],ax=ax)
        plt.tight_layout()
        #plt.savefig('absorption_conv_plot.png')
        plt.show()
    print('done!')

if __name__ == "__main__":
    #parse options
    parser = argparse.ArgumentParser(description='Test the yambopy script.')
    parser.add_argument('-r', '--run',     action="store_true",  help='run BSE convergence calculation')
    parser.add_argument('-a', '--analyse', action="store_true",  help='analyse results data')
    parser.add_argument('-p', '--plot',    action="store_true",  help='plot the results')
    parser.add_argument('-e', '--epsilon', action="store_true",  help='converge epsilon parameters')
    parser.add_argument('-b', '--bse',     action="store_true",  help='converge bse parameters')
    parser.add_argument('-u', '--nohup',   action="store_true",  help='run the commands with nohup')
    parser.add_argument('-t', '--threads', default=1, type=int,  help='number of threads to use')

    args = parser.parse_args()
    threads = args.threads
    nohup = args.nohup

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    if args.bse:
        what = 'bse'
    else:
        what = 'screening'

    create_save()
    if args.run:     bse_convergence(what=what,threads=threads,nohup=nohup)
    if args.analyse: analyse()
    if args.plot:    plot(what)
