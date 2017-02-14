#
#
# Tutorial Yambo School. April 2017 Lausanne 
# Convergence GW on hexagonal BN
#
#
from __future__ import print_function
import sys
from yambopy     import *
from qepy        import *
from schedulerpy import *
import argparse

yambo = 'yambo'
p2y = 'p2y'
prefix = 'bn'

# Scheduler ==>> bash
bash = Scheduler.factory

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
        shell = bash() 
        shell.add_command('mkdir -p database')
        shell.add_command('cd nscf/%s.save; %s; %s'%(prefix, p2y, yambo))
        shell.add_command('mv SAVE  ../../database/')
        shell.run()
        shell.clean()

def gw_convergence():
    #create the folder to run the calculation
    if not os.path.isdir('gw_conv'):
        shell = bash() 
        shell.add_command('mkdir -p gw_conv')
        shell.add_command('cp -r database/SAVE gw_conv/')
        shell.run()
        shell.clean()

    #create the yambo input file

    # GW calculation. Exact Dynamical Screening. Newton method
    y = YamboIn('%s -d -g n -V all'%yambo,folder='gw_conv')

    y['FFTGvecs'] = [2,'Ha']                     # Global Cutoff
    y['EXXRLvcs'] = [20,'Ha']                     # Self-energy. Exchange
    y['NGsBlkXd'] = [1,10]                       # Screening. Number of bands
    y['NGsBlkXd'] = [0,'mHa']                    # Cutoff Screening
    y['GbndRnge'] = [[1,10],'']                  # Self-energy. Number of bands
    y['QPkrange'] = [ [7,7,4,5], '']
    #y.arguments.append('ExtendOut')

    conv = { 'FFTGvecs': [[2,2,5,10,15,20],'Ha'],
             'NGsBlkXd': [[0,0,500,1000,1500,2000], 'mHa'],
             'BndsRnXd': [[[1,5],[1,10],[1,20],[1,30],[1,40],[1,50]],''] ,
             'GbndRnge': [[[1,5],[1,10],[1,20],[1,30],[1,40],[1,50]],''] }
  
             #'EXXRLvcs': [[2,2,5,10,15,20],'Ha'],
    def run(filename):
        """ Function to be called by the optimize function """
        folder = filename.split('.')[0]
        print(filename,folder)
        shell = bash() 
        shell.add_command('cd gw_conv; %s -F %s -J %s -C %s 2> %s.log'%(yambo,filename,folder,folder,folder))
        print(shell)
        shell.run()
        shell.clean()

    y.optimize(conv,run=run,ref_run=False)

def plot_convergence():
    #pack the files in .json files
    pack_files_in_folder('gw_conv')

    print('Select the converged value for each variable')

    shell = bash() 
    shell.add_command('python analyse_gw.py -bc 5 -kc 7 -bv 4 -kv 7 gw_conv FFTGvecs')
    shell.add_command('python analyse_gw.py -bc 5 -kc 7 -bv 4 -kv 7 gw_conv NGsBlkXd')
    shell.add_command('python analyse_gw.py -bc 5 -kc 7 -bv 4 -kv 7 gw_conv BndsRnXd')
    shell.add_command('python analyse_gw.py -bc 5 -kc 7 -bv 4 -kv 7 gw_conv GbndRnge')
    shell.run()
    shell.clean()

def gw():
    #create the folder to run the calculation
    if not os.path.isdir('gw'):
        shell = bash() 
        shell.add_command('mkdir -p gw')
        shell.add_command('cp -r database/SAVE gw/')
        shell.run()
        shell.clean()

    # GW calculation. Exact Dynamical Screening. Newton method
    y = YamboIn('%s -d -g n -V all'%yambo,folder='gw_conv')

    y['FFTGvecs'] = [2,'Ha']            # Global Cutoff
    y['EXXRLvcs'] = [20,'Ha']           # Self-energy. Exchange
    y['NGsBlkXd'] = [[1,40],'']         # Screening. Number of bands
    y['NGsBlkXd'] = [1500,'mHa']        # Cutoff Screening
    y['GbndRnge'] = [[1,20],'']         # Self-energy. Number of bands
    y['QPkrange'] = [ [1,7,4,5], '']
   
    y.write('gw/yambo_gw.in')

    shell = bash() 
    shell.add_command('cd gw; %s -F yambo_gw.in -J gw -C gw' % yambo)
    shell.run()
    shell.clean()

def plot_gw():
    #pack the files in .json files
    pack_files_in_folder('gw')

    #plot the results using yambm analyser
    ya = YamboAnalyser('gw')
    print('plot all qpoints')
    ya.plot_gw('qp')
    #ya.plot_gw('qp')

    print('plot along a path')
    path = [[[0,   0,   0],'$\Gamma$'],
            [[0.5, 0,   0],'M'],
            [[0.3333,0.3333, 0.0],'K'],
            [[0.0, 0.0, 0.0],'$\Gamma$']]
    ya.plot_gw_path('qp',path)

def xi():
    #create the folder to run the calculation
    if not os.path.isdir('gw'):
        shell = bash() 
        shell.add_command('mkdir -p gw')
        shell.add_command('cp -r database/SAVE gw/')
        shell.run()
        shell.clean()

    cohsex = YamboIn('%s -p c -g n -V all'%yambo,folder='gw')
    cohsex['FFTGvecs'] = [2,'Ha']            # Global Cutoff
    cohsex['EXXRLvcs'] = [20,'Ha']           # Self-energy. Exchange
    cohsex['NGsBlkXs'] = [[1,40],'']         # Screening. Number of bands
    cohsex['NGsBlkXs'] = [1500,'mHa']        # Cutoff Screening
    cohsex['GbndRnge'] = [[1,20],'']         # Self-energy. Number of bands
    cohsex['QPkrange'] = [ [1,7,4,5], '']
    cohsex.write('gw/yambo_cohsex.in')
    shell = bash() 
    shell.add_command('cd gw; %s -F yambo_cohsex.in -J cohsex' % yambo)
    shell.run()
    shell.clean()

    ppa = YamboIn('%s -p p -g n -V all'%yambo,folder='gw')
    ppa['FFTGvecs'] = [2,'Ha']            # Global Cutoff
    ppa['EXXRLvcs'] = [20,'Ha']           # Self-energy. Exchange
    ppa['NGsBlkXp'] = [[1,40],'']         # Screening. Number of bands
    ppa['NGsBlkXp'] = [1500,'mHa']        # Cutoff Screening
    ppa['GbndRnge'] = [[1,20],'']         # Self-energy. Number of bands
    ppa['QPkrange'] = [ [1,7,4,5], '']
    ppa.write('gw/yambo_ppa.in')
    shell = bash() 
    shell.add_command('cd gw; %s -F yambo_ppa.in -J ppa' % yambo)
    shell.run()
    shell.clean()

    ra = YamboIn('%s -d -g n -V all'%yambo,folder='gw')
    ra['FFTGvecs'] = [2,'Ha']            # Global Cutoff
    ra['EXXRLvcs'] = [20,'Ha']           # Self-energy. Exchange
    ra['NGsBlkXd'] = [[1,40],'']         # Screening. Number of bands
    ra['NGsBlkXd'] = [1500,'mHa']        # Cutoff Screening
    ra['GbndRnge'] = [[1,20],'']         # Self-energy. Number of bands
    ra['QPkrange'] = [ [1,7,4,5], '']
    ra.write('gw/yambo_ra.in')
    shell = bash() 
    shell.add_command('cd gw; %s -F yambo_ra.in -J ra' % yambo)
    shell.run()
    shell.clean()

if __name__ == "__main__":
    #parse options
    parser = argparse.ArgumentParser(description='GW convergence')
    parser.add_argument('-c' ,'--convergence',  action="store_true", help='Run convergence calculations')
    parser.add_argument('-p' ,'--plot', action="store_true",         help='Pack into json files and plot the convergence results')
    parser.add_argument('-g' ,'--gw', action="store_true",           help='Run a single GW calculation')
    parser.add_argument('-r' ,'--results', action="store_true",      help='Pack into json files and plot a single GW calculation')
    parser.add_argument('-x' ,'--xi', action="store_true",      help='Pack into json files and plot a single GW calculation')

    args = parser.parse_args()

    create_save()
    if args.convergence:    gw_convergence()
    if args.plot:           plot_convergence()
    if args.gw:             gw()
    if args.results:        plot_gw()
    if args.xi:            xi()
