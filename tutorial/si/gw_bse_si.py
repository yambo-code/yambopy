#
# Author: Henrique Pereira Coutada Miranda
# Run a GW+BSE calculation using Yambo
#
from __future__ import print_function
from yambopy import *
from qepy import *
from schedulerpy import *
import argparse
import sys

yambo =  'yambo'
p2y = 'p2y'
prefix = 'si'
folder = 'gw_bse'

# scheduler
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
        shell.add_command('pushd nscf/%s.save; %s; %s'%(prefix,p2y,yambo))
        shell.add_command('popd')
        shell.add_command('mkdir -p database')
        shell.add_command('mv nscf/%s.save/SAVE database'%prefix)
        shell.run()

    #create the folder to run the calculation
    if not os.path.isdir(folder):
        shell = scheduler()
        shell.add_command('mkdir -p %s'%folder)
        shell.add_command('cp -r database/SAVE %s/'%folder)
        shell.run()

def run(nthreads=1):
    """
    run gw+bse calculation using yambo
    """
    y = YamboIn('%s -d -p p -g n -V all'%yambo,folder=folder)
    QPKrange,_ = y['QPkrange']
    startk,endk,startb,endb = QPKrange
    y['QPkrange'] = [startk,endk,3,6]
    y['FFTGvecs'] = [30,'Ry']
    y['NGsBlkXs'] = [1,'Ry']
    y['BndsRnXs'] = [1,30]
    y.write('%s/yambo_run.in'%folder)

    print('running gw')
    shell = scheduler()
    shell.add_command('cd %s; mpirun -np %d %s -F yambo_run.in -J yambo'%(folder,nthreads,yambo))
    shell.run()

    #creathe the bse input file
    y = YamboIn('%s -b -o b -k sex -y d -V all'%yambo,folder=folder)
    y['FFTGvecs'] = [30,'Ry']
    y['NGsBlkXs'] = [1,'Ry']
    y['BndsRnXs'] = [1,30]
    y['BSEBands'] = [3,6]
    y['BEnSteps'] = 500
    y['BEnRange'] = [[0,8],'eV']
    y['KfnQPdb'] = 'E < yambo/ndb.QP'
    y.write('%s/yambo_run.in'%folder)

    #run the bse calculation using the dielectric function from gw
    print('running bse')
    shell = scheduler()
    shell.add_command('cd %s; mpirun -np %d %s -F yambo_run.in -J yambo'%(folder,nthreads,yambo))
    shell.run()

def analyse():
    """
    plot the absorption spectra
    """
    #read yambo output file
    yo = YamboOut(folder)
    yo.pack('%s/%s'%(folder,folder))

    #analyse the data
    ya = YamboAnalyser(folder)
    print(ya)

    print('plot BSE')
    ya.plot_bse('eel')
    ya.plot_bse('eps')

    print('plot GW')
    path = [[[0.5,   0,   0],'L'],
            [[  0,   0,   0],'$\Gamma$'],
            [[  0, 0.5, 0.5],'X'],
            [[1.0, 1.0, 1.0],'$\Gamma$']]
    ya.plot_gw('qp')
    ya.plot_gw_path('qp',path)

if __name__ == "__main__":

    #parse options
    parser = argparse.ArgumentParser(description='Run BSE calculations on BN.')
    parser.add_argument('-r', '--run',        action="store_true", help='Run BSE calculation')
    parser.add_argument('-c', '--cut',        action="store_true", help='Use coulomb truncation')
    parser.add_argument('-a', '--analyse',    action="store_true", help='plot the results')
    parser.add_argument('-t' ,'--nthreads',                        help='Number of threads', default=1)
    args = parser.parse_args()
    nthreads = int(args.nthreads)

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    cut = args.cut
    create_save()
    if args.run:     run(nthreads)
    if args.analyse: analyse()
