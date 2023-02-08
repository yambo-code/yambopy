from __future__ import print_function
from yambopy import *
from schedulerpy import *
import argparse
import os

layer_separation = 12
folder = 'gw+bse'
yambo  = 'yambo'
p2y    = 'p2y'
prefix = 'bn'
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

def run(nthreads=1,cut=False):
    """
    run gw+bse calculation using yambo
    """
    y = YamboIn.from_runlevel('%s -X f -g n -V all'%yambo,folder=folder)

    if cut:
        y['CUTGeo'] = 'box z'
        y['CUTBox'] = [0,0,layer_separation-1]

    QPKrange,_ = y['QPkrange']
    startk,endk,startb,endb = QPKrange
    y['QPkrange'] = [startk,endk,3,6]
    y['FFTGvecs'] = [30,'Ry']
    y['EXXRLvcs'] = [80,'Ry']
    y['NGsBlkXd'] = [3,'Ry']
    y['BndsRnXd'] = [1,30]
    y['GbndRnge'] = [1,30]
    y.arguments.append('em1d')
    y.write('%s/yambo_gw.in'%folder)

    print('running gw')
    shell = scheduler()
    shell.add_command('cd %s; mpirun -np %d %s -F yambo_gw.in -J yambo'%(folder,nthreads,yambo))
    shell.run()

    #
    #create the bse input file
    y = YamboIn.from_runlevel('%s -X s -o b -k sex -y d -V all'%yambo,folder=folder)

    if cut:
        y['CUTGeo'] = 'box z'
        y['CUTBox'] = [0,0,layer_separation-1]

    y['FFTGvecs'] = [30,'Ry']
    y['NGsBlkXs'] = [3,'Ry']
    y['BndsRnXs'] = [1,30]
    y['BSEBands'] = [3,6]
    y['BEnSteps'] = 500
    y['BEnRange'] = [[0,8],'eV']
    y['KfnQPdb'] = 'E < yambo/ndb.QP' #Include previously computed quasiparticle energies
    y.arguments.append('WRbsWF')

    y.write('%s/yambo_bse.in'%folder)

    #run the bse calculation using the dielectric function from gw
    print('running bse')
    shell = scheduler()
    shell.add_command('cd %s; mpirun -np %d %s -F yambo_bse.in -J yambo'%(folder,nthreads,yambo))
    shell.run()

if __name__ == "__main__":

    #parse options
    parser = argparse.ArgumentParser(description='Run GW+BSE calculations on BN.')
    parser.add_argument('-r', '--run',        action="store_true", help='Run BSE calculation')
    parser.add_argument('-c', '--cut',        action="store_true", help='Use coulomb truncation')
    parser.add_argument('-t' ,'--nthreads',                        help='Number of threads', default=1)
   
    args = parser.parse_args()
    nthreads = int(args.nthreads)

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    cut = args.cut

    create_save()
    if args.run:    run(nthreads,cut)
    
