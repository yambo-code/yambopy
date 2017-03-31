#
# Author: Henrique Pereira Coutada Miranda
# Run a BSE calculation using yambo
#
from __future__ import print_function
import sys
from yambopy import *
from qepy import *
from schedulerpy import *
import argparse
import shutil

yambo = "yambo"
p2y = "p2y"
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

def bse_convergence(what='dielectric'):
    #create the yambo input file
    y = YamboIn('%s -b -o b -k sex -y d -V all'%yambo,folder=folder)

    #common variables
    y['BSEBands'] = [4,5]
    y['BEnSteps'] = 500
    y['BEnRange'] = [[2.0,12.0],'eV']
    y['KfnQP_E']  = [2.91355133,1.0,1.0] #some scissor shift
    y['DBsIOoff'] = 'BS' #turn off writting BSE hamiltonian DB (better performance)

    print(what)

    if what == 'dielectric':
        #list of variables to optimize the dielectric screening
        conv = { 'FFTGvecs': [[10,15,20,30],'Ry'],
                 'NGsBlkXs': [[1,2,3,5,6], 'Ry'],
                 'BndsRnXs': [[1,10],[1,20],[1,30],[1,40]] }
    else:
        # converged parameters for epsilon
        y['FFTGvecs'] = [30,'Ry'] 
        y['NGsBlkXs'] = [2,'Ry']
        y['BndsRnXs'] = [[1,40],'Ry']

        # choose a large number of Bands in BSE
        # BSEEhEny will choose which transitions to include
        y['BSEBands'] = [1,10]

        #list of variables to optimize the BSE
        conv = { 'BSEEhEny': [[[1,10],[1,11],[1,12]],'eV'],
                 'BSENGBlk': [[0,1,2], 'Ry'],
                 'BSENGexx': [[10,15,20],'Ry']}

    def run(filename):
        """
        Function to be called by the optimize function
        """
        path = filename.split('.')[0]
        print(filename, path)
        shell = scheduler()        
        shell.add_command('cd %s; yambo -F %s -J %s -C %s 2> %s.log'%(folder,filename,path,path,path))
        shell.run()
        
        # if the variables are not exclusive to the BSE
        # save the dielectric screening
        if True:
            if any([word not in filename for word in ['FFTGvecs','NGsBlkXs','BndsRnXs']]):
                #check if em1s was calculated
                em1s_dir = "%s/%s"%(folder,path) 
                em1s     = "%s/ndb.em1s"%(em1s_dir)
                if os.path.isfile(em1s):
                    #copy all the files
                    files = [f for f in os.listdir(em1s_dir) if 'ndb.em1s' in f]
                    for f in files:
                        shutil.copy("%s/%s"%(em1s_dir,f),'%s/SAVE/'%folder)

    y.optimize(conv,run=run)

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
        a = YamboBSEAbsorptionSpectra(path,path=folder)
        excitons = a.get_excitons(min_intensity=0.0005,max_energy=7,Degen_Step=0.01)
        print( "nexcitons: %d"%len(excitons) )
        print( "excitons:" )
        print( excitons )
        a.get_wavefunctions(Degen_Step=0.01,repx=range(-1,2),repy=range(-1,2),repz=range(1))
        a.write_json(path)

    #plot the results using yambo analyser
    y = YamboAnalyser(folder)
    print(y)
    y.plot_bse('eps')
    print('done!')


if __name__ == "__main__":
    #parse options
    parser = argparse.ArgumentParser(description='Test the yambopy script.')
    parser.add_argument('-r', '--run',        action="store_true",  help='run BSE convergence calculation')
    parser.add_argument('-a', '--analyse',    action="store_true",  help='plot the results')
    parser.add_argument('-e', '--epsilon',    action="store_true",  help='converge epsilon parameters')
    parser.add_argument('-b', '--bse',        action="store_true",  help='converge bse parameters')
    parser.add_argument('-f', '--folder',     default="bse_run",    help='choose folder to put the results')
    parser.add_argument('-p', '--p2y',        default="store_true", help='p2y executable')
    parser.add_argument('-y', '--yambo',      default="store_true", help='yambo executable')

    args = parser.parse_args()
    folder = args.folder

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    if args.bse:
        what = 'bse'
    else:
        what = 'dielectric'

    create_save()
    if args.run:     bse_convergence(what=what)
    if args.analyse: analyse()