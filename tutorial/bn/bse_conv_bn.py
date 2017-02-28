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
        shell.add_command('cd nscf/%s.save; %s > p2y.log'%(prefix,p2y))
        shell.add_command('cd nscf/%s.save; %s > yambo.log'%(prefix,yambo))
        shell.add_command('mv nscf/%s.save/SAVE database'%prefix)
        shell.run()

    #create the folder to run the calculation
    if not os.path.isdir(folder):
        shell = scheduler()
        shell.add_command('mkdir -p %s'%folder)
        shell.add_command('cp -r database/SAVE %s/'%folder)
        shell.run()

def bse_convergence():
    #create the yambo input file
    y = YamboIn('%s -b -o b -k sex -y d -V all'%yambo,folder=folder)

    #common variables
    y['BSEBands'] = [4,5]
    y['BEnSteps'] = 500
    y['BEnRange'] = [[2.0,12.0],'eV']
    y['KfnQP_E']  = [2.91355133,1.0,1.0] #some scissor shift

    #list of variables to optimize and the values they might take
    conv = { 'FFTGvecs': [[10,15,20,30],'Ry'],
             'NGsBlkXs': [[1,2,3,5,6], 'Ry'],
             'BndsRnXs': [[1,10],[1,20],[1,30],[1,40]] }

    def run(filename):
        """
        Function to be called by the optimize function
        """
        path = filename.split('.')[0]
        print(filename, path)
        shell = scheduler()        
        shell.add_command('cd %s; yambo -F %s -J %s -C %s 2> %s.log'%(folder,filename,path,path,path))
        shell.run()

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
    parser.add_argument('-f', '--folder',     default="bse_run",    help='choose folder to put the results')
    parser.add_argument('-p', '--p2y',        default="store_true", help='p2y executable')
    parser.add_argument('-y', '--yambo',      default="store_true", help='yambo executable')

    args = parser.parse_args()
    folder = args.folder

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    create_save()
    if args.run:     bse_convergence()
    if args.analyse: analyse()
