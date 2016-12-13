#
# Author: Henrique Pereira Coutada Miranda
# Run a BSE calculation using yambo
#
from __future__ import print_function
import sys
from yambopy import *
from qepy import *
import argparse

#parse options
parser = argparse.ArgumentParser(description='Test the yambopy script.')
parser.add_argument('-dg','--doublegrid', action="store_true", help='Use double grid')
parser.add_argument('-r', '--run',        action="store_true", help='Run BSE calculation')
parser.add_argument('-a', '--analyse',    action="store_true", help='plot the results')
parser.add_argument('-f', '--folder',     default="bse_run",   help='choose folder where to put results')
args = parser.parse_args()
folder = args.folder

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

yambo = "yambo"

if not os.path.isdir('database'):
    os.mkdir('database')

#check if the nscf cycle is present
if os.path.isdir('nscf/bn.save'):
    print('nscf calculation found!')
else:
    print('nscf calculation not found!')
    exit()

#check if the SAVE folder is present
if not os.path.isdir('database/SAVE'):
    print('preparing yambo database')
    os.system('cd nscf/bn.save; p2y > p2y.log')
    os.system('cd nscf/bn.save; yambo > yambo.log')
    os.system('mv nscf/bn.save/SAVE database')

#check if the SAVE folder is present
if args.doublegrid:
    if not os.path.isdir('database_double/SAVE'):
        print('preparing yambo database')
        os.system('cd nscf_double/bn.save; p2y > p2y.log')
        os.system('cd nscf_double/bn.save; yambo > yambo.log')
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
    os.system('cd %s; ypp')

if args.run:
    #create the yambo input file
    y = YamboIn('%s -b -o b -k sex -y d -V all'%yambo,folder=folder)

    #common variables
    y['BSEBands'] = [4,5]
    y['BEnSteps'] = 500
    y['BEnRange'] = [[2.0,12.0],'eV']
    y['KfnQP_E']  = [2.91355133,1.0,1.0] #some scissor shift

    #list of variables to optimize and the values they might take
    conv = { 'FFTGvecs': [[10,15,20],'Ry'],
             'NGsBlkXs': [[1,2,3], 'Ry'],
             'BndsRnXs': [[1,10],[1,20],[1,30]] }

    def run(filename):
        """ Function to be called by the optimize function """
        path = filename.split('.')[0]
        print(filename, path)
        os.system('cd %s; yambo -F %s -J %s -C %s 2> %s.log'%(folder,filename,path,path,path))

    y.optimize(conv,run=run)

if args.analyse:
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
