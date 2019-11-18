#
# Author: Henrique Pereira Coutada Miranda
# Run a GW calculation using yambo
#
from __future__ import print_function
from builtins import range
from yambopy import *
from qepy import *
import argparse

yambo = "yambo"
p2y = "p2y"
folder='bse'

def doublegrid():
    global folder
    folder = "%s_dbg"%folder
    database()

    #check if the nscf cycle is present
    if os.path.isdir('nscf_double/mos2.save'):
        print('nscf_double calculation found!')
    else:
        print('nscf_double calculation not found!')
        exit()

    #check if the SAVE folder is present
    if not os.path.isdir('database_double/SAVE'):
        if not os.path.isdir('database_double'):
            os.mkdir('database_double')
        print('preparing yambo database')
        # we don't need to read the wavefunctions for the double grid
        os.system('cd nscf_double/mos2.save; %s -w > p2y.log'%p2y)
        os.system('cd nscf_double/mos2.save; %s > yambo.log'%yambo)
        os.system('mv nscf_double/mos2.save/SAVE database_double')

    #copy databases
    if not os.path.isdir(folder):
        os.mkdir(folder)
        os.system('cp -r database/SAVE %s'%folder)

    #initialize the double grid
    print("creating double grid")
    f = open('%s/ypp.in'%folder,'w')
    f.write("""kpts_map
    %DbGd_DB1_paths
    "../database_double"
    %""")
    f.close()
    os.system('cd %s; ypp'%folder)

def database():
    #check if the nscf cycle is present
    if os.path.isdir('nscf/mos2.save'):
        print('nscf calculation found!')
    else:
        print('nscf calculation not found!')

    #check if the SAVE folder is present
    if not os.path.isdir('database/SAVE'):
        if not os.path.isdir('database'):
            os.mkdir('database')
        print('preparing yambo database')
        # we don't need to read the wavefunctions for the double grid
        os.system('cd nscf/mos2.save; %s > p2y.log'%p2y)
        os.system('cd nscf/mos2.save; %s > yambo.log'%yambo)
        os.system('mv nscf/mos2.save/SAVE database')

    #copy databases
    if not os.path.isdir(folder):
        os.mkdir(folder)
        os.system('cp -r database/SAVE %s'%folder)

def run():
    database()

    #check if the SAVE folder is present
    if not os.path.isdir('database/SAVE'):
        if not os.path.isdir('database'):
            os.mkdir('database')
        print('preparing yambo database')
        os.system('cd nscf/mos2.save; %s > p2y.log'%p2y)
        os.system('cd nscf/mos2.save; %s > yambo.log'%yambo)
        os.system('mv nscf/mos2.save/SAVE database')

    #create the yambo input file
    y = YamboIn('%s -b -o b -k sex -y d -V all'%yambo,folder=folder)

    y['FFTGvecs'] = [20,'Ry']
    y['NGsBlkXs'] = [1,'Ry']
    y['BndsRnXs'] = [1,40]
    y['BSEBands'] = [8,11]
    y['BEnSteps'] = [500,'']
    y['BEnRange'] = [[0.0,6.0],'eV']

    y.arguments.append('WRbsWF')
    y.write('%s/yambo_run.in'%folder)

    print('running yambo')
    os.system('cd %s; %s -F yambo_run.in -J yambo'%(folder,yambo))

def analyse():
    #pack in a json file
    y = YamboOut('bse')
    y.pack()

    #get the absorption spectra
    a = YamboBSEAbsorptionSpectra('yambo',path='bse')
    excitons = a.get_excitons(min_intensity=0.5,max_energy=5,Degen_Step=0.001)
    print( "nexcitons: %d"%len(excitons) )
    print( "excitons:" )
    print( excitons )
    a.get_wavefunctions(Degen_Step=0.001,repx=list(range(-1,2)),repy=list(range(-1,2)),repz=list(range(1)),
                        Cells=[13,13,1],Hole=[0,0,9+.5], FFTGvecs=10,wf=True)
    a.write_json()

if __name__ == '__main__':
    #parse options
    parser = argparse.ArgumentParser(description='Test the yambopy script.')
    parser.add_argument('-r'  ,'--run',       action="store_true", help='Use double grid')
    parser.add_argument('-dg' ,'--doublegrid', action="store_true", help='Use double grid')
    parser.add_argument('-a', '--analyse',    action="store_true", help='plot the results')
    args = parser.parse_args()

    if args.doublegrid:
        doublegrid()        

    if args.run: run()
    if args.analyse: analyse()
