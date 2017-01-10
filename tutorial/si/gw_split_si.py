#
# Author: Henrique Pereira Coutada Miranda
# Run a GW calculation using Yambo
#
from __future__ import print_function
import sys
from yambopy import *
from qepy import *
import argparse

yambo = 'yambo'
p2y = 'p2y'
prefix = 'si'
folder = 'gw_split'

def create_save():
    #check if the nscf cycle is present
    if os.path.isdir('nscf/%s.save'%prefix):
        print('nscf calculation found!')
    else:
        print('nscf calculation not found!')
        exit() 

    #check if the SAVE folder is present
    if not os.path.isdir('database'):
        os.mkdir('database')
    if not os.path.isdir('database/SAVE'):
        print('preparing yambo database')
        os.system('cd nscf/%s.save; %s'%(prefix,p2y))
        os.system('cd nscf/%s.save; %s'%(prefix,yambo))
        os.system('mv nscf/%s.save/SAVE database'%prefix)

def gw_convergence():
    #create the folder to run the calculation
    if not os.path.isdir(folder):
        os.mkdir(folder)
    if not os.path.isdir('%s/SAVE'%folder):
        os.system('cp -r database/SAVE %s'%folder)

    #create the yambo input file
    y = YamboIn('%s -p p -g n -V all'%yambo,folder=folder)
    #get the information about the kpoints and bands
    startk, endk, startb, endb = [int(i) for i in y['QPkrange'][0]]
    y['FFTGvecs'] = [20,'Ry']
    y['NGsBlkXp'] = [1, 'Ry']
    y['BndsRnXp'] = [1,20]
    y['GbndRnge'] = [1,20]
    conv = { 'QPkrange': [[nk,nk,4,5] for nk in xrange(endk+1)] }

    def run(filename):
        """ Function to be called by the optimize function """
        path = filename.split('.')[0]
        print(filename,path)
        os.system('cd %s; %s -F %s -C %s -J %s 2> %s.log'%(folder,yambo,filename,path,path,path))

    y.optimize(conv,run=run)

def plot_convergence():
    #pack the files in .json files
    pack_files_in_folder(folder)

    #plot the results using yambm analyser
    ya = YamboAnalyser(folder)
    print(ya)
    print('plot all qpoints')
    ya.plot_gw('qp')
    print('plot along a path')

    path = [[[0.5,   0,   0],'L'],
            [[  0,   0,   0],'$\Gamma$'],
            [[  0, 0.5, 0.5],'X'],
            [[1.0, 1.0, 1.0],'$\Gamma$']]
    ya.plot_gw_path('qp',path)

if __name__ == "__main__":
    #parse options
    parser = argparse.ArgumentParser(description='GW convergence')
    parser.add_argument('-r' ,'--run',  action="store_true", help='Run the calculation')
    parser.add_argument('-p' ,'--plot', action="store_true", help='Pack into json files and plot the results')
    args = parser.parse_args()

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    if not os.path.isdir('database'): os.mkdir('database')

    create_save()
    if args.run:    gw_convergence()
    if args.plot:   plot_convergence()
