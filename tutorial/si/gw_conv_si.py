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
    if not os.path.isdir('gw_conv'):
        os.mkdir('gw_conv')
    if not os.path.isdir('gw_conv/SAVE'):
        os.system('cp -r database/SAVE gw_conv')

    #create the yambo input file
    y = YamboIn('%s -p p -g n -V all'%yambo,folder='gw_conv')
    y['GbndRnge'] = [[1,15],'']
    y['QPkrange'][0][2:4] = [2,6]
    conv = { 'FFTGvecs': [[5,10,15],'Ry'],
             'NGsBlkXp': [[1,2,3], 'Ry'],
             'BndsRnXp': [[1,10],[1,20],[1,30]] }

    def run(filename):
        """ Function to be called by the optimize function """
        folder = filename.split('.')[0]
        print(filename,folder)
        os.system('cd gw_conv; %s -F %s -J %s -C %s 2> %s.log'%(yambo,filename,folder,folder,folder))

    y.optimize(conv,run=run)

def plot_convergence(show=True):
    #pack the files in .json files
    pack_files_in_folder('gw_conv')

    #plot the results using yambm analyser
    ya = YamboAnalyser('gw_conv')
    print(ya)
    print('plot all qpoints')
    ya.plot_gw(show=show)
    print('plot along a path')

    path = Path([ [[1.0,1.0,1.0],'G'],
               [[0.0,0.5,0.5],'X'],
               [[0.0,0.0,0.0],'G'],
               [[0.5,0.0,0.0],'L']], [20,20,20])
    ya.plot_gw(path=path,show=show)

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
