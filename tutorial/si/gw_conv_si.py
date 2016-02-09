#
# Author: Henrique Pereira Coutada Miranda
# Run a GW calculation using Yambo
#
from __future__ import print_function
from yambopy.inputfile import *
from yambopy.outputfile import *
from yambopy.analyse import *
from pwpy.inputfile import *
from pwpy.outputxml import *
import subprocess

yambo = 'yambo'

if not os.path.isdir('database'):
    os.mkdir('database')

#check if the nscf cycle is present
if os.path.isdir('nscf/si.save'):
    print('nscf calculation found!')
else:
    print('nscf calculation not found!')
    exit() 

#check if the SAVE folder is present
if not os.path.isdir('database/SAVE'):
    print('preparing yambo database')
    os.system('cd nscf/si.save; p2y')
    os.system('cd nscf/si.save; yambo')
    os.system('mv nscf/si.save/SAVE database')

if not os.path.isdir('gw_conv'):
    os.mkdir('gw_conv')
    os.system('cp -r database/SAVE gw_conv')

#create the yambo input file
y = YamboIn('%s -d -g n -V all'%yambo,folder='gw_conv')
y['QPkrange'][0][2:4] = [6,10]
conv = { 'FFTGvecs': [[10,15,20],'Ry'],
         'NGsBlkXd': [[1,2,5], 'Ry'],
         'BndsRnXd': [[1,10],[1,20],[1,30]] }

def run(filename):
    """ Function to be called by the optimize function """
    folder = filename.split('.')[0]
    print(filename,folder)
    os.system('cd gw_conv; yambo -F %s -J %s -C %s 2> %s.log'%(filename,folder,folder,folder))

y.optimize(conv,run=run)

#pack the files in .json files
for dirpath,dirnames,filenames in os.walk('gw_conv'):
    #check if there are some output files in the folder
    if ([ f for f in filenames if 'o-' in f ]):
        y = YamboOut(dirpath,save_folder='gw_conv')
        y.pack()

#plot the results using yambmo analyser
y = YamboAnalyser('gw_conv')
print(y)
y.plot_gw('qp')
print('done!')
