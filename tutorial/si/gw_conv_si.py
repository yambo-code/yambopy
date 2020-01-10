#
# Run a GW convergence calculation using Yambo
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
    y = YamboIn.from_runlevel('%s -d -p p -g n -V all'%yambo,folder='gw_conv')
    y['GbndRnge'] = [[1,15],'']
    y['QPkrange'][0][2:4] = [2,6]
    conv = { 'FFTGvecs': [[1,1,2,5,10],'Ry'],
             'NGsBlkXp': [[0,0,500,1000,2000], 'mRy'],
             'BndsRnXp': [[1,10],[1,10],[1,20],[1,30],[1,40]] ,
             'GbndRnge': [[1,10],[1,10],[1,20],[1,30],[1,40]] }

    def run(filename):
        """ Function to be called by the optimize function """
        folder = filename.split('.')[0]
        print(filename,folder)
        os.system('cd gw_conv; %s -F %s -J %s -C %s 2> %s.log'%(yambo,filename,folder,folder,folder))

    y.optimize(conv,folder='gw_conv',run=run)

create_save()
gw_convergence()
