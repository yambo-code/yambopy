import os
from yambopy import *
from schedulerpy import *
import argparse

"""
Script to produce a yambo SAVE folder.
Usable to quickly initialize a standard yambo calculation.

Inputs:
 1. --nscf_dir='path/to/nscf/save/folder'
 2. --yambo='path/to/yambo/executables' [OPTIONAL]

NB:
 The SAVE folder is created in the directory where the script is launched!

"""

def generate_save(database,qe_save,y_dir,scheduler,noinit=False):
    """
    Generate SAVE folder from QE nscf calculation
    """
    #check if the nscf cycle is present
    if os.path.isdir(qe_save):
        print('nscf calculation found!')
    else:
        print('nscf calculation not found!')
        exit()

    if y_dir!="": y_dir=y_dir+"/"
    p2y = y_dir+"p2y"
    yambo = y_dir+"yambo"

    #check if the SAVE folder is present
    if os.path.isdir('%s/SAVE'%database):
        print('SAVE database found!')
    if not os.path.isdir('%s/SAVE'%database):
        print('preparing yambo database')

        p2y_run = scheduler()
        p2y_run.add_command('mkdir -p %s'%database)
        p2y_run.add_command('cd %s; %s > p2y.log ; cd -'%(qe_save,p2y))
        if not noinit: p2y_run.add_command('cd %s; %s > yambo.log ; cd -'%(qe_save,yambo))
        p2y_run.add_command('mv %s/SAVE %s'%(qe_save,database))
        p2y_run.run()
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate SAVE folder')
    parser.add_argument('-nscf','--nscf_dir', type=str,help='<Required> Path to nscf save folder', required=True)
    parser.add_argument('-y','--yambo_dir', type=str, default="", help='<Optional> Path to yambo executables')
    args = parser.parse_args()

    nscf_dir = args.nscf_dir
    yambo_dir = args.yambo_dir

    database = './'
    scheduler = Scheduler.factory

    generate_save(database,nscf_dir,yambo_dir,scheduler)

