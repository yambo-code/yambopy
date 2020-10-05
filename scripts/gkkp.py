import os
from yambopy import *
from schedulerpy import *
import argparse

"""
Script to produce a yambo SAVE folder along with *unexpanded* gkkp matrix elements.
Usable to quickly initialize a yambo_ph (equilibrium) calculation.

Inputs:
 1. --nscf_dir='path/to/nscf/save/folder'
 2. --elph_dir='path/to/dfpt/elph_dir/folder'
 3. --yambo='path/to/yambo/executables' [OPTIONAL]

NB:
 The SAVE folder is created in the directory where the script is launched!

"""

def initialize_SAVE(database,qe_save,y_dir,scheduler,out):
    """
    Generate SAVE folder from QE nscf calculation
    """
    #check if the nscf cycle is present
    if os.path.isdir(qe_save):
        out.msg('nscf calculation found!')
    else:
        out.msg('nscf calculation not found!')
        exit()

    p2y = "%s/p2y"%y_dir
    yambo = "%s/yambo"%y_dir

    #check if the SAVE folder is present
    if os.path.isdir('%s/SAVE'%database):
        out.msg('SAVE database found!')
    if not os.path.isdir('%s/SAVE'%database):
        out.msg('preparing yambo database')

        p2y_run = scheduler()
        p2y_run.add_command('mkdir -p %s'%database)
        p2y_run.add_command('cd %s; %s > p2y.log ; cd -'%(qe_save,p2y))
        p2y_run.add_command('cd %s; %s > yambo.log ; cd -'%(qe_save,yambo))
        p2y_run.add_command('cd %s; mv SAVE %s ; cd -'%(qe_save,database))
        p2y_run.run()

def initialize_gkkp(database,elph_save,y_dir,scheduler):
    """
    Read gkkp from dfpt calculation
    """
    #check if gkkp databases are already present
    if os.path.isfile('%s/SAVE/ndb.elph_gkkp'%database) or os.path.isfile('%s/SAVE/ndb.elph_gkkp_expanded'%database):
        out.msg('gkkp found!')
    else:
        #check if the elph_dir folder is present
        if not os.path.isfile('%s/s.dbph_000001'):
            out.msg('problem with dbph databases at %s'%elph_dir)
            exit()
        else:
            out.msg('reading gkkp')

            yambo_ph = "%s/yambo_ph"%y_dir
            ypp_ph = "%s/ypp_ph"%y_dir
            filnm1 = 'setup.in'
            filnm2 = 'gkkp.in'

            y1 = YamboIn.from_runlevel('-i -V RL',executable=yambo_ph,filename=filnm1,folder=database)
            y1.arguments.append('BSEscatt')
            y1.write('%s/%s'%(database,filnm1))
            yamboph_run = scheduler()
            if not os.path.islink('%s/elph_dir'%database): yamboph_run.add_command('cd %s ; ln -s %s . ; cd -'%(database,elph_save))
            yamboph_run.add_command('cd %s ; %s -F %s -J ./elph_dir ; cd -'%(database,yambo_ph,filnm1))
            yamboph_run.run()

            yph = YamboIn.from_runlevel('-gkkp',executable=ypp_ph,filename=filnmph,folder=database)
            #yph.arguments.append('GkkpExpand')
            yph['DBsPATH'] = "./elph_dir"
            yph.write('%s/%s'%(database,filnm2))
            yppph_run = self.scheduler()

            yppph_run.add_command('cd %s ; %s -F %s; cd -'%(database,ypp_ph,filnm2))
            yppph_run.run()
            if not os.path.isfile('%s/SAVE/ndb.elph_gkkp'%database):
                out.msg('[ERROR] ndb.elph_gkkp databases not created. Check the logs.')
        
if __name__ == "__main__"
    parser = argparse.ArgumentParser(description='Generate SAVE folder including gkkp databases')
    parser.add_argument('-nscf','--nscf_dir', type=str,help='<Required> Path to nscf save folder', required=True)
    parser.add_argument('-elph','--elph_dir', type=str,help='<Required> Path to elph_dir folder',required=True)
    parser.add_argument('-y','--yambo_dir', type=str, default="", help='<Optional> Path to yambo executables')
    args = parser.parse_args()

    nscf_dir = args.nscf_dir
    elph_dir = args.elph_dir
    yambo_dir = args.yambo_dir

    database = './'
    scheduler = Scheduler.factory

    # Start IO
    yf = YamboIO(out_name='YAMBOPY_GKKPsetup.log',out_path=database,print_to_shell=True)
    yf.IO_start()

    initialize_SAVE(database,nscf_dir,yambo_dir,scheduler,yf)
    initialize_gkkp(database,elph_dir,yambo_dir,scheduler,yf)

    # End IO
    yf.IO_close()
