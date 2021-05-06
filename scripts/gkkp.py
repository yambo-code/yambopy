import os
from yambopy import *
from scripts import generate_save
from schedulerpy import *
import argparse

"""
Script to produce a yambo SAVE folder along with gkkp matrix elements.
Usable to quickly initialize a yambo_ph (equilibrium) calculation.

Inputs:
 1. --nscf_dir='path/to/nscf/save/folder' [OPTIONAL]
 2. --elph_dir='path/to/dfpt/elph_dir/folder'
 3. --yambo='path/to/yambo/executables' [OPTIONAL]
 4. --expand                            [OPTIONAL]

NB:
 The SAVE folder is created in the directory where the script is launched!

"""

def generate_gkkp(database,qe_save,elph_save,y_dir,expand,scheduler):
    """
    Read gkkp from dfpt calculation
    """
    # Generate SAVE folder if nscf_path is given
    if qe_save != "": generate_save.generate_save(database,qe_save,y_dir,scheduler,noinit=True)
    
    #check if gkkp databases are already present
    if os.path.isfile('%s/SAVE/ndb.elph_gkkp'%database) or os.path.isfile('%s/SAVE/ndb.elph_gkkp_expanded'%database):
        print('gkkp found!')
    else:
        #check if the elph_dir folder is present
        if not os.path.isfile('%s/s.dbph_000001'%elph_save):
            print('problem with dbph databases at %s'%elph_save)
            exit()
        else:
            print('reading gkkp')

            if y_dir!="": y_dir=y_dir+"/"
            yambo_ph = y_dir+"yambo_ph"
            ypp_ph = y_dir+"ypp_ph"
            filnm1 = 'setup.in'
            filnm2 = 'gkkp.in'

            y1 = YamboIn.from_runlevel('-i -V RL',executable=yambo_ph,filename=filnm1,folder=database)
            y1.arguments.append('BSEscatt')
            y1.write('%s/%s'%(database,filnm1))
            yamboph_run = scheduler()
            if not os.path.islink('%s/elph_dir'%database): yamboph_run.add_command('cd %s ; ln -s %s . ; cd -'%(database,elph_save))
            yamboph_run.add_command('cd %s ; %s -F %s -J ./elph_dir ; cd -'%(database,yambo_ph,filnm1))
            yamboph_run.run()

            yph = YamboIn.from_runlevel('-gkkp',executable=ypp_ph,filename=filnm2,folder=database)
            if expand:
                yph.arguments.append('GkkpExpand')
                print('    expanding gkkp in the full BZ')
            yph['DBsPATH'] = "./elph_dir"
            if os.path.isfile('%s/s.dbph_bare_000001'%elph_save):
                print('    reading also bare gkkp')
                yph.arguments.append('GkkpReadBare')
            yph.write('%s/%s'%(database,filnm2))

            yppph_run = scheduler()
            yppph_run.add_command('cd %s ; %s -F %s; cd -'%(database,ypp_ph,filnm2))
            yppph_run.run()
            if ( not os.path.isfile('%s/SAVE/ndb.elph_gkkp'%database) ) and ( not os.path.isfile('%s/SAVE/ndb.elph_gkkp_expanded'%database) ):
                print('[ERROR] ndb.elph_gkkp databases not created. Check the logs.')
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate SAVE folder including gkkp databases')
    parser.add_argument('-nscf','--nscf_dir', type=str, default="", help='<Optional> Path to nscf save folder')
    parser.add_argument('-elph','--elph_dir', type=str,help='<Required> Path to elph_dir folder',required=True)
    parser.add_argument('-y','--yambo_dir', type=str, default="", help='<Optional> Path to yambo executables')
    parser.add_argument('-e','--expand', action="store_true", help="Expand GKKP")
    args = parser.parse_args()

    nscf_dir = args.nscf_dir
    elph_dir = args.elph_dir
    yambo_dir = args.yambo_dir
    expand = args.expand

    database = './'
    scheduler = Scheduler.factory

    generate_gkkp(database,nscf_dir,elph_dir,yambo_dir,expand,scheduler)
