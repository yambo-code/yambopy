import os
from yambopy import *
from yambocommandline.commands import generate_save
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
def get_q_k_size(database,elph_save):
    """
    Get number of k and q points in the (IBZ) grids
    """
    # kpoints
    db = Dataset(database+"/SAVE/ns.db1")
    Nk = len(db.variables['K-POINTS'][:].T)
    db.close()
    # qpoints
    Nq = len(glob('./elph_dir/s.dbph_0*'))
    return Nq,Nk

def generate_gkkp(database,qe_save,elph_save,y_dir,expand,scheduler):
    """
    Read gkkp from dfpt calculation
    """

    def run_ypp_ph(UseQindxB=False):
        """
        Run ypp_ph and do checks
        """
        yph = YamboIn.from_runlevel('-gkkp',executable=ypp_ph,filename=filnm2,folder=database)
        # Apparently, now gkkp_db must be always specified
        yph.arguments.append('gkkp_db')
        if expand:
            #if Nq!=Nk: yph.arguments.append('gkkp_db')
            yph.arguments.append('GkkpExpand')
            if UseQindxB: yph.arguments.append('UseQindxB')
            print('    expanding gkkp in the full BZ')
        yph['DBsPATH'] = "./elph_dir"
        if os.path.isfile('%s/s.dbph_bare_000001'%elph_save):
            print('    reading also bare gkkp')
            yph.arguments.append('GkkpReadBare')
        yph.write('%s/%s'%(database,filnm2)) 

        yppph_run = scheduler()
        yppph_run.add_command('cd %s ; %s -F %s; cd -'%(database,ypp_ph,filnm2))
        yppph_run.run()

    def dbs_are_not_there():
        """
        Check if elph databases were produced succesfully
        """
        return ( not os.path.isfile('%s/SAVE/ndb.elph_gkkp'%database) ) and \
               ( not os.path.isfile('%s/SAVE/ndb.elph_gkkp_expanded'%database) )

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
            y1['K_grids'] = "BSC"
            y1.write('%s/%s'%(database,filnm1))
            yamboph_run = scheduler()
            if not os.path.islink('%s/elph_dir'%database): yamboph_run.add_command('cd %s ; ln -s %s . ; cd -'%(database,elph_save))
            # Check if a custom grid is used
            Nk,Nq = get_q_k_size(database,elph_save)
            # Regular grid
            if Nk==Nq: yamboph_run.add_command('cd %s ; %s -F %s -J ./elph_dir ; cd -'%(database,yambo_ph,filnm1))
            # Custom grid
            else: yamboph_run.add_command('cd %s ; %s -F %s ; cd -'%(database,yambo_ph,filnm1))
            yamboph_run.run()

            # Run ypp_ph
            run_ypp_ph()
            
            dbs_are_not_there = ( not os.path.isfile('%s/SAVE/ndb.elph_gkkp'%database) ) and \
                                ( not os.path.isfile('%s/SAVE/ndb.elph_gkkp_expanded'%database) )
            if dbs_are_not_there:
                print("[WARNING] First attempt didn't work. Retrying with UseQindxB")
                run_ypp_ph(UseQindxB=True)
                #if dbs_are_not_there(): print('[ERROR] ndb.elph_gkkp databases not created. Check the logs.')
        
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
