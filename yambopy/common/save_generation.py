import os
from yambopy import *
from schedulerpy import *

class Create_YamboSave():
    """
    Driver for the generation of the yambo SAVE folder

    Must be run outside the folder where the nscf calculation took place.

    It contains:                                                                  CASES:      
        - get_SAVE function to generate basic p2y+yambo SAVE -------------------> save_type='simple' [default]
        - logic to generate more advanced SAVEs as :

             -- fixsymm (for real-time calculations)  --------------------------> save_type='fixsymm' [turned on by field_dir]
             -- elph    (for electron-phonon calculations) ---------------------> save_type='elph'    [turned on by elph_path]
             -- combinations of the above, including:

                 --- expansion of elph matrix elements -------------------------> save_type='expanded_elph'
                 --- fixsymm with elph matrix elements -------------------------> save_type='fixsymm+elph'
                 --- reading of bare elph matrix elements (always when found)

    Example of use:

    TO BE FILLED

    Generate a SAVE file with reduced symmetries:

        .. code-block:: python
        
            YamboRTSetup(FIELD_direction,QE_prefix,nscf=nscf_path,database=save_path)
    
    Include electron-phonon matrix elements:
    
        .. code-block:: python
        
        YamboRTSetup(FIELD_direction,QE_prefix,nscf=nscf_path,database=save_path,elph=elph_path)

    """
    def __init__(self,prefix,nscf='./nscf',database='./database',save_type='simple',field_dir=None,elph_path=None,MaxGvecs=None,yambo_exec_path=''):

        list_of_possibilities = ['simple','elph','expanded_elph','fixsymm','fixsymm+elph']

        if not os.path.isdir(database): os.mkdir(database)

        #Global variables        
        self.scheduler = Scheduler.factory
        self.MaxGvecs  = MaxGvecs
        self.prefix    = prefix
        
        #Yambo executables
        if yambo_exec_path != '': yambo_exec_path+='/'
        self.yambo     = yambo_exec_path + 'yambo'
        self.yambo_rt  = yambo_exec_path + 'yambo_rt'
        self.yambo_ph  = yambo_exec_path + 'yambo_ph' 
        self.p2y       = yambo_exec_path + 'p2y'
        self.ypp       = yambo_exec_path + 'ypp'
        self.ypp_ph    = yambo_exec_path + 'ypp_ph'
        
        # Logic-important variables
        self.elph_path = elph_path
        self.field_dir = field_dir
        self.save_type = save_type
        
        #Start IO
        self.yf = YamboIO(out_name='YAMBOPY_SAVEsetup.log',out_path=database,print_to_shell=True)
        self.yf.IO_start()

        # [1] Logic to determine which SAVE to generate
        self.which_SAVE()
        # Up to this point we have save_type, field_dir and elph_path set correctly
        
        # [2] Generate SAVE
        if self.save_type=='simple' or self.save_type not in list_of_possibilities: 
            # [2a] Standard SAVE already initialized
            self.get_SAVE(nscf,database)
        else:
            # [2b] Standard SAVE not initialised
            self.get_SAVE(nscf,database,noinit=True)
            
            # [3] Additional operations
            if self.save_type=='elph'
                # [3a] SAVE including ndb.elph_gkkp databases
                self.get_SAVE_elph()
            elif self.save_type=='expanded_elph'
                # [3b] SAVE including expanded ndb.elph_gkkp_expanded databases
                self.get_SAVE_elph(expand=True)
            elif self.save_type=='fixsymm'
                # [3c] SAVE with symmetries removed according to an external field
                self.get_SAVE_fixsymm(database)
            elif self.save_type=='fixsymm+elph'
                # [3d] SAVE in cases [3c]+[3b]
                self.get_SAVE_fixsymm(database,with_elph=True)
            else: 
                raise ValueError('No suitable save_type found. See docstring for information on avalaible save types.')
        
        #End IO
        self.yf.IO_close()

    def which_SAVE():
        """
        Set logicals for execution
        """
        if self.save_type[:3]=='fix' and self.field_dir is None: 
            self.yf.msg('Field direction not specified, setting it to [1,0,0].')
            self.field_dir = [1,0,0]
            
        if self.field_dir is not None and self.save_type[:3]!='fix':
            self.yf.msg('External field found: fixsymm turned on.')
            self.save_type='fixsymm'
            
        if self.save_type[-4:]=='elph' and self.elph_path is None:
            try: os.isdir('./elph_dir') or os.isfile('./elph_dir')
            except: raise FileNotFoundError('elph_dir directory not specified.')
            
        if self.elph_path is not None and self.save_type[-4:]!='elph':
            self.yf.msg('elph_dir path found: elph turned on.')
            if self.save_type=='fixsymm': self.save_type='fixsymm+elph'
            else: self.save_type='elph'

    def get_SAVE(self,nscf,database,noinit=False):
        """
        Generate SAVE folder from QE nscf calculation
        """
        qe_save = '%s/%s.save'%(nscf,self.prefix)
        #check if the nscf cycle is present
        if os.path.isdir(qe_save): 
            self.yf.msg('nscf calculation found!')
        else: 
            self.yf.msg('nscf calculation not found!')
            raise FileNotFoundError('nscf calculation not found!')

        #check if the SAVE folder is present
        if os.path.isdir('%s/SAVE'%database):
            self.yf.msg('SAVE database found!')
            
        #generate SAVE
        if not os.path.isdir('%s/SAVE'%database):
            self.yf.msg('preparing yambo RT database')
            p2y_run = self.scheduler()
            p2y_run.add_command('mkdir -p %s'%database)
            p2y_run.add_command('cd %s; %s > p2y.log ; cd -'%(qe_save,self.p2y))
            if not noinit: p2y_run.add_command('cd %s; %s > yambo.log ; cd -'%(qe_save,self.yambo))
            p2y_run.add_command('cd %s; mv SAVE %s ; cd -'%(qe_save,database))
            p2y_run.run()

    def get_SAVE_elph(self,database,expand=False,filnm1='setup.in'):
        """
        Add gkkp from dfpt calculation to original SAVE folder
        """
        
        gkkp_dbs = 'ndb.elph_gkkp'
        if expand: gkkp_dbs+='_expanded'

        #check if gkkp databases are already present 
        if os.path.isfile('%s/SAVE/%s'%(database,gkkp_dbs)):
            self.yf.msg('gkkp already present.')
        else:
            #check if the elph_dir folder is present
            try: os.path.isfile('%s/elph_dir/s.dbph_000001'%self.elph_path)
            except: raise FileNotFoundError('problem with dbph databases at %s/elph_dir'%self.elph_path)
            else:
                self.yf.msg('Reading gkkp')
            
                filnm2 = 'gkkp.in'

                # Initialisation using yambo_ph
                y1 = YamboIn.from_runlevel('-i -V RL',executable=self.yambo_ph,filename=filnm1,folder=database)
                y1.arguments.append('BSEscatt')
                if self.MaxGvecs is not None: y1['MaxGvecs'] = self.MaxGvecs
                y1.write('%s/%s'%(database,filnm1))
                yamboph_run = scheduler()
                if not os.path.islink('%s/elph_dir'%database): yamboph_run.add_command('cd %s ; ln -s %s . ; cd -'%(database,self.elph_path))
                yamboph_run.add_command('cd %s ; %s -F %s -J ./elph_dir ; cd -'%(database,self.yambo_ph,filnm1))
                yamboph_run.run()

                # Reading and expansion of gkkp
                yph = YamboIn.from_runlevel('-gkkp',executable=self.ypp_ph,filename=filnm2,folder=database)
                if expand:
                    yph.arguments.append('GkkpExpand')
                    self.yf.msg('    expanding gkkp in the full BZ')
                yph['DBsPATH'] = "./elph_dir"
                
                # Check if bare gkkp are present
                if os.path.isfile('%s/elph_dir/s.dbph_bare_000001'%self.elph_path):
                    self.yf.msg('    reading also bare gkkp')
                    yph.arguments.append('GkkpReadBare')
                yph.write('%s/%s'%(database,filnm2))
                
                yppph_run = scheduler()
                yppph_run.add_command('cd %s ; %s -F %s; cd -'%(database,ypp_ph,filnm2))
                yppph_run.run()
                    
                if not os.path.isfile('%s/SAVE/%s'%(database,gkkp_dbs)):
                    self.yf.msg('[ERROR] %s databases not created. Check the logs.'%gkkp_dbs)
                    raise FileNotFoundError('%s databases not created. Check the logs.'%gkkp_dbs)

    def get_SAVE_fixsymm(self,database,with_elph=False):
        """
        Generate SAVE folder with reduced symmetries starting from original SAVE
        """
        filnm1 = 'setup.in'
        filnm2 = 'fixsymm.in'
        #check if symmetries have been removed
        if os.path.isdir('%s/FixSymm'%database):
            self.yf.msg('FixSymm folder found!')
        if not os.path.isdir('%s/FixSymm'%database):
            
            if not with_elph:
                # Initialisation using yambo_rt
                y1 = YamboIn.from_runlevel('-i -V RL',executable=self.yambo_rt,filename=filnm1,folder=database)
                if self.MaxGvecs is not None: y1['MaxGvecs'] = self.MaxGvecs
                y1.write('%s/%s'%(database,filnm1))
                yambort_run = self.scheduler()
                yambort_run.add_command('cd %s ; %s -F %s; cd -'%(database,self.yambo_rt,filnm1))
                yambort_run.run()
            else:
                # Get expanded gkkp
                self.get_SAVE_elph(database,expand=True,filnm1=filnm1)

            # Removing symmetries
            self.yf.msg('Removing symmetries')
            y2 = YamboIn.from_runlevel('-y',executable=self.ypp,filename=filnm2,folder=database)
            y2['Efield1']=self.field_dir
            y2.arguments.append('RmTimeRev')
            y2.write('%s/%s'%(database,filnm2))
            ypp_run = self.scheduler()
            ypp_run.add_command('cd %s ; %s -F %s ; cd -'%(database,self.ypp,filnm2))
            
            # Re-initialisation
            if not with_elph:
                ypp_run.add_command('cd %s/FixSymm ; %s ; cd -'%(database,self.yambo_rt))
            else:
                ypp_run.add_command('cd %s/FixSymm/SAVE ; mv ../../SAVE/ndb.elph_gkkp_expanded* . ; cd -'%database)
                ypp_run.add_command('cd %s/FixSymm ; cp ../%s . ; cd -'%(database,filnm1))
                ypp_run.add_command('cd %s/FixSymm ; %s -F %s -J SAVE ; cd -'%(database,self.yambo_ph,filnm1))
            ypp_run.run()
