import os
from yambopy import *
from schedulerpy import *

class YamboRTSetup():
    """
    Class to run the setup for RT calculations.

    Must be run outside the folder where the nscf calculation took place.    

    Example of use:

    Generate a SAVE file with reduced symmetries:

        .. code-block:: python
        
            YamboRTSetup(FIELD_direction,QE_prefix,nscf=nscf_path,database=save_path)
    
    Include electron-phonon matrix elements:
    
        .. code-block:: python
        
        YamboRTSetup(FIELD_direction,QE_prefix,nscf=nscf_path,database=save_path,elph=elph_path)
        

    TO DO: 
      - make it a command-line tool
      - add double grid support
    """
    def __init__(self,field_dir,prefix,nscf='nscf',database='database',MaxGvecs=None,elph_path=None,yambo_rt='yambo_rt',p2y='p2y',ypp='ypp',yambo_ph='yambo_ph',ypp_ph='ypp_ph'):

        self.scheduler = Scheduler.factory
        self.field_dir = field_dir
        self.MaxGvecs  = MaxGvecs
        self.prefix    = prefix
        self.yambo_rt  = yambo_rt
        self.yambo_ph  = yambo_ph 
        self.p2y       = p2y
        self.ypp       = ypp
        self.ypp_ph    = ypp_ph
        
        #Start IO
        self.yf = YamboIO(out_name='YAMBOPY_RTsetup.log',out_path=database,print_to_shell=True)
        self.yf.IO_start()

        self.initialize_SAVE(nscf,database)
        if elph_path is None: self.FixSymm(database)
        else: self.FixSymm_with_elph(database,elph_path) 
        
        self.yf.IO_close()

    def initialize_SAVE(self,nscf,database):
        """
        Generate SAVE folder from QE nscf calculation
        """
        qe_save = '%s/%s.save'%(nscf,self.prefix)
        #check if the nscf cycle is present
        if os.path.isdir(qe_save):
            self.yf.msg('nscf calculation found!')
        else:
            self.yf.msg('nscf calculation not found!')
            exit()

        #check if the SAVE folder is present
        if os.path.isdir('%s/SAVE'%database):
            self.yf.msg('SAVE database found!')
        if not os.path.isdir('%s/SAVE'%database):
            self.yf.msg('preparing yambo RT database')
            if os.path.isfile('%s/data-file.xml'%qe_save): qe_xml = 'data-file.xml'
            if os.path.isfile('%s/data-file-schema.xml'%qe_save): qe_xml = 'data-file-schema.xml'
            p2y_run = self.scheduler()
            p2y_run.add_command('mkdir -p %s'%database)
            p2y_run.add_command('cd %s; %s -F %s > p2y.log ; cd -'%(qe_save,self.p2y,qe_xml))
            p2y_run.add_command('cd %s; mv SAVE %s ; cd -'%(qe_save,database))
            p2y_run.run()

    def FixSymm(self,database):
        """
        Generate SAVE folder with reduced symmetries starting from original SAVE
        """
        filnm1 = 'setup.in'
        filnm2 = 'fixsymm.in'
        #check if symmetries have been removed
        if os.path.isdir('%s/FixSymm'%database):
            self.yf.msg('FixSymm folder found!')
        if not os.path.isdir('%s/FixSymm'%database):
            self.yf.msg('Removing symmetries')
            y1 = YamboIn.from_runlevel('-i -V RL',executable=self.yambo_rt,filename=filnm1,folder=database)
            if self.MaxGvecs is not None:
                y1['MaxGvecs'] = self.MaxGvecs
                y1.write('%s/%s'%(database,filnm1))
            yambort_run = self.scheduler()
            yambort_run.add_command('cd %s ; %s -F %s; cd -'%(database,self.yambo_rt,filnm1))
            yambort_run.run()

            y2 = YamboIn.from_runlevel('-y',executable=self.ypp,filename=filnm2,folder=database)
            y2['Efield1']=self.field_dir
            y2.arguments.append('RmTimeRev')
            y2.write('%s/%s'%(database,filnm2))

            ypp_run = self.scheduler()
            ypp_run.add_command('cd %s ; %s -F %s ; cd -'%(database,self.ypp,filnm2))
            ypp_run.add_command('cd %s/FixSymm ; %s ; cd -'%(database,self.yambo_rt))
            ypp_run.run()
    
    def FixSymm_with_elph(self,database,elph_path):
        """
        Generate SAVE folder with reduced symmetries starting from original SAVE
        and adding expanded gkkp matrix elements
        """
        filnm1 = 'setup.in'
        filnm2 = 'fixsymm.in'
        filnmph = 'gkkp.in'
        #check if symmetries have been removed
        if os.path.isdir('%s/FixSymm'%database):
            self.yf.msg('FixSymm folder found!')
        if not os.path.isdir('%s/FixSymm'%database):
        
            if os.path.isfile('%s/SAVE/ndb.elph_gkkp_expanded'):
                self.yf.msg('gkkp already expanded.')
            if not os.path.isfile('%s/SAVE/ndb.elph_gkkp_expanded'):
                self.yf.msg('Reading and expanding gkkp')
                y1 = YamboIn.from_runlevel('-i -V RL',executable=self.yambo_ph,filename=filnm1,folder=database)
                y1.arguments.append('BSEscatt')
                if self.MaxGvecs is not None:
                    y1['MaxGvecs'] = self.MaxGvecs
                y1.write('%s/%s'%(database,filnm1))  
                yamboph_run = self.scheduler()
                if not os.path.islink('elph_dir'): yamboph_run.add_command('cd %s ; ln -s %s/elph_dir . ; cd -'%(database,elph_path))
                yamboph_run.add_command('cd %s ; %s -F %s; cd -'%(database,self.yambo_ph,filnm1))
                yamboph_run.run()
            
                yph = YamboIn.from_runlevel('-gkkp',executable=self.ypp_ph,filename=filnmph,folder=database)
                yph.arguments.append('GkkpExpand')
                yph['DBsPATH'] = "'./elph_dir'"
                yph.write('%s/%s'%(database,filnmph))          
                yppph_run = self.scheduler()
                yppph_run.add_command('cd %s ; %s -F %s; cd -'%(database,self.ypp_ph,filnmph))
                yppph_run.run()            
            
            self.yf.msg('Removing symmetries')
            y2 = YamboIn.from_runlevel('-y',executable=self.ypp,filename=filnm2,folder=database)
            y2['Efield1']=self.field_dir
            y2.arguments.append('RmTimeRev')
            y2.write('%s/%s'%(database,filnm2))
            ypp_run = self.scheduler()
            ypp_run.add_command('cd %s ; %s -F %s ; cd -'%(database,self.ypp,filnm2))
            ypp_run.add_command('cd %s/FixSymm/SAVE ; cp ../../SAVE/ndb.elph_gkkp_expanded* . ; cd -'%database)
            ypp_run.add_command('cd %s/FixSymm ; %s ; cd -'%(database,self.yambo_rt))
            ypp_run.run()
