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
        
            RTSetup(FIELD_direction,QE_prefix,nscf=nscf_path,database=save_path)

    TO DO: make it a command-line tool
    """
    def __init__(self,field_dir,prefix,nscf='nscf',database='database',MaxGvecs=None,yambo_rt='yambo_rt',p2y='p2y',ypp='ypp'):

        self.scheduler = Scheduler.factory
        self.field_dir = field_dir
        self.MaxGvecs  = MaxGvecs
        self.prefix    = prefix
        self.yambo_rt  = yambo_rt
        self.p2y       = p2y
        self.ypp       = ypp

        self.initialize_SAVE(nscf,database)
        self.FixSymm(database)

    def initialize_SAVE(self,nscf,database):
        """
        Generate SAVE folder from QE nscf calculation
        """
        qe_save = '%s/%s.save'%(nscf,self.prefix)
        #check if the nscf cycle is present
        if os.path.isdir(qe_save):
            print('nscf calculation found!')
        else:
            print('nscf calculation not found!')
            exit()

        #check if the SAVE folder is present
        if os.path.isdir('%s/SAVE'%database):
            print('SAVE database found!')
        if not os.path.isdir('%s/SAVE'%database):
            print('preparing yambo RT database')
            if os.path.isfile('%s/data-file.xml'%qe_save): qe_xml = 'data-file.xml'
            if os.path.isfile('%s/data-file-schema.xml'%qe_save): qe_xml = 'data-file-schema.xml'
            p2y_run = self.scheduler()
            p2y_run.add_command('mkdir -p %s'%database)
            p2y_run.add_command('cd %s; %s -F %s > p2y.log'%(qe_save,self.p2y,qe_xml))
            p2y_run.add_command('mv SAVE ../../%s'%database)
            p2y_run.run()

    def FixSymm(self,database):
        """
        Generate SAVE folder with reduced symmetries starting from original SAVE
        """
        filnm1 = 'setup.in'
        filnm2 = 'fixsymm.in'
        #check if symmetries have been removed
        if os.path.isdir('%s/FixSymm'%database):
            print('FixSymm folder found!')
        if not os.path.isdir('%s/FixSymm'%database):
            print('Removing symmetries')
            y1 = YamboIn.from_runlevel('-i -V RL',executable=self.yambo_rt,filename=filnm1,folder=database)
            if self.MaxGvecs is not None:
                y1['MaxGvecs'] = self.MaxGvecs
                y1.write('%s/%s'%(database,filnm1))
            yambort_run = self.scheduler()
            yambort_run.add_command('cd %s ; %s -F %s; cd ../'%(database,self.yambo_rt,filnm1))
            yambort_run.run()

            y2 = YamboIn.from_runlevel('-y',executable=self.ypp,filename=filnm2,folder=database)
            y2['Efield1']=self.field_dir
            y2.arguments.append('RmTimeRev')
            y2.write('%s/%s'%(database,filnm2))

            ypp_run = self.scheduler()
            ypp_run.add_command('cd %s ; %s -F %s ; cd ../'%(database,self.ypp,filnm2))
            ypp_run.add_command('cd %s/FixSymm ; %s ; cd ../../'%(database,self.yambo_rt))
            ypp_run.run()
