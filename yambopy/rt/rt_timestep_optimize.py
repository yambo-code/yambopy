from yambopy import *
from schedulerpy import *
import os

class YamboRTStep_Optimize():
    """ 
    Class to run convergence tests for the RT time step.

    Example of use:

        .. code-block:: python
    
            RTStep_Optimize(input_path,SAVE_path,RUN_path)

    SO FAR: creation of folder structure and running of the TD simulations
    TO DO: (1) output reading; 
           (2) option to produce figures/plot for analysis in specific folders
           (3) calculation of optimal time step(s)
           (4) dynamic convergence runs
    """

    def __init__(self,input_path='./yambo.in',SAVE_path='./SAVE',RUN_path='./RT_time-step_optimize',yambo_rt='yambo_rt'):
    
        self.scheduler = Scheduler.factory
        input_path, input_name = input_path.rsplit('/',1)
        self.yin = YamboIn.from_file(filename=input_name,folder=input_path)
        self.RUN_path = RUN_path
        self.yambo_rt = yambo_rt

        self.create_folder_structure(SAVE_path)

        conv = self.FIND_values()
        self.RUN_convergence(conv)

    def create_folder_structure(self,SAVE_path):
        
        if not os.path.isdir(self.RUN_path):
            shell = self.scheduler()
            shell.add_command('mkdir -p %s'%self.RUN_path)
            shell.add_command('cd %s ; ln -s ../%s . ; cd ..'%(self.RUN_path,SAVE_path))
            shell.run()
            shell.clean()

        if not os.path.islink('%s/SAVE'%self.RUN_path):
            shell = self.scheduler()
            shell.add_command('cd %s ; ln -s ../%s . ; cd ..'%(self.RUN_path,SAVE_path))
            shell.run()
            shell.clean()

    def FIND_values(self):
        """ 
        Determine time step values to be run.
        """
        conv = { 'RTstep': [[1,5,10,11,12,15],'as'] } #Hardcoded. TO DO: dynamical.
        return conv

    def COMPUTE_dipoles(self):
        """
        Compute the dipoles once and for all
        """
        #Create input; take dipoles variables (i.e. bands and parall) from self.yin; run in folder

    def RUN_convergence(self,conv):
        """
        Run the yambo_rt calculations
        """

        def run(filename):
            """ Function to be called by the optimize function """
            folder = filename.split('.')[0]
            folder = folder + conv.get('RTstep')[1] #Add time step units
            print(filename,folder)
            shell = self.scheduler()
            shell.add_command('cd %s'%self.RUN_path)
            #THIS must be replaced by a more advanced submission method
            shell.add_command('%s -F %s -J %s -C %s 2> %s.log'%(self.yambo_rt,filename,folder,folder,folder))
            shell.run()
            shell.clean()

        #ADD an initial calculation of the DIPOLES in order not to recompute them (use -J "X/Y" within the runs
        self.yin.optimize(conv,folder=self.RUN_path,run=run,ref_run=False)
