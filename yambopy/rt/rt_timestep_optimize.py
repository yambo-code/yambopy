from yambopy import *
from schedulerpy import *
import os

def int_round(value):
    return int(round(value))

class YamboRTStep_Optimize():
    """ 
    Class to run convergence tests for the RT time step.

    Note: time steps must be given in as units.    

    Example of use:

        .. code-block:: python
    
            YamboRTStep_Optimize(input_path,SAVE_path,RUN_path)

    SO FAR: creation of folder structure and running of the TD simulations
    TO DO: (1) output reading; DONE 
           (2) option to produce figures/plot for analysis in specific folders; DONE
           (3) database_manager: separate class that is called if db_manager is True; prints the __str__ function of this class (to be coded) along with timestamp in specific file that is checked. If printed info are present, calculations move to a different folder, if not, file is added. If db_manager is False, this doesn't happen. Database manager can also be called with 'retrieve_info' function to display contents of file (maybe with yamboin format)
           (3) calculation of optimal time step(s)
           (4) dynamic convergence runs
    """

    def __init__(self,input_path='./yambo.in',SAVE_path='./SAVE',RUN_path='./RT_time-step_optimize',yambo_rt='yambo_rt',TSteps_min_max=[5,25],NSimulations=5,db_manager=True):
    
        self.scheduler = Scheduler.factory
        input_path, input_name = input_path.rsplit('/',1)
        self.yin = YamboIn.from_file(filename=input_name,folder=input_path)
        self.RUN_path = RUN_path
        self.yambo_rt = yambo_rt

        self.ref_time = 30. #Simulation duration (fs) after field ends.
        self.TSteps_min_max = TSteps_min_max
        self.NSimulations = NSimulations

        self.create_folder_structure(SAVE_path)
        
        self.COMPUTE_dipoles()
        conv = self.FIND_values()
        self.RUN_convergence(conv)

        self.PLOT_output()

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
        Determine time step values to be run and simulation lengths.
        """

        #Check which laser is used
        if self.yin['Field1_kind']=="DELTA":
            print("Field kind: DELTA")
            FieldTime = 0.

        if self.yin['Field1_kind']=="QSSIN":
            print("Field kind: QSSIN")
            if 'Field1_FWHM' in self.yin.variables.keys():
                if self.yin['Field1_FWHM']==0.:
                    print("Please use the variable Field1_FWHM to set field width (not Field1_kind)")
                    print("Exiting...")
                    exit()
            else:
                print("Please use the variable Field1_FWHM to set field width (not Field1_Width)")
                print("Exiting...")
                exit()
            print("with FWHM: %f %s"%(self.yin['Field1_FWHM'][0],self.yin['Field1_FWHM'][1]))
            FieldTime = 6.*self.yin['Field1_FWHM'][0]            

        #Set simulations duration            
        NETime = FieldTime + self.ref_time 
        self.yin['NETime'] = [ NETime, 'fs' ]
        self.NETime = NETime
        print("Total duration of simulation set to: %f fs"%NETime)

        #Set time steps
        nstep = (self.TSteps_min_max[1]-self.TSteps_min_max[0])/(self.NSimulations-1.)
        self.time_steps=[int_round(self.TSteps_min_max[0]+nstep*i) for i in range(self.NSimulations)]
        if self.time_steps[-1] != self.TSteps_min_max[-1]: 
            self.TSteps_min_max[-1]==self.time_step[-1]
        conv = { 'RTstep': [ [self.time_steps[0]]+self.time_steps,'as'] }
 
        return conv

    def COMPUTE_dipoles(self,DIP_folder='dipoles'):
        """
        Compute the dipoles once and for all:
        In order for the dipoles to be compatible with a negf run 
        [a default optics run does not produce compatible dipoles], 
        the 'negf' argument is appended which causes the calculation to crash
        *after* the dipoles are computed.
        """
        if not os.path.isfile('%s/%s/ndb.dipoles'%(self.RUN_path,DIP_folder)):
            ydipoles = YamboIn()
            ydipoles.arguments.append('dipoles')
            ydipoles.arguments.append('negf')
            ydipoles['DIP_ROLEs'] = self.yin['DIP_ROLEs']
            ydipoles['DIP_CPU'] = self.yin['DIP_CPU']
            ydipoles['DipBands'] = self.yin['DipBands']
            ydipoles.write('%s/dipoles.in'%self.RUN_path)
            print("Running dipoles...")
            shell = self.scheduler()
            shell.add_command('cd %s'%self.RUN_path)
            #THIS must be replaced by a more advanced submission method
            shell.add_command('%s -F dipoles.in -J %s -C %s 2> %s.log'%(self.yambo_rt,DIP_folder,DIP_folder,DIP_folder))
            shell.run()
            shell.clean() 
        else:
            print("Dipoles found.")

        self.DIP_folder = DIP_folder

    def RUN_convergence(self,conv):
        """
        Run the yambo_rt calculations
        """
        print("Running RT time step convergence...")

        RToutput = []
        def run(filename):
            """ Function to be called by the optimize function """
            folder = filename.split('.')[0]
            folder = folder + conv.get('RTstep')[1] #Add time step units
            print(filename,folder)
            shell = self.scheduler()
            shell.add_command('cd %s'%self.RUN_path)
            #THIS must be replaced by a more advanced submission method
            shell.add_command('%s -F %s -J %s,%s -C %s 2> %s.log'%(self.yambo_rt,filename,folder,self.DIP_folder,folder,folder))
            shell.run()
            shell.clean()

            out_dir = '%s/%s'%(self.RUN_path,folder)
            RToutput.append(YamboRTDB(calc=out_dir)) #Read output

        print("Running %d simulations for time steps from %d to %d as"%(self.NSimulations,self.TSteps_min_max[0],self.TSteps_min_max[1]))

        self.yin.optimize(conv,folder=self.RUN_path,run=run,ref_run=False)
        self.RToutput = RToutput

    #def ANALYSE_output():
    #"""
    #Analyse carriers and polarization data to get optimal time step
    #"""
    def PLOT_output(self,save_dir='plots'):
        """
        Generic plots generated by default, to be accessed by the user
        """
        import matplotlib.pyplot as plt

        print("Plotting results.")
        out_dir = '%s/%s'%(self.RUN_path,save_dir)
        if not os.path.isdir(out_dir): 
            shell = self.scheduler()
            shell.add_command('mkdir -p %s'%out_dir)
            shell.run()
            shell.clean()

        time_steps = self.time_steps
        lwidth=0.8
        ts_colors = plt.cm.gist_rainbow(np.linspace(0.,1.,num=self.NSimulations))

        # Plot for each time step
        for ts in range(len(self.RToutput)):
        
            pol   = self.RToutput[ts].polarization
            pol_sq = pol[0]*pol[0] + pol[1]*pol[1] + pol[2]*pol[2]
            times = np.linspace(0.,self.NETime,num=pol.shape[1])
            f, (axes) = plt.subplots(4,1,sharex=True)
            axes[0].plot(times, pol[0], '-', lw=lwidth, color='blue',  label='pol-x')
            axes[1].plot(times, pol[1], '-', lw=lwidth, color='green', label='pol-y') 
            axes[2].plot(times, pol[2], '-', lw=lwidth, color='red',   label='pol-z')
            axes[3].plot(times, pol_sq, '-', lw=lwidth, color='orange',label='|pol|^2')
            for ax in axes:
                ax.axhline(0.,lw=0.5,color='gray',zorder=-5)
                ax.legend(loc='upper left')
            f.tight_layout()
             
            plt.savefig('%s/polarizations_%das.png'%(out_dir,self.time_steps[ts]),format='png',dpi=150)
        # Plot for all time steps
        f, (axes) = plt.subplots(4,1,sharex=True)
        for ts in range(len(self.RToutput)):

            pol   = self.RToutput[ts].polarization
            pol_sq = pol[0]*pol[0] + pol[1]*pol[1] + pol[2]*pol[2]
            times = np.linspace(0.,self.NETime,num=pol.shape[1])
            axes[0].plot(times, pol[0], '-', lw=lwidth, color=ts_colors[ts], label=time_steps[ts])
            axes[1].plot(times, pol[1], '-', lw=lwidth, color=ts_colors[ts], label=time_steps[ts])
            axes[2].plot(times, pol[2], '-', lw=lwidth, color=ts_colors[ts], label=time_steps[ts])
            axes[3].plot(times, pol_sq, '-', lw=lwidth, color=ts_colors[ts], label=time_steps[ts])
        for ax in axes:
            ax.axhline(0.,lw=0.5,color='gray',zorder=-5)
            #ax.legend(loc='upper left')
        f.tight_layout()

        plt.savefig('%s/polarizations_comparison.png'%out_dir,format='png',dpi=150) 
