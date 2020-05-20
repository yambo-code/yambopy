from yambopy import *
from schedulerpy import *
import os
overflow = 1e8

def int_round(value):
    return int(round(value))

def is_exactly_divisible(num,den,step=0):
    num_test = num - step
    if num_test % den ==0:
        return num_test,step
    else:
        step_new = step+1
        return is_exactly_divisible(num,den,step=step_new)

def relative_error(values,reference,tol):
    """
    Computes relative error (as deviation from reference) of values.
    """    
    err = 0.
    for iv,value in enumerate(values):
        err_tmp = ( (value-reference[iv])/tol )**2.
        err += err_tmp
    global_error = np.sqrt( err/len(values) )
    return global_error

class YamboRTStep_Optimize():
    """ 
    Class to run convergence tests for the RT time step.

    Note: time steps must be given in as units.    

    Example of use:

        .. code-block:: python
    
            YamboRTStep_Optimize(input_path,SAVE_path,RUN_path,TSteps_min_max,NSimulations)

        TO DO:
           (2) New tests for polarization
    """

    def __init__(self,input_path='./yambo.in',SAVE_path='./SAVE',RUN_path='./RT_time-step_optimize',yambo_rt='yambo_rt',TSteps_min_max=[5,50],NSimulations=5,db_manager=True,tol_eh=1e-4,tol_pol=1e-4):
    
        self.scheduler = Scheduler.factory
        input_path, input_name = input_path.rsplit('/',1)
        self.yin = YamboIn.from_file(filename=input_name,folder=input_path)
        self.RUN_path = RUN_path
        self.yambo_rt = yambo_rt

        self.ref_time = 30. #Simulation duration (fs) after field ends. HARDCODED.
        self.TSteps_min_max = TSteps_min_max
        self.NSimulations = NSimulations
        self.tol_eh = tol_eh
        self.tol_pol= tol_pol

        self.create_folder_structure(SAVE_path)

        self.yf = YamboIO(out_name='YAMBOPY_RTStepConvergence.log',out_path=self.RUN_path,print_to_shell=True)
        self.yf.IO_start()
        
        self.COMPUTE_dipoles()
        conv = self.FIND_values()
        self.RUN_convergence(conv)

        self.ANALYSE_output()
        self.PLOT_output()

        self.yf.IO_close()

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
            self.yf.msg("Field kind: DELTA")
            FieldTime = 0.

        if self.yin['Field1_kind']=="QSSIN":
            self.yf.msg("Field kind: QSSIN")
            if 'Field1_FWHM' in self.yin.variables.keys():
                if self.yin['Field1_FWHM']==0.:
                    self.yf.msg("Please use the variable Field1_FWHM to set field width (not Field1_kind)")
                    self.yf.msg("Exiting...")
                    exit()
            else:
                self.yf.msg("Please use the variable Field1_FWHM to set field width (not Field1_Width)")
                self.yf.msg("Exiting...")
                exit()
            self.yf.msg("with FWHM: %f %s"%(self.yin['Field1_FWHM'][0],self.yin['Field1_FWHM'][1]))
            FieldTime = 6.*self.yin['Field1_FWHM'][0]            

        self.yf.msg("Field direction: %s"%(str(self.yin['Field1_Dir'][0])))

        #Set simulations duration            
        NETime = FieldTime + self.ref_time 
        self.yin['NETime'] = [ NETime, 'fs' ]
        self.NETime = NETime
        self.yf.msg("Total duration of simulation set to: %f fs"%NETime)

        #Set time steps
        nstep = (self.TSteps_min_max[1]-self.TSteps_min_max[0])/(self.NSimulations-1.)
        self.time_steps=[int_round(self.TSteps_min_max[0]+nstep*i) for i in range(self.NSimulations)]
        if self.time_steps[-1] != self.TSteps_min_max[-1]: 
            self.TSteps_min_max[-1]==self.time_steps[-1]
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
            self.yf.msg("Running dipoles...")
            shell = self.scheduler()
            shell.add_command('cd %s'%self.RUN_path)
            #THIS must be replaced by a more advanced submission method
            shell.add_command('%s -F dipoles.in -J %s -C %s 2> %s.log'%(self.yambo_rt,DIP_folder,DIP_folder,DIP_folder))
            shell.run()
            shell.clean() 
        else:
            self.yf.msg("Dipoles found.")

        self.DIP_folder = DIP_folder

    def RUN_convergence(self,conv):
        """
        Run the yambo_rt calculations
        """
        self.yf.msg("Running RT time step convergence...")

        RToutput = []
        NaN_check = []
        def run(filename):
            """ Function to be called by the optimize function """
            folder = filename.split('.')[0]
            folder = folder + conv.get('RTstep')[1] #Add time step units
            self.yf.msg('%s %s'%(filename,folder))
            shell = self.scheduler()
            shell.add_command('cd %s'%self.RUN_path)
            #THIS must be replaced by a more advanced submission method
            shell.add_command('%s -F %s -J %s,%s -C %s 2> %s.log'%(self.yambo_rt,filename,folder,self.DIP_folder,folder,folder))
            shell.run()
            shell.clean()

            out_dir = '%s/%s'%(self.RUN_path,folder)
            #Read output
            RTDB = YamboRTDB(calc=out_dir)
            RToutput_no_nan, NaN_test = self.check_nan(RTDB)
            NaN_check.append(NaN_test)
            RToutput.append(RToutput_no_nan) 

        self.yf.msg("Running %d simulations for time steps from %d to %d as"%(self.NSimulations,self.TSteps_min_max[0],self.TSteps_min_max[1]))

        self.yin.optimize(conv,folder=self.RUN_path,run=run,ref_run=False)
        self.RToutput = RToutput
        self.NaN_check = NaN_check

    def check_nan(self,RTDB):
        """ 
        Check computed polarizations for NaN values.
        """
        NaN_test = True
        # Check for NaN
        if np.isnan(RTDB.polarization).any() or np.isnan(RTDB.diff_carriers).any(): 
            RTDB.polarization = np.nan_to_num(RTDB.polarization) #Set to zero for plots
            NaN_test = False 
            self.yf.msg("[WARNING] Yambo produced NaN values during this run")
        # Check for +/-Infinity
        if np.greater(np.abs(RTDB.polarization),overflow).any():
            RTDB.polarization[np.abs(RTDB.polarization)>overflow] = 0. #Set to zero for plots
            NaN_test = False
            self.yf.msg("[WARNING] Yambo produced Infinity values during this run")           
 
        return RTDB, NaN_test

    def ANALYSE_output(self):
        """
        Driver to analyse output and provide a suggestion for an optimal time step.
        - There are two values of tolerance: one for carriers, one for polarization
        - Four increasingly stringent checks are performed: 
            [1] NaN and overflow check to exclude botched runs
            [2] Conservation of electron number check 
            [3] Error check of |pol|^2 (assuming lowest time step as reference)
            [4] Error check of pol along the field direction
        """
        self.yf.msg("---------- ANALYSIS ----------")
        list_passed = [ts for ts,sim in enumerate(self.NaN_check) if sim]
        self.yf.msg("[1] NaN and overflow test:")
        self.yf.msg("    Passed by %d out of %d."%(len(list_passed),self.NSimulations))
        eh_check = self.electron_conservation_test(list_passed) 
        list_passed = [ts for ts,sim in enumerate(eh_check) if sim]
        self.yf.msg("[2] Conservation of electron number test (tol=%.0e):"%self.tol_eh)
        self.yf.msg("    Passed by %d out of %d."%(len(list_passed),self.NSimulations))
        pol_sq_check = self.pol_error_test(list_passed,which_pol='pol_sq')
        list_passed = [ts for ts,sim in enumerate(pol_sq_check) if sim]
        self.yf.msg("[3] Error in |pol|^2 test (tol=%.0e):"%self.tol_pol)
        self.yf.msg("    Passed by %d out of %d."%(len(list_passed),self.NSimulations))
        pol_x_check = self.pol_error_test(list_passed,which_pol='pol_along_field')
        list_passed = [ts for ts,sim in enumerate(pol_x_check) if sim]
        self.yf.msg("[4] Error in pol along field test (tol=%.0e):"%self.tol_pol)
        self.yf.msg("    Passed by %d out of %d."%(len(list_passed),self.NSimulations))
        self.yf.msg(" ")
        self.yf.msg("Based on the analysis, the suggested time step is: ")
        if len(list_passed)==0.:
            self.yf.msg("### NONE of the time step values considered passed the tests!")
            self.yf.msg("### Consider decreasing the time steps and trying again.")
        else:
            self.yf.msg("### %d as ###"%self.time_steps[list_passed[-1]])
        self.yf.msg("------------------------------")

    def electron_conservation_test(self,sims):
        """
        Tests if elements of ratio_carriers are greater than tolerance.
        If any of them is, then the simulation in question has not passed the eh_test.
        """
        eh_check = []
        for ts in sims:
            eh_test = True
            eh_carriers = np.greater(self.RToutput[ts].ratio_carriers,self.tol_eh)
            if any(eh_carriers): eh_test = False
            eh_check.append(eh_test)
        return eh_check

    def pol_error_test(self,sims,which_pol):
        """
        Bins the pol array into 5 fs intervals after the laser, computes the means
        of the various bins, and compares them to the means of the reference time step.
        The test is passed if all the means are within tolerance.
        """
        #Set up binning for different lasers, time steps
        ratio_laser = self.ref_time/self.NETime
        n_bins = int(self.ref_time/5) #[WARNING] self.ret_time must be in fs and divisible by 5
        bins_average = np.zeros([len(sims),n_bins])
        for i,ts in enumerate(sims):
            pol   = self.RToutput[ts].polarization    
            if which_pol == 'pol_sq':  # Test for |pol|^2
                pol_analyse = pol[0]*pol[0] + pol[1]*pol[1] + pol[2]*pol[2] 
            if which_pol == 'pol_along_field': # Test for pol along field
                dr, _ = self.pol_along_field()
                pol_analyse = pol[dr]
            l_tmp = len(pol_analyse)
            cut_pol = int_round((1.-ratio_laser)*l_tmp)
            l_pol = len(pol_analyse[cut_pol:])
            new_l, steps = is_exactly_divisible(l_pol,n_bins)

            #Do the binning and compute the mean
            bins_ts = pol_analyse[cut_pol+steps:].reshape(-1,n_bins)
            bins_average[i]=bins_ts.mean(axis=0)
 
        #Perform the test
        pol_check = []
        for i in range(len(sims)):
            error = relative_error(bins_average[i],bins_average[0],self.tol_pol)
            if error< 1.: pol_check.append(True)
            else:         pol_check.append(False)        

        return pol_check
 
    def pol_along_field(self):
        field = self.yin['Field1_Dir']
        if field[0]!=0.:   dr,pol_label=[0,'pol-x']
        elif field[1]!=0.: dr,pol_label=[0,'pol-y']
        elif field[2]!=0.: dr,pol_label=[0,'pol-z']
        else:              dr,pol_label=[0,'pol-x']
        return dr,pol_label

    def PLOT_output(self,save_dir='plots'):
        """
        Generic plots generated by default, to be accessed by the user
        """
        import matplotlib.pyplot as plt

        self.yf.msg("Plotting results.")
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
        for ts in range(self.NSimulations):
        
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
        for ts in range(self.NSimulations):

            label = '%das'%time_steps[ts]
            pol   = self.RToutput[ts].polarization
            pol_sq = pol[0]*pol[0] + pol[1]*pol[1] + pol[2]*pol[2]
            times = np.linspace(0.,self.NETime,num=pol.shape[1])
            axes[0].plot(times, pol[0], '-', lw=lwidth, color=ts_colors[ts], label=label)
            axes[1].plot(times, pol[1], '-', lw=lwidth, color=ts_colors[ts], label=label)
            axes[2].plot(times, pol[2], '-', lw=lwidth, color=ts_colors[ts], label=label)
            axes[3].plot(times, pol_sq, '-', lw=lwidth, color=ts_colors[ts], label=label)
            handles, labels = axes[3].get_legend_handles_labels()
        for ax in axes: ax.axhline(0.,lw=0.5,color='gray',zorder=-5)
        
        f.legend(handles, labels, loc='center right')
        f.tight_layout()

        plt.savefig('%s/polarizations_comparison.png'%out_dir,format='png',dpi=150) 

        # Plot for all time steps along field direction
        dr, pol_label = self.pol_along_field()
        f, (axes) = plt.subplots(self.NSimulations,1,sharex=True)
        for ts in range(self.NSimulations):
    
            pol = self.RToutput[ts].polarization
            times = np.linspace(0.,self.NETime,num=pol.shape[1])
            pol_ts_label = "%s_%das"%(pol_label,time_steps[ts])
            axes[ts].plot(times, pol[dr], '-', lw=lwidth, color=ts_colors[ts], label=pol_ts_label) 
        for ax in axes:
            ax.axhline(0.,lw=0.5,color='gray',zorder=-5)
            ax.legend(loc='upper left')
        f.tight_layout()        

        plt.savefig('%s/polarizations_field_direction.png'%out_dir,format='png',dpi=150)
            
