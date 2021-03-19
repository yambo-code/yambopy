from yambopy import *
from schedulerpy import *
import time
import os
from copy import deepcopy
overflow = 1e8

def integerize(number):
    """
    Check if number is integer, if False make it integer
    """
    if isinstance(number,int): return number
    if isinstance(number,float):
        if number.is_integer(): return int(number)
        else: return integerize(number*10)

class YamboRTStep_Optimize():
    """ 
    Class to run convergence tests for the RT time step.

    Note: time steps must be given in as units.    

    - Needs an initialised RT SAVE
    - Needs an RT input 
    - Optional arguments: directory paths, max time step, time step increase, max number of runs
    - Optional: specify "yscheduler" as an instance of schedulerpy to run on clusters. This needs dynamical managing of submitted jobs.
    - Optional: add convergence loop with respect to field intensity

    Example of use:
        .. code-block:: python
    
            YamboRTStep_Optimize(input_path,SAVE_path,RUN_path,ref_time,TStep_MAX,TStep_increase,NSimulations)
    
        .. code-block:: python
    
            YamboRTStep_Optimize(input_path,SAVE_path,RUN_path,ref_time,TStep_MAX,TStep_increase,NSimulations,FieldInt)

    """

    def __init__(self,input_path='./yambo.in',SAVE_path='./SAVE',RUN_path='./RT_time-step_optimize',yambo_rt='yambo_rt',\
                 ref_time=60,TStep_MAX=30,TStep_increase=5,NSimulations=6,FieldInt=None,yscheduler=None,\
                 tol_eh=1e-4,tol_pol=5e-3,Tpoints_min=30,plot_results=True):
        #Configuring schedulers
        self.frontend = Scheduler.factory(scheduler="bash")
        if yscheduler is not None: #Here we use, e.g., slurm 
            self.jobrun = yscheduler
            self.wait_up = True #This will enable calls to a function that will make the code wait upon completion of previous submitted job
        else: 
            self.jobrun = Scheduler.factory(scheduler="bash")
            self.wait_up = False
        #Setting global variables
        input_path, input_name = input_path.rsplit('/',1)
        self.yin = YamboIn.from_file(filename=input_name,folder=input_path)
        self.RUN_path = RUN_path
        self.yambo_rt = yambo_rt

        self.ref_time = ref_time #Simulation duration (fs) after field ends.
        self.TStep_MAX = TStep_MAX
        self.TStep_increase = TStep_increase
        self.NSimulations = NSimulations
        self.Tpoints_min = Tpoints_min
        self.tol_eh = tol_eh
        self.tol_pol= tol_pol
        self.FieldInt = FieldInt
        self.time_odm = 1 # Important: this needs to be an integer
        #Generate directories
        self.create_folder_structure(SAVE_path)
        #Start IO
        if self.wait_up: prnt_io = False
        else: prnt_io = True
        self.yf = YamboIO(out_name='YAMBOPY_RTStepConvergence.log',out_path=self.RUN_path,print_to_shell=prnt_io)
        self.yf.IO_start()
        #
        if self.wait_up: self.yf.msg("The workflow is run using job submission through a scheduler.")
        else: self.yf.msg("The workflow is run locally.")
        #Check for consistent input parameters
        if integerize(self.TStep_MAX) % integerize(self.TStep_increase) !=0: #Here RaiseError may be used
            self.yf.msg("The polarization is computed at discrete times.")
            self.yf.msg("In order to compare efficiently results with different time steps,")
            self.yf.msg("please select a time increment that divides exactly the max time step.")
            self.yf.msg("Exiting...")
            exit() 
        if self.NSimulations < 2:
            self.yf.msg("NSimulations is too small to perform convergence tests.")
            self.yf.msg("Exiting...")
            exit()
        #Compute the dipoles, then prepare RT input
        self.FIND_values()
        self.COMPUTE_dipoles()
        #Run RT simulations and analyse data
        self.RUN_convergence()
        #[OPTIONAL] plot results
        if plot_results: self.PLOT_output()

        self.yf.IO_close()

    def create_folder_structure(self,SAVE_path):
        
        if not os.path.isdir(self.RUN_path):
            shell = self.frontend
            shell.add_command('mkdir -p %s'%self.RUN_path)
            shell.add_command('cd %s ; ln -s ../%s . ; cd ..'%(self.RUN_path,SAVE_path))
            shell.run()
            shell.clean()

        if not os.path.islink('%s/SAVE'%self.RUN_path):
            shell = self.frontend
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
                if self.yin['Field1_FWHM']==0.: # Here RaiseError may be used
                    self.yf.msg("Please use the variable Field1_FWHM to set field width (not Field1_kind)")
                    self.yf.msg("Exiting...")
                    exit()
            else:
                self.yf.msg("Please use the variable Field1_FWHM to set field width (not Field1_Width)")
                self.yf.msg("Exiting...")
                exit()
            fwhm = self.yin['Field1_FWHM'][0]
            sigma = fwhm/(2*np.sqrt(2*np.log(2))) #Standard deviation from FWHM from normal distr.
            sigma = np.around(sigma,2)
            self.yf.msg("with FWHM:    %f %s"%(fwhm,self.yin['Field1_FWHM'][1]))
            self.yf.msg("with dev.st.: %f %s"%(sigma,self.yin['Field1_FWHM'][1]))
            FieldTime = 6.*sigma      

        self.yf.msg("Field direction: %s"%(str(self.yin['Field1_Dir'][0])))

        #Set field intensity if given in input
        if self.FieldInt is not None: self.yin['Field1_Int']=[ self.FieldInt, 'kWLm2' ]
        else: self.FieldInt=self.yin['Field1_Int'][0] 

        #Set time steps
        time_steps = [ self.TStep_MAX - i*self.TStep_increase for i in range(self.NSimulations)]
        self.time_steps = np.array( [ ts for ts in time_steps if ts>0 ] )
        self.NSimulations = len(self.time_steps)
        self.TSteps_min_max=[self.TStep_MAX,self.TStep_MAX-(self.NSimulations-1)*self.TStep_increase]
        # We have to consider the case of non-integer below-as time steps
        if self.TStep_increase % 1 > 0: self.time_odm = 1000

        #Set simulations time settings (field time + lcm(time_steps) + hardcoded duration to analyse)
        ts_integerized = np.array( [ int('{0:g}'.format(ts*self.time_odm)) for ts in self.time_steps ] ) 
        # Grievous floating-point issues here for non-integer time steps... solved with above line
        ts_lcm = float(np.lcm.reduce(ts_integerized))/(1000.*self.time_odm) # in fs
        if self.ref_time/ts_lcm<self.Tpoints_min:
            self.yf.msg("[ERR] less than %d time points for polarization."%self.Tpoints_min)
            self.yf.msg("Exiting...")
            exit()
        self.yin['Field1_Tstart'] = [ts_lcm, 'fs']
        NETime = ts_lcm + FieldTime + self.ref_time
        self.yin['NETime'] = [ NETime, 'fs' ]
        self.NETime = NETime
        self.yf.msg("Total duration of simulations set to: %f fs"%NETime)
        self.yin['IOCachetime'] = [[ts_lcm,5*ts_lcm],'fs']

    def COMPUTE_dipoles(self,DIP_folder='dipoles'):
        """
        Compute the dipoles once and for all:
        In order for the dipoles to be compatible with a negf run 
        [a default optics run does not produce compatible dipoles], 
        the 'negf' argument is appended which causes the calculation to crash
        *after* the dipoles are computed.
        """
        if not os.path.isfile('%s/%s/ndb.dipoles'%(self.RUN_path,DIP_folder)):
            # Input...
            ydipoles = YamboIn()
            ydipoles.arguments.append('dipoles')
            ydipoles.arguments.append('negf')
            ydipoles['DIP_ROLEs'] = self.yin['DIP_ROLEs']
            ydipoles['DIP_CPU'] = self.yin['DIP_CPU']
            ydipoles['DipBands'] = self.yin['DipBands']
            ydipoles.write('%s/dipoles.in'%self.RUN_path)
            # Running...
            self.yf.msg("Running dipoles...")
            shell = deepcopy(self.jobrun)
            shell.name = 'dipoles'
            shell.add_mpirun_command('%s -F dipoles.in -J %s -C %s 2> %s.log'%(self.yambo_rt,DIP_folder,DIP_folder,DIP_folder))
            shell.run(filename='%s/rt.sh'%self.RUN_path)
            if self.wait_up: wait_for_job(shell,self.RUN_path) #OLD VERSION: self.wait_for_job(shell)
            shell.clean() 
        else:
            self.yf.msg("Dipoles found.")

        self.DIP_folder = DIP_folder

    def input_to_run(self,param,value,units):
        """
        Generate input for a specific run
        """
        yrun = deepcopy(self.yin)
        yrun[param] = [ value, units]
        return yrun

    def check_if_run_has_been_already_done(self,folder):
        """
        Skips running ts if:
        - ndb.output_polarization is found (sloppy check), AND
        - The string "Game over" appears in the report file (strict check)
        """
        skip1=False
        skip2=False
        ndb_pol_file='%s/%s/ndb.output_polarization'%(self.RUN_path,folder)
        report_file='%s/%s/r-%s_dipoles_negf'%(self.RUN_path,folder,folder)

        if os.path.isfile(ndb_pol_file): skip1=True
        if os.path.isfile(report_file):
            with open(report_file) as report:
                if 'Game Over' in report.read(): skip2=True

        return skip1*skip2

    def RUN_convergence(self,param='RTstep',units='as'):
        """
        Run the yambo_rt calculations flow.
        """        
        self.yf.msg("Running RT time step convergence...")
        RToutput  =    []
        NaN_check =    []
        eh_check  =    []
        pol_sq_check = []
        pol_x_check  = []
        time_steps = self.time_steps
        for i,ts in enumerate(time_steps):
            self.yf.msg("Running simulation for time step: %s %s"%('{0:g}'.format(ts),units))

            # Part 1: file preparation and run
            if self.time_odm==1: filename = '%s_%05d%s.in'%(param,ts,units)
            if self.time_odm==1000: filename = '%s_%05d%s.in'%(param,ts*self.time_odm,'zs')
            folder   = filename.split('.')[0]
            #Skip execution if output found:
            if self.check_if_run_has_been_already_done(folder):
                self.yf.msg("Found output for time step: %s %s"%('{0:g}'.format(ts),units))
            else:
                yrun = self.input_to_run(param,float(ts),units)
                yrun.write('%s/%s'%(self.RUN_path,filename))
                shell = deepcopy(self.jobrun)
                shell.name = '%s_%s_%s'%('{:.0E}'.format(self.FieldInt).replace("E+0", "E"),'{0:g}'.format(ts),shell.name)
                shell.add_mpirun_command('%s -F %s -J %s,%s -C %s 2> %s.log'%(self.yambo_rt,filename,folder,self.DIP_folder,folder,folder))
                shell.run(filename='%s/rt.sh'%self.RUN_path)
                if self.wait_up: wait_for_job(shell,self.RUN_path) #OLD VERSION: self.wait_for_job(shell)
                shell.clean()

            # Part 2: perform single-run analysis and store output
            out_dir = '%s/%s'%(self.RUN_path,folder)
            #Read output
            RTDB = YamboRTDB(calc=out_dir) #Read output
            RToutput_no_nan, NaN_test = self.nan_test(RTDB)              #[TEST1] NaN and overflow
            RToutput.append(RToutput_no_nan)
            if NaN_test: eh_test = self.electron_conservation_test(RTDB) #[TEST2] Electron number
            else:        eh_test = False
            NaN_check.append(NaN_test) 
            eh_check.append(eh_test) 

            # Part 3: perform polarization tests between subsequent runs
            if i==0: passed_counter = 0
            if i>0: 
                pol_sq_test, pol_x_test, passed_counter = self.ANALYSE_pol(RToutput,eh_check,passed_counter) #[TEST3],[TEST4] Polarization squared and along field direction
                pol_sq_check.append(pol_sq_test)
                pol_x_check.append(pol_x_test)

            # Part 4: decide if convergence was reached or we have to keep going
            if passed_counter==2:
                TStep_passed = self.time_steps[i-1]
                break
            
        if passed_counter==2: self.TStep_passed = TStep_passed
        if passed_counter==1: self.TStep_passed = self.time_steps[-1]
        if passed_counter==0: self.TStep_passed = None
            
        self.NSimulations = len(RToutput)
        self.RToutput = RToutput
        self.ANALYSE_output(NaN_check,eh_check,pol_sq_check,pol_x_check,passed_counter)

    def ANALYSE_output(self,NaN,eh,pol2,polx,passed,units='as'):
        """
        Output information and suggestion for an optimal time step.
        - There are two values of tolerance: one for carriers, one for polarization
        - Four increasingly stringent checks are performed:
            [1] NaN and overflow check to exclude botched runs
            [2] Conservation of electron number check
            [3] Error check of |pol|^2 (assuming lowest time step as reference)
            [4] Error check of pol along the field direction
        """
        self.yf.msg("---------- ANALYSIS ----------")

        NaN_passed = sum(NaN)
        self.yf.msg("[1] NaN and overflow test:")
        self.yf.msg("    Passed by %d out of %d."%(NaN_passed,self.NSimulations))

        eh_passed = sum(eh)
        self.yf.msg("[2] Conservation of electron number test (tol=%.0e):"%self.tol_eh)
        self.yf.msg("    Passed by %d out of %d."%(eh_passed,self.NSimulations))

        pol2_passed = sum(pol2)
        self.yf.msg("[3] Error in |pol|^2 test (tol=%.0e):"%self.tol_pol)
        self.yf.msg("    Passed by %d out of %d."%(pol2_passed,self.NSimulations-1))

        polx_passed = sum(polx)
        self.yf.msg("[4] Error in pol along field test (tol=%.0e):"%self.tol_pol)
        self.yf.msg("    Passed by %d out of %d."%(polx_passed,self.NSimulations-1))

        if passed == 1:
            self.yf.msg(" ")
            self.yf.msg("[WARNING] The lowest time step passed all the tests, but the")
            self.yf.msg("          additional safety run with a reduced step was")
            self.yf.msg("          not done due to NSimulations limit being reached.")

        if self.NSimulations == 2 or self.NSimulations == 3:
            self.yf.msg(" ")
            self.yf.msg("[WARNING] The largest time step already looks converged.")

        tp = self.TStep_passed
        self.yf.msg(" ")
        self.yf.msg("Based on the analysis, the suggested time step is: ")
        if tp is not None: self.yf.msg("### %s %s ###"%('{0:g}'.format(tp),units))
        else: self.yf.msg("[ERR] NSimulations limit reached before converged value was found.")
        self.yf.msg("------------------------------")        
 
    def ANALYSE_pol(self,RToutput,eh_check,passed):
        """
        Driver with the logical structure to manage polarization tests
        """
        if eh_check[-1]==True and eh_check[-2]==True:
            pol_sq_test = self.pol_error_test(RToutput,which_pol='pol_sq')
            pol_x_test  = self.pol_error_test(RToutput,which_pol='pol_along_field')
            
            if pol_sq_test and pol_x_test: passed = passed + 1
            else: passed = 0

        else: 
            pol_sq_test = False
            pol_x_test  = False

        return pol_sq_test, pol_x_test, passed

    def nan_test(self,RTDB):
        """ 
        Check computed polarizations for NaN values.
        """
        NaN_test = True
        # Check for NaN
        if np.isnan(RTDB.polarization).any() or np.isnan(RTDB.diff_carriers).any(): 
            RTDB.polarization = np.nan_to_num(RTDB.polarization) #Set to zero for plots
            NaN_test = False 
            #self.yf.msg("[WARNING] Yambo produced NaN values during this run")
        # Check for +/-Infinity
        if np.greater(np.abs(RTDB.polarization),overflow).any():
            RTDB.polarization[np.abs(RTDB.polarization)>overflow] = 0. #Set to zero for plots
            NaN_test = False
            #self.yf.msg("[WARNING] Yambo produced Infinity values during this run")           
 
        return RTDB, NaN_test

    def electron_conservation_test(self,RTDB):
        """
        Tests if elements of ratio_carriers are greater than tolerance.
        If any of them is, then the simulation in question has not passed the eh_test.
        """
        eh_carriers = np.greater(np.abs(RTDB.ratio_carriers),self.tol_eh)
        if any(eh_carriers): eh_test = False
        else:                eh_test = True
        return eh_test

    def pol_error_test(self,RTout,which_pol):
        """
        Computes the relative errors of the polarizations for each cached time.
        The cached times coincide for different runs.
        """
        pol_analyse= []
        pol_n1 = RTout[-1].polarization   
        pol_n0 = RTout[-2].polarization 
        if which_pol == 'pol_sq':  # Test for |pol|^2
            pol_analyse_n1 = pol_n1[0]*pol_n1[0] + pol_n1[1]*pol_n1[1] + pol_n1[2]*pol_n1[2] 
            pol_analyse_n0 = pol_n0[0]*pol_n0[0] + pol_n0[1]*pol_n0[1] + pol_n0[2]*pol_n0[2] 
        if which_pol == 'pol_along_field': # Test for pol along field
            pol_analyse_n1, _ = self.pol_along_field(pol_n1)
            pol_analyse_n0, _ = self.pol_along_field(pol_n0)            
 
        #Perform the test
        rel_err_pol = (pol_analyse_n1-pol_analyse_n0)/self.tol_pol
        error = np.greater(rel_err_pol,1.).any()
        if error: pol_test = False
        else:     pol_test = True      

        return pol_test
 
    def pol_along_field(self,pol):
        field = self.yin['Field1_Dir'][0]
        n_dirs = np.count_nonzero(field)
        pol_x = field[0]*pol[0]+field[1]*pol[1]+field[2]*pol[2]
        pol_x = pol_x/float(n_dirs)
        if field[0]!=0.:   pol_label='pol-x'
        elif field[1]!=0.: pol_label='pol-y'
        elif field[2]!=0.: pol_label='pol-z'
        else:              pol_label='pol-xyz'
        return pol_x,pol_label

    def PLOT_output(self,save_dir='plots'):
        """
        Generic plots generated by default, to be accessed by the user
        """
        import matplotlib.pyplot as plt

        plots_dir = 'runs'
        odm = self.time_odm
        
        if odm==1: units='as'
        if odm==1000: units='zs'
        for ts in self.time_steps: plots_dir += '_%s'%('{0:g}'.format(ts*odm))

        self.yf.msg("Plotting results.")
        out_dir = '%s/%s'%(self.RUN_path,save_dir)
        shell = self.frontend
        if not os.path.isdir(out_dir): shell.add_command('mkdir -p %s'%out_dir)
        shell.add_command('mkdir -p %s/%s'%(out_dir,plots_dir))
        shell.run()
        shell.clean()

        plots_dir = '%s/%s'%(out_dir,plots_dir)
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
             
            plt.savefig('%s/polarizations_%d%s.png'%(plots_dir,self.time_steps[ts]*odm,units),format='png',dpi=150)

        # Plot for all time steps
        f, (axes) = plt.subplots(4,1,sharex=True)
        for ts in range(self.NSimulations):

            label = '%d%s'%(time_steps[ts]*odm,units)
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

        plt.savefig('%s/polarizations_comparison.png'%plots_dir,format='png',dpi=150) 

        # Plot for all time steps |pol|^2
        f, (axes) = plt.subplots(self.NSimulations,1,sharex=True)
        for ts in range(self.NSimulations):

            pol = self.RToutput[ts].polarization
            pol_sq = pol[0]*pol[0] + pol[1]*pol[1] + pol[2]*pol[2]
            times = np.linspace(0.,self.NETime,num=pol.shape[1])
            pol_ts_label = "%d%s"%(time_steps[ts]*odm,units)
            axes[ts].plot(times, pol_sq, '-', lw=lwidth, color=ts_colors[ts], label=pol_ts_label)
        for ax in axes:
            ax.axhline(0.,lw=0.5,color='gray',zorder=-5)
            ax.legend(loc='upper left')
        f.tight_layout()

        plt.savefig('%s/polarizations_squared.png'%plots_dir,format='png',dpi=150)

        # Plot for all time steps along field direction
        f, (axes) = plt.subplots(self.NSimulations,1,sharex=True)
        for ts in range(self.NSimulations):
    
            pol = self.RToutput[ts].polarization
            pol_x,pol_label = self.pol_along_field(pol)
            times = np.linspace(0.,self.NETime,num=pol.shape[1])
            pol_ts_label = "%s_%d%s"%(pol_label,time_steps[ts]*odm,units)
            axes[ts].plot(times, pol_x, '-', lw=lwidth, color=ts_colors[ts], label=pol_ts_label) 
        for ax in axes:
            ax.axhline(0.,lw=0.5,color='gray',zorder=-5)
            ax.legend(loc='upper left')
        f.tight_layout()        

        plt.savefig('%s/polarizations_field_direction.png'%plots_dir,format='png',dpi=150)

    # OLD_VERSION: IN-WORKFLOW wait_for_job
    # def wait_for_job(self,shell,time_step=10.):
    #     """
    #     Let the python execution sleep until job completion
    #     """
    #     job_status = shell.check_job_status(self.RUN_path)
    #     condition = job_status=='R' or job_status=='PD' or job_status=='CG'
    #     while condition:
    #         time.sleep(time_step)
    #         job_status = shell.check_job_status(self.RUN_path) 
    #         condition = job_status=='R' or job_status=='PD' or job_status=='CG'
