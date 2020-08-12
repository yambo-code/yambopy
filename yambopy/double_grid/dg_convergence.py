from yambopy import *
from qepy import *
from schedulerpy import *
import os
import subprocess

class YamboDG_Optimize():
    """ 
    Class to generate and run convergence tests for the yambo double grid. 

    - Needs a quantum espresso scf save folder
    - Needs nscf (qe) and independent-particles (ip, yambo) inputs
    - Needs lists of coarse grids (CG) and fine grids (FG); [NOTE] Only random FG presently implemented. 
    - Additional arguments: directory paths, file names, experimental laser energy, etc.
    - The workflow is divided in FOUR STEPS that can be executed separately:
        1. nscf CG [STEPS='1']
        2. nscf FG and ip CG [STEPS='2']
        3. ip FG [STEPS='3']
        4. plot results [STEPS='4']

    Example of use:
        .. code-block:: python
    
            YamboDG_Optimize(cg_grids,fg_grids,prefix,scf_path,pseudo_path,...,STEP='all')

    """
    
    def __init__(self,cg_grids,fg_grids,prefix,nscf_input,ip_input,scf_save_path,pseudo_path,RUN_path='./',nscf_out="nscf",y_out_dir="results",qe_run_script=None,yambo_run_script=None,E_laser=0.,pw='pw.x',yambo='yambo',ypp='ypp',p2y='p2y',STEPS='all'):
        
        #Setting global variables
        self.scheduler = Scheduler.factory
        self.cg_grids = cg_grids
        self.cg_strings = []
        for cg in cg_grids: self.cg_strings.append("%dx%dx%d"%(cg[0],cg[1],cg[2]))
        self.fg_grids = fg_grids
        self.fg_strings = []
        for fg in self.fg_grids:
                temp_ls = []
                for i in range(len(fg)): temp_ls.append(str(fg[i])+'_fg')
                self.fg_strings.append(temp_ls)
        self.prefix = prefix
        self.scf_save_path = scf_save_path
        self.pseudo_path = pseudo_path
        self.RUN_path = RUN_path
        
        # Path of nscf and ip calculations and final plots
        self.nscf_dir = '%s/nscf_grids'%RUN_path
        self.ip_dir = '%s/ip_grids'%RUN_path
        self.plot_dir = '%s/plots'%RUN_path
        if not os.path.isdir(self.nscf_dir): os.mkdir(self.nscf_dir)
        if not os.path.isdir(self.ip_dir): os.mkdir(self.ip_dir)
        
        #Inputs
        self.nscf_inp = nscf_input
        self.ip_inp   = ip_input

        #Start IO
        self.yf = YamboIO(out_name='YAMBOPY_double-grid_Optimize.log',out_path=self.RUN_path,print_to_shell=True)
        self.yf.IO_start()

        self.yf.msg("So far so good.")
        
        self.yf.IO_close()
        
    # def create_folder_structure(self,SAVE_path):
    # 
    #     if not os.path.isdir(self.RUN_path):
    #         shell = self.scheduler()
    #         shell.add_command('mkdir -p %s'%self.RUN_path)
    #         shell.add_command('cd %s ; ln -s ../%s . ; cd ..'%(self.RUN_path,SAVE_path))
    #         shell.run()
    #         shell.clean()
    # 
    #     if not os.path.islink('%s/SAVE'%self.RUN_path):
    #         shell = self.scheduler()
    #         shell.add_command('cd %s ; ln -s ../%s . ; cd ..'%(self.RUN_path,SAVE_path))
    #         shell.run()
    #         shell.clean()
    # 
    # def FIND_values(self):
    #     """ 
    #     Determine time step values to be run and simulation lengths.
    #     """
    # 
    #     #Check which laser is used
    #     if self.yin['Field1_kind']=="DELTA":
    #         self.yf.msg("Field kind: DELTA")
    #         FieldTime = 0.
    # 
    #     if self.yin['Field1_kind']=="QSSIN":
    #         self.yf.msg("Field kind: QSSIN")
    #         if 'Field1_FWHM' in self.yin.variables.keys():
    #             if self.yin['Field1_FWHM']==0.: # Here RaiseError may be used
    #                 self.yf.msg("Please use the variable Field1_FWHM to set field width (not Field1_kind)")
    #                 self.yf.msg("Exiting...")
    #                 exit()
    #         else:
    #             self.yf.msg("Please use the variable Field1_FWHM to set field width (not Field1_Width)")
    #             self.yf.msg("Exiting...")
    #             exit()
    #         fwhm = self.yin['Field1_FWHM'][0]
    #         sigma = fwhm/(2*np.sqrt(2*np.log(2))) #Standard deviation from FWHM from normal distr.
    #         sigma = np.around(sigma,2)
    #         self.yf.msg("with FWHM:    %f %s"%(fwhm,self.yin['Field1_FWHM'][1]))
    #         self.yf.msg("with dev.st.: %f %s"%(sigma,self.yin['Field1_FWHM'][1]))
    #         FieldTime = 6.*sigma      
    # 
    #     self.yf.msg("Field direction: %s"%(str(self.yin['Field1_Dir'][0])))
    # 
    #     #Set time steps
    #     time_steps = [ self.TStep_MAX - i*self.TStep_increase for i in range(self.NSimulations)]
    #     self.time_steps = [ ts for ts in time_steps if ts>0 ]
    #     self.NSimulations = len(self.time_steps)
    #     self.TSteps_min_max=[self.TStep_MAX,self.TStep_MAX-(self.NSimulations-1)*self.TStep_increase]
    # 
    #     #Set simulations time settings (field time + lcm(time_steps) + hardcoded duration to analyse)
    #     ts_lcm = float(np.lcm.reduce(self.time_steps))/1000. # in fs
    #     if self.ref_time/ts_lcm<self.Tpoints_min:
    #         self.yf.msg("[ERR] less than %d time points for polarization."%self.Tpoints_min)
    #         self.yf.msg("Exiting...")
    #         exit()
    #     self.yin['Field1_Tstart'] = [ts_lcm, 'fs']
    #     NETime = ts_lcm + FieldTime + self.ref_time
    #     self.yin['NETime'] = [ NETime, 'fs' ]
    #     self.NETime = NETime
    #     self.yf.msg("Total duration of simulations set to: %f fs"%NETime)
    #     self.yin['IOCachetime'] = [[ts_lcm,ts_lcm],'fs']
    # 
    # def COMPUTE_dipoles(self,DIP_folder='dipoles'):
    #     """
    #     Compute the dipoles once and for all:
    #     In order for the dipoles to be compatible with a negf run 
    #     [a default optics run does not produce compatible dipoles], 
    #     the 'negf' argument is appended which causes the calculation to crash
    #     *after* the dipoles are computed.
    #     """
    #     if not os.path.isfile('%s/%s/ndb.dipoles'%(self.RUN_path,DIP_folder)):
    #         ydipoles = YamboIn()
    #         ydipoles.arguments.append('dipoles')
    #         ydipoles.arguments.append('negf')
    #         ydipoles['DIP_ROLEs'] = self.yin['DIP_ROLEs']
    #         ydipoles['DIP_CPU'] = self.yin['DIP_CPU']
    #         ydipoles['DipBands'] = self.yin['DipBands']
    #         ydipoles.write('%s/dipoles.in'%self.RUN_path)
    #         self.yf.msg("Running dipoles...")
    #         shell = self.scheduler()
    #         shell.add_command('cd %s'%self.RUN_path)
    #         #THIS must be replaced by a more advanced submission method
    #         shell.add_command('%s -F dipoles.in -J %s -C %s 2> %s.log'%(self.yambo_rt,DIP_folder,DIP_folder,DIP_folder))
    #         shell.run()
    #         shell.clean() 
    #     else:
    #         self.yf.msg("Dipoles found.")
    # 
    #     self.DIP_folder = DIP_folder
    # 
    # def input_to_run(self,param,value,units):
    #     """
    #     Generate input for a specific run
    #     """
    #     from copy import deepcopy
    #     yrun = deepcopy(self.yin)
    #     yrun[param] = [ value, units]
    #     return yrun
    # 
    # def RUN_convergence(self,param='RTstep',units='as'):
    #     """
    #     Run the yambo_rt calculations flow.
    #     """        
    #     self.yf.msg("Running RT time step convergence...")
    #     RToutput  =    []
    #     NaN_check =    []
    #     eh_check  =    []
    #     pol_sq_check = []
    #     pol_x_check  = []
    #     time_steps = self.time_steps
    #     for i,ts in enumerate(time_steps):
    #         self.yf.msg("Running simulation for time step: %d as"%ts)
    # 
    #         # Part 1: file preparation and run
    #         filename = '%s_%05d%s.in'%(param,ts,units)
    #         folder   = filename.split('.')[0]
    #         #self.yf.msg('%s %s'%(filename,folder))
    #         yrun = self.input_to_run(param,ts,units)
    #         yrun.write('%s/%s'%(self.RUN_path,filename))
    #         shell = self.scheduler()
    #         shell.add_command('cd %s'%self.RUN_path)
    #         #THIS must be replaced by a more advanced submission method
    #         shell.add_command('%s -F %s -J %s,%s -C %s 2> %s.log'%(self.yambo_rt,filename,folder,self.DIP_folder,folder,folder))
    #         shell.run()
    #         shell.clean()
    # 
    #         # Part 2: perform single-run analysis and store output
    #         out_dir = '%s/%s'%(self.RUN_path,folder)
    #         #Read output
    #         RTDB = YamboRTDB(calc=out_dir) #Read output
    #         RToutput_no_nan, NaN_test = self.nan_test(RTDB)              #[TEST1] NaN and overflow
    #         RToutput.append(RToutput_no_nan)
    #         if NaN_test: eh_test = self.electron_conservation_test(RTDB) #[TEST2] Electron number
    #         else:        eh_test = False
    #         NaN_check.append(NaN_test) 
    #         eh_check.append(eh_test) 
    # 
    #         # Part 3: perform polarization tests between subsequent runs
    #         if i==0: passed_counter = 0
    #         if i>0: 
    #             pol_sq_test, pol_x_test, passed_counter = self.ANALYSE_pol(RToutput,eh_check,passed_counter) #[TEST3],[TEST4] Polarization squared and along field direction
    #             pol_sq_check.append(pol_sq_test)
    #             pol_x_check.append(pol_x_test)
    # 
    #         # Part 4: decide if convergence was reached or we have to keep going
    #         if passed_counter==2:
    #             TStep_passed = self.time_steps[i-2]
    #             break
    # 
    #     if passed_counter==2: self.TStep_passed = TStep_passed
    #     if passed_counter==1: self.TStep_passed = self.time_steps[-2]
    #     if passed_counter==0: self.TStep_passed = None
    # 
    #     self.NSimulations = len(RToutput)
    #     self.RToutput = RToutput
    #     self.ANALYSE_output(NaN_check,eh_check,pol_sq_check,pol_x_check,passed_counter)
    # 
    # def ANALYSE_output(self,NaN,eh,pol2,polx,passed):
    #     """
    #     Output information and suggestion for an optimal time step.
    #     - There are two values of tolerance: one for carriers, one for polarization
    #     - Four increasingly stringent checks are performed:
    #         [1] NaN and overflow check to exclude botched runs
    #         [2] Conservation of electron number check
    #         [3] Error check of |pol|^2 (assuming lowest time step as reference)
    #         [4] Error check of pol along the field direction
    #     """
    #     self.yf.msg("---------- ANALYSIS ----------")
    # 
    #     NaN_passed = sum(NaN)
    #     self.yf.msg("[1] NaN and overflow test:")
    #     self.yf.msg("    Passed by %d out of %d."%(NaN_passed,self.NSimulations))
    # 
    #     eh_passed = sum(eh)
    #     self.yf.msg("[2] Conservation of electron number test (tol=%.0e):"%self.tol_eh)
    #     self.yf.msg("    Passed by %d out of %d."%(eh_passed,self.NSimulations))
    # 
    #     pol2_passed = sum(pol2)
    #     self.yf.msg("[3] Error in |pol|^2 test (tol=%.0e):"%self.tol_pol)
    #     self.yf.msg("    Passed by %d out of %d."%(pol2_passed,self.NSimulations-1))
    # 
    #     polx_passed = sum(polx)
    #     self.yf.msg("[4] Error in pol along field test (tol=%.0e):"%self.tol_pol)
    #     self.yf.msg("    Passed by %d out of %d."%(polx_passed,self.NSimulations-1))
    # 
    #     if passed == 1:
    #         self.yf.msg(" ")
    #         self.yf.msg("[WARNING] The lowest time step passed all the tests, but the")
    #         self.yf.msg("          additional safety run with a reduced step was")
    #         self.yf.msg("          not done due to NSimulations limit being reached.")
    # 
    #     if self.NSimulations == 2 or self.NSimulations == 3:
    #         self.yf.msg(" ")
    #         self.yf.msg("[WARNING] The largest time step already looks converged.")
    # 
    #     tp = self.TStep_passed
    #     self.yf.msg(" ")
    #     self.yf.msg("Based on the analysis, the suggested time step is: ")
    #     if tp is not None: self.yf.msg("### %d as ###"%tp)
    #     else: self.yf.msg("[ERR] NSimulations limit reached before converged value was found.")
    #     self.yf.msg("------------------------------")        
    # 
    # def ANALYSE_pol(self,RToutput,eh_check,passed):
    #     """
    #     Driver with the logical structure to manage polarization tests
    #     """
    #     if eh_check[-1]==True and eh_check[-2]==True:
    #         pol_sq_test = self.pol_error_test(RToutput,which_pol='pol_sq')
    #         pol_x_test  = self.pol_error_test(RToutput,which_pol='pol_along_field')
    # 
    #         if pol_sq_test and pol_x_test: passed = passed + 1
    #         else: passed = 0
    # 
    #     else: 
    #         pol_sq_test = False
    #         pol_x_test  = False
    # 
    #     return pol_sq_test, pol_x_test, passed
    # 
    # def nan_test(self,RTDB):
    #     """ 
    #     Check computed polarizations for NaN values.
    #     """
    #     NaN_test = True
    #     # Check for NaN
    #     if np.isnan(RTDB.polarization).any() or np.isnan(RTDB.diff_carriers).any(): 
    #         RTDB.polarization = np.nan_to_num(RTDB.polarization) #Set to zero for plots
    #         NaN_test = False 
    #         #self.yf.msg("[WARNING] Yambo produced NaN values during this run")
    #     # Check for +/-Infinity
    #     if np.greater(np.abs(RTDB.polarization),overflow).any():
    #         RTDB.polarization[np.abs(RTDB.polarization)>overflow] = 0. #Set to zero for plots
    #         NaN_test = False
    #         #self.yf.msg("[WARNING] Yambo produced Infinity values during this run")           
    # 
    #     return RTDB, NaN_test
    # 
    # def electron_conservation_test(self,RTDB):
    #     """
    #     Tests if elements of ratio_carriers are greater than tolerance.
    #     If any of them is, then the simulation in question has not passed the eh_test.
    #     """
    #     eh_carriers = np.greater(RTDB.ratio_carriers,self.tol_eh)
    #     if any(eh_carriers): eh_test = False
    #     else:                eh_test = True
    #     return eh_test
    # 
    # def pol_error_test(self,RTout,which_pol):
    #     """
    #     Computes the relative errors of the polarizations for each cached time.
    #     The cached times coincide for different runs.
    #     """
    #     pol_analyse= []
    #     pol_n1 = RTout[-1].polarization   
    #     pol_n0 = RTout[-2].polarization 
    #     if which_pol == 'pol_sq':  # Test for |pol|^2
    #         pol_analyse_n1 = pol_n1[0]*pol_n1[0] + pol_n1[1]*pol_n1[1] + pol_n1[2]*pol_n1[2] 
    #         pol_analyse_n0 = pol_n0[0]*pol_n0[0] + pol_n0[1]*pol_n0[1] + pol_n0[2]*pol_n0[2] 
    #     if which_pol == 'pol_along_field': # Test for pol along field
    #         dr, _ = self.pol_along_field()
    #         pol_analyse_n1 = pol_n1[dr]
    #         pol_analyse_n0 = pol_n0[dr]            
    # 
    #     #Perform the test
    #     rel_err_pol = (pol_analyse_n1-pol_analyse_n0)/self.tol_pol
    #     error = np.greater(rel_err_pol,1.).any()
    #     if error: pol_test = False
    #     else:     pol_test = True      
    # 
    #     return pol_test
    # 
    # def pol_along_field(self):
    #     field = self.yin['Field1_Dir']
    #     if field[0]!=0.:   dr,pol_label=[0,'pol-x']
    #     elif field[1]!=0.: dr,pol_label=[0,'pol-y']
    #     elif field[2]!=0.: dr,pol_label=[0,'pol-z']
    #     else:              dr,pol_label=[0,'pol-x']
    #     return dr,pol_label
    # 
    # def PLOT_output(self,save_dir='plots'):
    #     """
    #     Generic plots generated by default, to be accessed by the user
    #     """
    #     import matplotlib.pyplot as plt
    # 
    #     self.yf.msg("Plotting results.")
    #     out_dir = '%s/%s'%(self.RUN_path,save_dir)
    #     if not os.path.isdir(out_dir): 
    #         shell = self.scheduler()
    #         shell.add_command('mkdir -p %s'%out_dir)
    #         shell.run()
    #         shell.clean()
    # 
    #     time_steps = self.time_steps
    #     lwidth=0.8
    #     ts_colors = plt.cm.gist_rainbow(np.linspace(0.,1.,num=self.NSimulations))
    # 
    #     # Plot for each time step
    #     for ts in range(self.NSimulations):
    # 
    #         pol   = self.RToutput[ts].polarization
    #         pol_sq = pol[0]*pol[0] + pol[1]*pol[1] + pol[2]*pol[2]
    #         times = np.linspace(0.,self.NETime,num=pol.shape[1])
    #         f, (axes) = plt.subplots(4,1,sharex=True)
    #         axes[0].plot(times, pol[0], '-', lw=lwidth, color='blue',  label='pol-x')
    #         axes[1].plot(times, pol[1], '-', lw=lwidth, color='green', label='pol-y') 
    #         axes[2].plot(times, pol[2], '-', lw=lwidth, color='red',   label='pol-z')
    #         axes[3].plot(times, pol_sq, '-', lw=lwidth, color='orange',label='|pol|^2')
    #         for ax in axes:
    #             ax.axhline(0.,lw=0.5,color='gray',zorder=-5)
    #             ax.legend(loc='upper left')
    #         f.tight_layout()
    # 
    #         plt.savefig('%s/polarizations_%das.png'%(out_dir,self.time_steps[ts]),format='png',dpi=150)
    # 
    #     # Plot for all time steps
    #     f, (axes) = plt.subplots(4,1,sharex=True)
    #     for ts in range(self.NSimulations):
    # 
    #         label = '%das'%time_steps[ts]
    #         pol   = self.RToutput[ts].polarization
    #         pol_sq = pol[0]*pol[0] + pol[1]*pol[1] + pol[2]*pol[2]
    #         times = np.linspace(0.,self.NETime,num=pol.shape[1])
    #         axes[0].plot(times, pol[0], '-', lw=lwidth, color=ts_colors[ts], label=label)
    #         axes[1].plot(times, pol[1], '-', lw=lwidth, color=ts_colors[ts], label=label)
    #         axes[2].plot(times, pol[2], '-', lw=lwidth, color=ts_colors[ts], label=label)
    #         axes[3].plot(times, pol_sq, '-', lw=lwidth, color=ts_colors[ts], label=label)
    #         handles, labels = axes[3].get_legend_handles_labels()
    #     for ax in axes: ax.axhline(0.,lw=0.5,color='gray',zorder=-5)
    # 
    #     f.legend(handles, labels, loc='center right')
    #     f.tight_layout()
    # 
    #     plt.savefig('%s/polarizations_comparison.png'%out_dir,format='png',dpi=150) 
    # 
    #     # Plot for all time steps |pol|^2
    #     f, (axes) = plt.subplots(self.NSimulations,1,sharex=True)
    #     for ts in range(self.NSimulations):
    # 
    #         pol = self.RToutput[ts].polarization
    #         pol_sq = pol[0]*pol[0] + pol[1]*pol[1] + pol[2]*pol[2]
    #         times = np.linspace(0.,self.NETime,num=pol.shape[1])
    #         pol_ts_label = "%das"%time_steps[ts]
    #         axes[ts].plot(times, pol_sq, '-', lw=lwidth, color=ts_colors[ts], label=pol_ts_label)
    #     for ax in axes:
    #         ax.axhline(0.,lw=0.5,color='gray',zorder=-5)
    #         ax.legend(loc='upper left')
    #     f.tight_layout()
    # 
    #     plt.savefig('%s/polarizations_squared.png'%out_dir,format='png',dpi=150)
    # 
    #     # Plot for all time steps along field direction
    #     dr, pol_label = self.pol_along_field()
    #     f, (axes) = plt.subplots(self.NSimulations,1,sharex=True)
    #     for ts in range(self.NSimulations):
    # 
    #         pol = self.RToutput[ts].polarization
    #         times = np.linspace(0.,self.NETime,num=pol.shape[1])
    #         pol_ts_label = "%s_%das"%(pol_label,time_steps[ts])
    #         axes[ts].plot(times, pol[dr], '-', lw=lwidth, color=ts_colors[ts], label=pol_ts_label) 
    #     for ax in axes:
    #         ax.axhline(0.,lw=0.5,color='gray',zorder=-5)
    #         ax.legend(loc='upper left')
    #     f.tight_layout()        
    # 
    #     plt.savefig('%s/polarizations_field_direction.png'%out_dir,format='png',dpi=150)
    # 
