from yambopy import *
from qepy import *
from schedulerpy import *
import os
import subprocess
from glob import glob
from copy import deepcopy
import matplotlib.pyplot as plt

class YamboDG_Optimize():
    """ 
    Class to generate and run convergence tests for the yambo double grid. 
    
    ** This class is useful but complex, read the description well AND/OR check its tutorial! **

    - Needs a quantum espresso scf save folder
    - Needs nscf (qe) and yambo [desired calculation type] inputs PASSED AS PwIn AND YamboIn OBJECTS
         -- If converging the double grids, an independent-particles (ip) yambo input is required 
    - Needs lists of coarse grids (CG) and fine grids (FG); [NOTE] Only random FG presently implemented. 
    - Additional arguments: directory paths, file names, experimental laser energy [for absorption], etc.
    - The workflow is divided into FOUR STEPS that can be executed separately:
        1. nscf CG [STEPS='1']
        2. nscf FG and yambo CG [STEPS='2']
        3. yambo FG [STEPS='3']
        4. if 'converge_DG' is on (therefore with yambo--> ip):
             -- TODO: Analyis, report, plot results and give ip-converged value [STEPS='4']
             -- TODO: Move the double grid generation functions to a different submodule       
 
    - Scheme of the workflow:
        -- If job submissions are used, the workflow is better submitted in subsequent steps
        -- If planning a parallel traversal (each independent branch simultaneously) of this tree 
           with job submissions, see function branch_wise_flow
           
            NSCF                                       YAMBO
            |                                          |
    step 1  CG_1              CG_2 ... CG_N            |
            |                 |                        |
    step 2  FG_11 ... FG_M1   FG_12 ... FG_M2 ...      CG_1              CG_2 ... CG_N
                                                       |                 |
    step 3                                             FG_11 ... FG_M1   FG_12 ... FG_M2 ...
                                                        \         |      |         /
                                                         \        \      |        /
    step 4                                                 _________ PLOTS ______
    
    Some optional variables
      
      - E_laser: external laser energy (for RT checks)
      - STEPS: which workflow steps to execute
      - RUN: if False, no job is submitted
      - converge_DG: if True, enables automatic double grid convergence; requires IP yambo input.
      - nscf_out, y_out_dir: pw and yambo output directories
      - qe_scheduler,y_scheduler: SchedulerPy object for cluster submission (default: bash)
      - wait_all: if cluster submission is on, forces the master python process to wait for all sumbitted jobs to complete before exiting
      - yambo_calc_type: name the yambo calculations
      - yambo_exec: either yambo, yambo_ph or yambo_rt
      - save_type: simple, elph, expanded_elph, fixsymm, fixsymm+elph
            -- NOTE: if *elph is used: prepare a symlink to elph_dir in RUN_path
      - noyambo: if RUN is True, run only QE steps plus SAVEs, but not actual yambo calculations

    Example of use (frontend):           
        .. code-block:: python
    
            YamboDG_Optimize(cg_grids,fg_grids,prefix,scf_path,pseudo_path,...,STEPS='all')
            
    Example of use (job submission: each step dependent on the one before)
    
        .. code-block:: python 
            
            scheduler1 = Scheduler(...)
            scheduler2 = Scheduler(...,dependency=scheduler1)
            scheduler3 = Scheduler(...,dependency=scheduler2)
            
            scheduler1.add_command('python dg_script.py -steps 1')
            sheduler1.run()
            
            scheduler2.add_command('python dg_script.py -steps 2')
            sheduler2.run()
            
            scheduler3.add_command('python dg_script.py -steps 3')
            sheduler3.run()            
       
        ..
        
        .. code-block:: dg_script.py 
           
           YamboDG_Optimize(cg_grids,fg_grids,prefix,scf_path,pseudo_path,...,wait_all==True,STEPS=args.steps)
    
    TO DO:
      - Separate double grid generation and double grid convergence (simple option 'converge_DG' might suffice)
      - If automatic DG convergence assessment is on, then implement MOMENTA of the abs spectra as a method to check convergence
      - Change check_nscf_completed into common/check_qe_completed

    """
    
    def __init__(self,cg_grids,fg_grids,prefix,nscf_input,ya_input,E_laser=0.,STEPS='all',RUN=True, converge_DG=False,\
                 scf_save_path='./scf',pseudo_path='./pseudos',RUN_path='./',nscf_out="nscf",y_out_dir="results",\
                 qe_scheduler=None,y_scheduler=None,wait_all=True,pw_exec_path='',yambo_exec_path='',\
                 yambo_exec='yambo',save_type='simple',yambo_calc_type="yambo",noyambo=False):

        #Configuring schedulers
        self.frontend = Scheduler.factory(scheduler="bash")
        if y_scheduler is not None and qe_scheduler is not None: #Here we use, e.g., slurm 
            self.qejobrun, self.yjobrun = qe_scheduler, y_scheduler
            self.wait_up = True #This will enable calls to a function that will make the code wait upon completion of previous submitted job
            self.job_folders = []
            self.job_shells  = []
        elif y_scheduler is None and qe_scheduler is None: # Both schedulers must be present to activate job submission
            self.qejobrun, self.yjobrun = Scheduler.factory(scheduler="bash"), Scheduler.factory(scheduler="bash")
            self.wait_up = False
        else: raise UserWarning('Submission mode is on only for either yambo or qe')

        #Setting global variables
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
        self.yambo_calc_type = yambo_calc_type
        self.noyambo = noyambo
        self.E_laser = E_laser
        
        # Initialize JOBID lists (used only in submission mode)
        self.qe_id_cg = [ -1 for cg in self.cg_strings ]
        self.ya_id_cg = [ -1 for cg in self.cg_strings ]
        self.qe_id_fg = [ [ -1 for fg_i in fg ] for fg in self.fg_strings ]
        self.ya_id_fg = [ [ -1 for fg_i in fg ] for fg in self.fg_strings ]
        
        # Path of nscf and ip calculations and final plots
        if converge_DG: self.yambo_calc_type='ip'
        if not os.path.isdir(RUN_path): os.mkdir(RUN_path)
        self.nscf_dir = '%s/nscf_grids'%RUN_path
        self.yambo_dir = '%s/%s_grids'%(RUN_path,self.yambo_calc_type)
        self.plot_dir = '%s/plots'%RUN_path
        if not os.path.isdir(self.nscf_dir): os.mkdir(self.nscf_dir)
        if not os.path.isdir(self.yambo_dir): os.mkdir(self.yambo_dir)

        #Inputs
        self.nscf_inp = nscf_input
        self.ya_inp   = ya_input
        if converge_DG:
            yambo_exec = 'yambo'
            self.ip_input_is_there()
        
        #Executables
        if yambo_exec_path != '': yambo_exec_path_aux=yambo_exec_path+'/'
        if pw_exec_path != '': pw_exec_path+='/'
        self.pw = pw_exec_path + 'pw.x'
        self.ypp = yambo_exec_path_aux + 'ypp'
        self.yambo = yambo_exec_path_aux + yambo_exec
        
        # Automatically determine which SAVE to create (better to specify it explicitly)
        if yambo_exec == 'yambo': save_type='simple'
        elif yambo_exec == 'yambo_ph' and not save_type[-4:]=='elph': save_type='elph'
        elif yambo_exec == 'yambo_rt' and not save_type[3:]=='fix': save_type='fixsymm'
        # Deal with elph_path
        elph_path = None
        if save_type[-4:]=='elph': 
            if ( not os.path.isdir('%s/elph_dir'%self.RUN_path) ) and ( not os.path.isfile('%s/elph_dir'%self.RUN_path) ):
                raise FileNotFoundError('Please mv or symlink the elph_dir folder to the RUN_path of this workflow.')
            else:
                elph_path=self.RUN_path

        #Start IO
        self.yf = YamboIO(out_name='YAMBOPY_double-grid_Optimize.log',out_path=self.RUN_path,print_to_shell=True)
        self.yf.IO_start()
        
        if converge_DG: self.yf.msg('#### DOUBLE GRID CONVERGENCE WORKFLOW ####')
        else: self.yf.msg('#### DOUBLE GRID WORKFLOW FOR %s CALCULATIONS ####'%self.yambo_calc_type)
        
        self.driver(STEPS,RUN,nscf_out,y_out_dir,save_type,elph_path,yambo_exec_path,converge_DG)
        
        #End IO        
        self.yf.IO_close()
        
        # Check if python must exit immediately or after all submitted jobs have completed
        if self.wait_up:
            if wait_all: wait_for_all_jobs(self.job_shells,self.job_folders)
            for shell in self.job_shells: shell.clean()
    
    def driver(self,STEPS,RUN,nscf_out,y_out_dir,save_type,elph_path,yambo_exec_path,converge_DG):
        """
        Worflow driver.
        
        It runs the following:
            - setup functions
            - job submission functions
            - double grid convergence functions
        """
        if STEPS=='1' or STEPS=='all': 
            self.yf.msg("------  STEP 1 ------")
            self.setup_cg()
            if RUN: self.run_jobs(nscf_out,y_out_dir)
        
        if STEPS=='2' or STEPS=='all':
            self.yf.msg("------  STEP 2 ------")
            for ig,cg in enumerate(self.cg_strings):
                # NSCF status check
                calc_dir  = '%s/%s_coarse_grid'%(self.nscf_dir,cg)            
                calc = self.check_nscf_completed(calc_dir,nscf_out)
                if calc: self.yf.msg("     NSCF CG %s found."%cg)
                else: self.yf.msg("     NSCF CG %s NOT found."%cg)
                # YAMBO status check
                ycalc_dir = '%s/%s_coarse_grid'%(self.yambo_dir,cg)
                ycalc = self.yambo_output_is_NOT_there(ycalc_dir,y_out_dir)
                if ycalc: self.yf.msg("     YAMBO CG %s NOT found."%cg)
                else: self.yf.msg("     YAMBO CG %s found."%cg)      
                #          
                if calc and ycalc:
                    yambo_dir = '%s/%s_coarse_grid'%(self.yambo_dir,cg)
                    CreateYamboSave(self.prefix,save_type=save_type,nscf=calc_dir,elph_path=elph_path,\
                                    database=os.path.abspath(yambo_dir),yambo_exec_path=yambo_exec_path,printIO=False)
                    self.setup_fg(calc_dir,yambo_dir,self.fg_grids[ig],self.fg_strings[ig])
            if RUN: self.run_jobs(nscf_out,y_out_dir)
                    
        if STEPS=='3' or STEPS=='all':
            self.yf.msg("------  STEP 3 ------")
            for ig,cg in enumerate(self.cg_strings):
                for iff,fg in enumerate(self.fg_strings[ig]):
                    # NSCF status check
                    calc_dir  = '%s/%s_coarse_grid/%s'%(self.nscf_dir,cg,fg)
                    calc = self.check_nscf_completed(calc_dir,nscf_out)
                    if calc: self.yf.msg("     NSCF CG %s FG %s found."%(cg,fg))
                    else: self.yf.msg("     NSCF CG %s FG %s NOT found."%(cg,fg))
                    # YAMBO status check
                    ycalc_dir = '%s/%s_coarse_grid/%s'%(self.yambo_dir,cg,fg)
                    ycalc = self.yambo_output_is_NOT_there(ycalc_dir,y_out_dir)
                    if ycalc: self.yf.msg("     YAMBO CG %s FG %s NOT found."%(cg,fg))
                    else: self.yf.msg("     YAMBO CG %s FG %s found."%(cg,fg))   
                    #      
                    if calc and ycalc:
                        yambo_dir = '%s/%s_coarse_grid/%s'%(self.yambo_dir,cg,fg)
                        if not os.path.isfile('%s/SAVE/ndb.Double_Grid'%yambo_dir):
                            CreateYamboSave(self.prefix,save_type='simple',nscf=calc_dir,elph_path=None,\
                                            database="%s/dg_SAVE"%os.path.abspath(yambo_dir),yambo_exec_path=yambo_exec_path,printIO=False)
                            self.setup_yambo_fg(yambo_dir,self.fg_grids[ig][iff],y_out_dir)
            if RUN: self.run_jobs(nscf_out,y_out_dir)
        
        # This is a plotting routine in the ip case. It has to be updated to a full convergence analysis and report.
        # The involved functions can possibly be moved into another file
        # 
        if (STEPS=='4' or STEPS=='all') and converge_DG :
            self.yf.msg("----------  STEP 4 ----------")
            self.yf.msg("-- Double grid convergence --")
            yout = 'o-%s.eps_q1_ip'%y_out_dir
            ip_data = []
            ip_labels = []
            #
            for ig,cg in enumerate(self.cg_strings):
                ip_data_temp = []
                ip_labels_temp = []
                temp_path = '%s/%s_coarse_grid/%s/%s'%(self.yambo_dir,cg,y_out_dir,yout)
                if os.path.isfile(temp_path):
                    self.yf.msg("     IP CG %s found."%cg)
                    ip_data_temp.append(np.loadtxt(temp_path))
                    ip_labels_temp.append(cg)
                else:
                    self.yf.msg("     IP CG %s NOT found."%cg)
                for iff,fg in enumerate(self.fg_strings[ig]):
                    temp_path = '%s/%s_coarse_grid/%s/%s/%s'%(self.yambo_dir,cg,fg,y_out_dir,yout)
                    if os.path.isfile(temp_path):
                        self.yf.msg("     IP CG %s FG %s found."%(cg,fg))
                        ip_data_temp.append(np.loadtxt(temp_path))
                        ip_labels_temp.append('%s N_FG=%d'%(cg,self.fg_grids[ig][iff]))
                    else:
                        self.yf.msg("     IP CG %s FG %s NOT found."%(cg,fg))
                ip_data.append(ip_data_temp)
                ip_labels.append(ip_labels_temp)
            #
            ip_data = [ip for ip in ip_data if ip != []]
            ip_labels = [ip for ip in ip_labels if ip != []]
            self.plot_ip_spectra(ip_data,ip_labels,self.E_laser)
        
    def setup_cg(self):
        """ First step of the workflow: setup CG folder tree and CG nscf calculations
        """
        for ig,cg in enumerate(self.cg_grids):
            work_dir = "%s/%s_coarse_grid"%(self.nscf_dir,self.cg_strings[ig])
            if not os.path.isdir(work_dir):
                os.mkdir(work_dir)
                self.generate_inputfile(work_dir,cg)
                os.system('cp -r %s/%s.save %s'%(self.scf_save_path,self.prefix,work_dir))
            yambo_dir = "%s/%s_coarse_grid"%(self.yambo_dir,self.cg_strings[ig])
            if not os.path.isdir(yambo_dir):
                os.mkdir(yambo_dir)
                self.ya_inp.write('%s/%s.in'%(yambo_dir,self.yambo_calc_type))
                
    def setup_fg(self,nscf_cg_dir,yambo_cg_dir,fg_grid,fg_string):
        """ Second step of the workflow: setup FG folder tree and FG (CG) nscf (yambo) calculations
        """
        for iff,fg in enumerate(fg_grid):
            fg_nscf_inp = '%s_fg.nscf'%self.prefix
            rand_nm = 'random_kpt_list_%s.dat'%fg_string[iff]
            ypp_inp = 'ypp_fg.in'
            work_dir = "%s/%s"%(nscf_cg_dir,fg_string[iff])
            if not os.path.isdir(work_dir):
                os.mkdir(work_dir)
                os.system('cp -r %s/%s.save %s'%(self.scf_save_path,self.prefix,work_dir))
            yambo_dir = "%s/%s"%(yambo_cg_dir,fg_string[iff])
            if not os.path.isdir(yambo_dir):
                os.mkdir(yambo_dir)
                os.mkdir("%s/dg_SAVE"%yambo_dir)
                self.ya_inp.write('%s/%s.in'%(yambo_dir,self.yambo_calc_type))
            if not os.path.isfile('%s/%s'%(work_dir,fg_nscf_inp)):
                self.generate_ypp_input_random_grid(yambo_cg_dir,fg,ypp_inp)
                ypp_run = self.frontend
                ypp_run.add_command('cd %s; %s -F %s > ypp_fg.log'%(yambo_cg_dir,self.ypp,ypp_inp))
                ypp_run.add_command('mv o.random_k_pts %s'%rand_nm)
                ypp_run.add_command('cp %s %s'%(rand_nm,work_dir))
                ypp_run.run()
                kpts_rndm = np.loadtxt('%s/%s'%(work_dir,rand_nm))
                if len(kpts_rndm) != fg:
                    self.yf.msg("[WARNING] Actual fine grid number of kpts is %d instead of %d"%(len(kpts_rndm),fg))
                self.generate_inputfile(work_dir,fg,klist=kpts_rndm)

    def setup_yambo_fg(self,yambo_fg_dir,fg_num,yresults_dir,clean_dg_saves=True):
        """ Third step of the workflow: map FG to CG and FG yambo calculations
        """
        ypp_inp = 'ypp_map.in'
        os_run = self.frontend
        if os.path.isfile('%s/../%s/ndb.dipoles'%(yambo_fg_dir,yresults_dir)):
            os_run.add_command('cd %s; cp ../%s/ndb.dipoles* ../SAVE/ ; cp -r ../SAVE .'%(yambo_fg_dir,yresults_dir))
        else:
            os_run.add_command('cd %s; cp -r ../SAVE .'%yambo_fg_dir)
        os_run.run()
        self.generate_ypp_input_map_grid(yambo_fg_dir,fg_num,ypp_inp)
        ypp_run = self.frontend
        ypp_run.add_command('cd %s; %s -F %s > ypp_map.log'%(yambo_fg_dir,self.ypp,ypp_inp))
        ypp_run.run()
        if os.path.isfile('%s/SAVE/ndb.Double_Grid'%yambo_fg_dir):
            self.yf.msg("  -> Double Grid ready. <-")
            if clean_dg_saves: 
                os.system('rm -r %s/dg_SAVE/SAVE'%yambo_fg_dir)
                os.system('touch %s/dg_SAVE/DOUBLE_GRID_SAVE_REMOVED_TO_SAVE_DISK_SPACE'%yambo_fg_dir)
        else:
            self.yf.msg("  -> Double Grid NOT ready. <-")
    
    def generate_inputfile(self,folder,kpts,klist=None):
        """ Modify nscf input file in case of CG or FG kpoint grids
        """
        import copy
        qe_input = copy.deepcopy(self.nscf_inp)
        qe_input.control['pseudo_dir'] = "'%s'"%self.pseudo_path
        
        if klist is None:
            qe_input.kpoints = kpts
        else:
            qe_input.ktype = "crystal"
            qe_input.klist = klist
            
        qe_input.write('%s/%s.nscf'%(folder,self.prefix))
    
    def generate_ypp_input_random_grid(self,folder,fg_num,inp_nm):
        """ Create ypp input file for the generation of the FG coordinates [units: yambo "rlu" -> qe "crystal"]
        """
        yppin = YamboIn()
        yppin.arguments.append('bzgrids')
        yppin.arguments.append('Random_Grid')
        yppin['OutputAlat'] = self.nscf_inp.system['celldm(1)']
        yppin['cooOut'] = "rlu"
        yppin['BZ_random_Nk'] = fg_num

        yppin.write('%s/%s'%(folder,inp_nm))
        
    def generate_ypp_input_map_grid(self,folder,fg_num,inp_nm):
        """ Create ypp input file for the mapping of the FG grid to the CG grid
        """
        yppin = YamboIn()
        yppin.arguments.append('kpts_map')
        yppin['FineGd_mode']='unexpanded'
        yppin['BZ_DbGd_Nk']=fg_num
        yppin.arguments.append('SkipCheck')
        yppin['FineGd_DB1_paths'] = ['./dg_SAVE']

        yppin.write('%s/%s'%(folder,inp_nm))
            
    def run_jobs(self,out_qe,out_yambo):
        """ Workflow
        """
        #MODULARIZATION ISSUE: 
        #       remember that function check_nscf_completed has a dependency 
        #       on the name of the qe output file - it has to be '*.out' - and 
        #       hence it depends on how this name is given in function shell_run     
                
        for ig,cg in enumerate(self.cg_strings): # ---------- Outer COARSE GRID loop ----------
            # temp_dir: where qe calculations are run
            # save_dir: where yambo calculations are run   

            temp_dir = '%s/%s_coarse_grid'%(self.nscf_dir,cg)
            
            if not self.check_nscf_completed(temp_dir,out_qe):
                self.yf.msg("Running NSCF CG %s..."%cg)      ############## Run NSCF CG
                self.qe_id_cg[ig] = self.shell_run("qe_%s"%cg,temp_dir,out_qe,'qe')
                if self.wait_up: self.job_folders.append(temp_dir)

            save_dir = '%s/%s_coarse_grid'%(self.yambo_dir,cg)
            if os.path.isdir('%s/SAVE'%save_dir):

                if self.yambo_output_is_NOT_there(save_dir,out_yambo) and not self.noyambo:
                    self.yf.msg("Running YAMBO CG %s..."%cg)     ############ Run YAMBO CG
                    self.ya_id_cg[ig] = self.shell_run("ya_%s"%cg,save_dir,out_yambo,'y')    # depends on JOBID='%d'%self.qe_id_cg[ig])
                    if self.wait_up: self.job_folders.append(save_dir)

            try: self.fg_strings[ig]
            except IndexError as err: raise Exception('No fine grid(s) provided for CG %s.'%cg) from err 
                      
            for iff,fg in enumerate(self.fg_strings[ig]):      # ---------- Inner FINE GRID loop ----------

                temp_dir = '%s/%s_coarse_grid/%s'%(self.nscf_dir,cg,fg)
                if os.path.isdir(temp_dir):

                    if not self.check_nscf_completed(temp_dir,out_qe):
                        self.yf.msg("Running NSCF CG %s FG %s..."%(cg,fg)) ######### Run NSCF FG
                        self.qe_id_fg[ig][iff] = self.shell_run("qe_%s"%cg,temp_dir,out_qe,'qe') # depends on JOBID='%d'%self.qe_id_cg[ig])
                        if self.wait_up: self.job_folders.append(temp_dir)

                save_dir = '%s/%s_coarse_grid/%s'%(self.yambo_dir,cg,fg)
                if os.path.isdir('%s/SAVE'%save_dir) and os.path.isfile('%s/SAVE/ndb.Double_Grid'%save_dir):
                    
                    if self.yambo_output_is_NOT_there(save_dir,out_yambo) and not self.noyambo:
                        self.yf.msg("Running YAMBO CG %s FG %s..."%(cg,fg))     ############ Run YAMBO FG
                        self.ya_id_fg[ig][iff] = self.shell_run("ya_%s"%cg,save_dir,out_yambo,'y') # depends on JOBID='%d:%d'%(self.ya_id_cg[ig],self.qe_id_fg[ig][iff]))
                        if self.wait_up: self.job_folders.append(save_dir)

    def shell_run(self,jname,run_dir,out_dir,exec,JOBID=None):
        """ 
        Submit job
        
            exec:
                qe: runs pw.x
                y:  runs the yambo* executable of choice
            
            jname: job name
            JOBID: job id of simulation that the present job has a dependency on
            run_dir: where job is run
            out_dir: name of output (folder and log for yambo; '%s.out'%out_dir file for qe)  
            
            returns id of present submitted job          
        """
        if exec=='qe': shell = deepcopy(self.qejobrun)
        if exec=='y':  shell = deepcopy(self.yjobrun)
        shell.name = '%s_%s'%(jname,shell.name)
        
        # Add dependency if specified
        if self.wait_up and JOBID is not None:
            dependency='afterok:%s'%JOBID
            shell.dependency=dependency
        
        if exec=='qe': shell.add_mpirun_command('%s -inp %s.nscf > %s.out'%(self.pw,self.prefix,out_dir))
        if exec=='y': shell.add_mpirun_command('%s -F %s.in -J %s -C %s 2> %s.log'%(self.yambo,self.yambo_calc_type,out_dir,out_dir,out_dir))
        shell.run(filename='%s/%s.sh'%(run_dir,exec)) ### Specify run path
        
        # Manage submissions if specified
        #if self.wait_up: wait_for_job(shell,run_dir)
        if self.wait_up: this_job_id = shell.jobid
        else:            this_job_id = -1
        if self.wait_up: self.job_shells.append(shell) 
        else: shell.clean()
        
        return this_job_id
    
    def ip_input_is_there(self):
        """ Check if yambo ip input is correctly given in converge_DG case
        """
        condition = 'chi' in self.ya_inp.arguments and \
                    'optics' in self.ya_inp.arguments and \
                    self.ya_inp['Chimod']=='IP'
        if not condition:
            raise FileNotFoundError("IP input file not recognised: are you sure you specified 'chi' and 'optics' runlevels and Chimod='IP'?")

    def check_nscf_completed(self,folder,out_prefix):
        """ Check if nscf calculation has completed and proceed
        """
        status = True
        try:
            check = subprocess.check_output("grep JOB %s/%s*"%(folder,out_prefix), shell=True, stderr=subprocess.STDOUT)
            check = check.decode('utf-8')
            check = check.strip().split()[-1]
        except subprocess.CalledProcessError as e:
            check = ""
        if check != "DONE.": status = False
        return status

    def yambo_output_is_NOT_there(self,calc_dir,results_dir):
        """ Check if yambo produced non-empty outputs
        """
        condition = 0
        out_files = glob('%s/%s/o-%s.*'%(calc_dir,results_dir,results_dir))
        for out_file in out_files:
            test = ( not os.path.isfile(out_file) ) or ( os.path.isfile(out_file) and os.stat(out_file).st_size == 0 )
            if test: condition+=1
        if condition==len(out_files): return True # This means no output file has been produced: this calculation must be run
        elif condition==0: return False # This means all output files have been produced: this calculation must be skipped
        else: raise UserWarning('Some output files have been produced and some not -- check the above calculation.') 
    
    #TODO: move convergence-related functions to a new submodule
    def plot_ip_spectra(self,data,labels,w_laser,fig_name='ip_spectra',xlims=None):
        """ Plot results in a dynamic layout (chosen by FP)
        """
        if not os.path.isdir(self.plot_dir): os.mkdir(self.plot_dir)
        lwidth = 0.8
        f, (axes) = plt.subplots(1,len(data),sharey=True,figsize=(12,4))

        for ig,cg in enumerate(data):
            colors = plt.cm.gist_rainbow(np.linspace(0.,1.,num=len(cg)))
            axes[ig].set_ylim(0.,1.1*np.max(cg[0][:,1]))
            if xlims is not None: axes[ig].set_xlim(xlims[0],xlims[1])
            for iff,fg in enumerate(cg):
                axes[ig].plot(fg[:,0],fg[:,1], '-', lw = lwidth, color = colors[iff], label = labels[ig][iff])
            axes[ig].fill_between(cg[0][:,0], cg[-1][:,1], color = colors[-1], alpha = 0.2)

        for ax in axes:
            ax.axvline(w_laser,lw=0.8,color='black',zorder=-1)
            ax.legend()

        f.tight_layout()

        plt.savefig('%s/%s.png'%(self.plot_dir,fig_name),format='png',dpi=150)
    
    def clean_everything(self):
        """ Remove workflow tree
        """
        rm_run = self.frontend
        rm_run.add_command('rm -rf %s %s %s'%(self.nscf_dir,self.yambo_dir,self.plot_dir))
        rm_run.run()
        self.yf.msg("Workflow removed.")
    
    def branch_wise_flow(self):
        """
        The workflow dependencies are complicated.
        We understand them in two steps.
        
        First, how the actual job submissions depend on each other:
        
                qe_CG_1                     qe_CG_2              ...            qe_CG_n
                   |                          |                                    |
             ______|________________          |                                    |    
             |         |           |       (same in parallel)                   (same in parallel)
             |         |           |       
        qe_FG_1  ... qe_FG_n    y_CG_1
           |           |______ ... | ... _____
           |_______________________|_________|
                              |              |
                              |              |
                            y_FG_1   ...   y_FG_n      
                            
                            if converge_DG: REDUX ALL JOBS
                                                  |
                                                  |
                                                PLOTS  
                            
        All the separate branches of the workflow could be run sequentially using the dependency system of the scheduler.
        
        Second, how the workflow functions actually depend on each other (single tree branch is shown here):
        
                        qe_CG_1           
                           | 
                     ______|________________                    first barrier   
                     |         |           | 
                     |         |           |
                     |---------|--------- SAVE   
                     |         |           |                    second barrier
                setup_fg    setup_fg       |                (very small delay here)
                     |         |           |                    third barrier
                     |         |           |
                qe_FG_1  ... qe_FG_n    y_CG_1
                   |           |______ ... | ... _____
                   |_______________________|_________|          fourth barrier
                                      |              |
                                      |              | 
                                setup_yambo_fg  setup_yambo_fg
                                      |              |          fifth barrier
                                      |              |
                                    y_FG_1   ...   y_FG_n 
                                    
        For this reason, it's better to plan a submission 'STEP by STEP' of the workflow instead than 'branch by branch'. 
        However, the latter way would be more efficient if implemented. 
        """
