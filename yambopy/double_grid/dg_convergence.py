from yambopy import *
from qepy import *
from schedulerpy import *
import os
import subprocess
import matplotlib.pyplot as plt

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
        
    - Scheme of the workflow:
    
            NSCF                                       IP
            |                                          |
    step 1  CG_1              CG_2 ... CG_N            |
            |                 |                        |
    step 2  FG_11 ... FG_M1   FG_12 ... FG_M2 ...      CG_1              CG_2 ... CG_N
                                                       |                 |
    step 3                                             FG_11 ... FG_M1   FG_12 ... FG_M2 ...
                                                        \         |      |         /
                                                         \        \      |        /
    step 4                                                 _________ PLOTS ______
    
    Example of use:           
        .. code-block:: python
    
            YamboDG_Optimize(cg_grids,fg_grids,prefix,scf_path,pseudo_path,...,STEPS='all')

    """
    
    def __init__(self,cg_grids,fg_grids,prefix,nscf_input,ip_input,scf_save_path,pseudo_path,RUN_path='./',nscf_out="nscf",y_out_dir="results",qe_run_script=None,yambo_run_script=None,E_laser=0.,pw='pw.x',yambo='yambo',ypp='ypp',p2y='p2y',STEPS='all',RUN=True):
        
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
        
        #Executables
        self.pw = pw
        self.p2y = p2y
        self.ypp = ypp
        self.yambo = yambo

        #Start IO
        self.yf = YamboIO(out_name='YAMBOPY_double-grid_Optimize.log',out_path=self.RUN_path,print_to_shell=True)
        self.yf.IO_start()

        if STEPS=='1' or STEPS=='all': 
            self.setup_cg()
            if RUN: self.run_jobs(nscf_out,y_out_dir,qe_script=qe_run_script,yambo_script=yambo_run_script)
        
        if STEPS=='2' or STEPS=='all':
            for ig,cg in enumerate(self.cg_strings):
                calc_dir = '%s/%s_coarse_grid'%(self.nscf_dir,cg)
                calc = self.check_nscf_completed(calc_dir,nscf_out)
                if calc: self.yf.msg("NSCF CG %s found."%cg)
                else: self.yf.msg("NSCF CG %s NOT found."%cg)
                if calc:
                    yambo_dir = '%s/%s_coarse_grid'%(self.ip_dir,cg)
                    self.generate_SAVEs(calc_dir,yambo_dir)
                    self.setup_fg(calc_dir,yambo_dir,self.fg_grids[ig],self.fg_strings[ig])
            if RUN: self.run_jobs(nscf_out,y_out_dir,qe_script=qe_run_script,yambo_script=yambo_run_script)
                    
        if STEPS=='3' or STEPS=='all':
            for ig,cg in enumerate(self.cg_strings):
                for iff,fg in enumerate(self.fg_strings[ig]):
                    calc_dir = '%s/%s_coarse_grid/%s'%(self.nscf_dir,cg,fg)
                    calc = self.check_nscf_completed(calc_dir,nscf_out)
                    if calc: self.yf.msg("NSCF CG %s FG %s found."%(cg,fg))
                    else: self.yf.msg("NSCF CG %s FG %s NOT found."%(cg,fg))
                    if calc:
                        yambo_dir = '%s/%s_coarse_grid/%s'%(self.ip_dir,cg,fg)
                        if not os.path.isfile('%s/SAVE/ndb.Double_Grid'%yambo_dir):
                            self.generate_SAVEs(calc_dir,"%s/dg_SAVE"%yambo_dir)
                            self.setup_yambo_fg(yambo_dir,self.fg_grids[ig][iff],y_out_dir)
            if RUN: self.run_jobs(nscf_out,y_out_dir,qe_script=qe_run_script,yambo_script=yambo_run_script)
        
        if STEPS=='4' or STEPS=='all':
            yout = 'o-%s.eps_q1_ip'%y_out_dir
            ip_data = []
            ip_labels = []
            #
            for ig,cg in enumerate(self.cg_strings):
                ip_data_temp = []
                ip_labels_temp = []
                temp_path = '%s/%s_coarse_grid/%s/%s'%(self.ip_dir,cg,y_out_dir,yout)
                if os.path.isfile(temp_path):
                    self.yf.msg("IP CG %s found."%cg)
                    ip_data_temp.append(np.loadtxt(temp_path))
                    ip_labels_temp.append(cg)
                else:
                    self.yf.msg("IP CG %s NOT found."%cg)
                for iff,fg in enumerate(self.fg_strings[ig]):
                    temp_path = '%s/%s_coarse_grid/%s/%s/%s'%(self.ip_dir,cg,fg,y_out_dir,yout)
                    if os.path.isfile(temp_path):
                        self.yf.msg("IP CG %s FG %s found."%(cg,fg))
                        ip_data_temp.append(np.loadtxt(temp_path))
                        ip_labels_temp.append('%s N_FG=%d'%(cg,self.fg_grids[ig][iff]))
                    else:
                        self.yf.msg("IP CG %s FG %s NOT found."%(cg,fg))
                ip_data.append(ip_data_temp)
                ip_labels.append(ip_labels_temp)
            #
            ip_data = [ip for ip in ip_data if ip != []]
            ip_labels = [ip for ip in ip_labels if ip != []]
            self.plot_ip_spectra(ip_data,ip_labels,E_laser)
        
        self.yf.IO_close()
        
    def setup_cg(self):
        """ First step of the workflow: setup CG folder tree and CG nscf calculations
        """
        for ig,cg in enumerate(self.cg_grids):
            work_dir = "%s/%s_coarse_grid"%(self.nscf_dir,self.cg_strings[ig])
            if not os.path.isdir(work_dir):
                os.mkdir(work_dir)
                self.generate_inputfile(work_dir,cg)
                os.system('cp -r %s/%s.save %s'%(self.scf_save_path,self.prefix,work_dir))
            yambo_dir = "%s/%s_coarse_grid"%(self.ip_dir,self.cg_strings[ig])
            if not os.path.isdir(yambo_dir):
                os.mkdir(yambo_dir)
                self.ip_inp.write('%s/ip.in'%yambo_dir)
                
    def setup_fg(self,nscf_cg_dir,ip_cg_dir,fg_grid,fg_string):
        """ Second step of the workflow: setup FG folder tree and FG (CG) nscf (ip) calculations
        """
        for iff,fg in enumerate(fg_grid):
            fg_nscf_inp = '%s_fg.nscf'%self.prefix
            rand_nm = 'random_kpt_list_%s.dat'%fg_string[iff]
            ypp_inp = 'ypp_fg.in'
            work_dir = "%s/%s"%(nscf_cg_dir,fg_string[iff])
            if not os.path.isdir(work_dir):
                os.mkdir(work_dir)
                os.system('cp -r %s/%s.save %s'%(self.scf_save_path,self.prefix,work_dir))
            yambo_dir = "%s/%s"%(ip_cg_dir,fg_string[iff])
            if not os.path.isdir(yambo_dir):
                os.mkdir(yambo_dir)
                os.mkdir("%s/dg_SAVE"%yambo_dir)
                self.ip_inp.write('%s/ip.in'%yambo_dir)
            if not os.path.isfile('%s/%s'%(work_dir,fg_nscf_inp)):
                self.generate_ypp_input_random_grid(ip_cg_dir,fg,ypp_inp)
                ypp_run = self.scheduler()
                ypp_run.add_command('cd %s; %s -F %s > ypp_fg.log'%(ip_cg_dir,self.ypp,ypp_inp))
                ypp_run.add_command('mv o.random_k_pts %s'%rand_nm)
                ypp_run.add_command('cp %s %s'%(rand_nm,work_dir))
                ypp_run.run()
                kpts_rndm = np.loadtxt('%s/%s'%(work_dir,rand_nm))
                if len(kpts_rndm) != fg:
                    self.yf.msg("[WARNING] Actual fine grid number of kpts is %d instead of %d"%(len(kpts_rndm),fg))
                self.generate_inputfile(work_dir,fg,klist=kpts_rndm)

    def setup_yambo_fg(self,ip_fg_dir,fg_num,yresults_dir,clean_dg_saves=True):
        """ Third step of the workflow: map FG to CG and FG ip calculations
        """
        ypp_inp = 'ypp_map.in'
        os_run = self.scheduler()
        os_run.add_command('cd %s; cp ../%s/ndb.dipoles ../SAVE/ ; cp -r ../SAVE .'%(ip_fg_dir,yresults_dir))
        os_run.run()
        self.generate_ypp_input_map_grid(ip_fg_dir,fg_num,ypp_inp)
        ypp_run = self.scheduler()
        ypp_run.add_command('cd %s; %s -F %s > ypp_map.log'%(ip_fg_dir,self.ypp,ypp_inp))
        ypp_run.run()
        if os.path.isfile('%s/SAVE/ndb.Double_Grid'%ip_fg_dir):
            self.yf.msg("                    Double Grid ready.")
            if clean_dg_saves: os.system('rm -r %s/dg_SAVE/SAVE'%ip_fg_dir)
        else:
            self.yf.msg("                    Double Grid NOT ready.")
    
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
        yppin['FineGd_mode']='mixed'
        yppin['BZ_DbGd_Nk']=fg_num
        yppin.arguments.append('SkipCheck')
        yppin['FineGd_DB1_paths'] = ['./dg_SAVE']

        yppin.write('%s/%s'%(folder,inp_nm))
        
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
        
    def generate_SAVEs(self,folder,yambo_folder):
        """Generate yambo SAVE folders
        """
        if not os.path.isdir('%s/SAVE'%yambo_folder):
            self.yf.msg("        Generating SAVE...")
            p2y_run = self.scheduler()
            p2y_run.add_command('cd %s/%s.save; %s > p2y.log'%(folder,self.prefix,self.p2y))
            p2y_run.add_command('%s > yambo.log'%self.yambo)
            p2y_run.add_command('mv SAVE %s'%os.path.abspath(yambo_folder))
            p2y_run.run()
        else:
            self.yf.msg("        SAVE folder found.")
            
    def run_jobs(self,out_qe,out_yambo,qe_script=None,yambo_script=None):
        """ Automatic job submission
        """
        eps_out = 'o-%s.eps_q1_ip'%out_yambo
        
        for ig,cg in enumerate(self.cg_strings):

            temp_dir = '%s/%s_coarse_grid'%(self.nscf_dir,cg)
            check = self.check_nscf_completed(temp_dir,out_qe)
            if not check:
                self.yf.msg("Running NSCF CG %s..."%cg)      ############## Run NSCF CG
                if qe_script is None: self.submit_frontend(self.pw,temp_dir,out_qe,prefix=self.prefix)
                else: self.submit_cluster(temp_dir,qe_script)

            save_dir = '%s/%s_coarse_grid'%(self.ip_dir,cg)
            if os.path.isdir('%s/SAVE'%save_dir):
                eps_path = '%s/%s/%s'%(save_dir,out_yambo,eps_out)
                condition = ( not os.path.isfile(eps_path) ) or \
                            ( os.path.isfile(eps_path) and os.stat(eps_path).st_size == 0 )
                if condition:
                    self.yf.msg("Running IP CG %s..."%cg)     ############ Run IP CG
                    if yambo_script is None: self.submit_frontend(self.yambo,save_dir,out_yambo)
                    else: self.submit_cluster(save_dir,yambo_script)

            for fg in self.fg_strings[ig]:

                temp_dir = '%s/%s_coarse_grid/%s'%(self.nscf_dir,cg,fg)
                if os.path.isdir(temp_dir):
                    check = self.check_nscf_completed(temp_dir,out_qe)
                    if not check:
                        self.yf.msg("Running NSCF CG %s FG %s..."%(cg,fg)) ######### Run NSCF FG
                        if qe_script is None: self.submit_frontend(self.pw,temp_dir,out_qe,prefix=self.prefix)
                        else: self.submit_cluster(temp_dir,qe_script)

                save_dir = '%s/%s_coarse_grid/%s'%(self.ip_dir,cg,fg)
                if os.path.isdir('%s/SAVE'%save_dir) and os.path.isfile('%s/SAVE/ndb.Double_Grid'%save_dir):
                    eps_path = '%s/%s/%s'%(save_dir,out_yambo,eps_out)
                    condition = ( not os.path.isfile(eps_path) ) or \
                                ( os.path.isfile(eps_path) and os.stat(eps_path).st_size == 0 )
                    if condition:
                        self.yf.msg("Running IP CG %s FG %s..."%(cg,fg))     ############ Run IP FG
                        if yambo_script is None: self.submit_frontend(self.yambo,save_dir,out_yambo)
                        else: self.submit_cluster(save_dir,yambo_script)
    
    def submit_frontend(self,code,run_dir,out_file,prefix=None):
        """ Submit frontend calculation
        """
        run_job = self.scheduler()
        if code[-4:]=='pw.x': run_job.add_command('cd %s; %s -inp %s.nscf > %s.out'%(run_dir,code,self.prefix,out_file))
        if code[-4:]=='ambo': run_job.add_command('cd %s; %s -F ip.in -J %s -C %s'%(run_dir,code,out_file,out_file))
        run_job.run()

    def submit_cluster(self,run_dir,script_path):
        """ Submit job to cluster
        """
        script_name = script_path.rsplit('/',1)[-1]
        run_job = self.scheduler()
        if not os.path.isfile('%s/%s'%(run_dir,script_name)): run_job.add_command('cp %s %s/'%(script_path,run_dir))
        run_job.add_command('cd %s; sbatch %s'%(run_dir,script_name))
        run_job.run()
        
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
        rm_run = self.scheduler()
        rm_run.add_command('rm -rf %s %s %s'%(self.nscf_dir,self.ip_dir,self.plot_dir))
        rm_run.run()
        self.yf.msg("Workflow removed.")