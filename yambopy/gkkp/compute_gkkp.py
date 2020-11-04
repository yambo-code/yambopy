from yambopy import *
from qepy import *
from schedulerpy import *
import os
from copy import deepcopy

class YamboGkkpCompute():
    """
    Class to obtain yambo ndb.elph* databases starting from scratch.
    
    It runs the necessary pw.x and ph.x simulations, followed by the yambo setup.
    
    Inputs needed:
        - work_dir: directory where flow is run and yambo SAVE will appear
        - scf_input:   pw scf input file
        - dvscf_input: ph dvscf input file
        - gkkp_input:  ph gkkp input file
        - expand:      whether to expand the matrix elements
        - pw_exec_path, yambo_exec_path: paths to executables
        - qe_scheduler: optional scheduler for cluster submission
        - wait_up: if True, during cluster submission python process remains active and creates SAVE at the end
                   if False, python process exits after submitting jobs and class must be called again to generate SAVE later
        
    TODO: allow for random qpoints in ph calculations
    """
    
    def __init__(self,scf_input,dvscf_input,gkkp_input,work_dir='.',pw_exec_path='',yambo_exec_path='',qe_scheduler=None,expand=False,wait_up=False):

        self.RUN_path = work_dir

        #Configuring schedulers
        self.frontend = Scheduler.factory(scheduler="bash")
        if qe_scheduler is not None: #Here we use, e.g., slurm 
            self.qejobrun = qe_scheduler
            self.wait_up = wait_up
        else: 
            self.qejobrun = Scheduler.factory(scheduler="bash")
            self.wait_up = False

        #Executables
        if yambo_exec_path != '': yambo_exec_path+='/'
        self.yambo     = yambo_exec_path + 'yambo'
        self.yambo_ph  = yambo_exec_path + 'yambo_ph' 
        self.p2y       = yambo_exec_path + 'p2y'
        self.ypp_ph    = yambo_exec_path + 'ypp_ph'
        if pw_exec_path != '': pw_exec_path+='/'
        self.pw = pw_exec_path + 'pw.x'
        self.ph = pw_exec_path + 'ph.x'

        # Inputs
        self.scf_input   = scf_input
        self.dvscf_input = dvscf_input
        self.gkkp_input  = gkkp_input
        self.prefix = self.scf_input['prefix']
        
        # Output names
        self.out_scf   = 'scf.out'
        self.out_nscf  = 'nscf.out'
        self.out_dvscf = 'dvscf.out'
        self.out_gkkp  = 'gkkp.out'
        
        #Start IO
        self.yf = YamboIO(out_name='YAMBOPY_gkkp_calculation.log',out_path=self.RUN_path,print_to_shell=True)
        self.yf.IO_start()
        
        self.yf.msg('#### GKKP WORKFLOW ####')

        # Create folder structure
        self.setup_calculations()

        # Run jobs
        if not self.scf_status:
            self.run_scf()
            self.yf.msg('Running scf.')
        
        if not self.dvscf_status:
            self.run_dvscf()
            self.yf.msg('Running dvscf.')

        if not self.gkkp_status:
            self.run_gkkp()
            self.yf.msg('Running gkkp.')
        
        if not self.nscf_status:
            self.run_nscf()
            self.yf.msg('Running nscf.')
            
        # Create SAVE
        if not self.are_gkkp_there:
            if self.is_SAVE_there: os.rmdir('%s/SAVE'%self.RUN_path)
            if expand: save_type = 'expanded_elph'
            else:      save_type = 'elph'
            
            CreateYamboSave(self.prefix,save_type=save_type,nscf=self.nscf_dir,elph_path=self.gkkp_dir,database=self.RUN_path,\
                            yambo_exec_path=yambo_exec_path,printIO=False)

        #End IO        
        self.yf.IO_close()
        
    def setup_calculations(self):
        """
        Generate workflow tree
        """
        
        # Directory names (hardcoded)
        self.dft_dir  = '%s/dft'%self.RUN_path
        self.scf_dir  = '%s/scf'%dft
        self.gkkp_dir = '%s/gkkp'%dft
        self.nscf_dir = '%s/nscf'%dft
        
        # Logicals
        self.gkkp_status  = False
        self.dvscf_status = False
        self.scf_status   = False
        self.nscf_status  = False
        
        # Check if SAVE and/or gkkp dbs are there already
        if os.path.isdir('%s/SAVE'%self.RUN_path):
            self.yf.msg('SAVE folder found!')
            self.is_SAVE_there = True
            if os.path.isfile('%s/SAVE/ndb.elph_gkkp') or os.path.isfile('%s/SAVE/ndb.elph_gkkp_expanded'):
                self.yf.msg('ndb.elph databases already found!')
                self.are_gkkp_there = True
            else:
                self.are_gkkp_there = False
        else:
            self.is_SAVE_there  = False
            self.are_gkkp_there = False
            
        if not self.are_gkkp_there:
            if not os.path.isdir(self.dft_dir):  os.mkdir(self.dft_dir)
            if not os.path.isdir(self.scf_dir):  os.mkdir(self.scf_dir)
            if not os.path.isdir(self.gkkp_dir): os.mkdir(self.gkkp_dir)
            if not os.path.isdir(self.nscf_dir): os.mkdir(self.nscf_dir)
            
            # Check if any qe calculations have been already done
            self.gkkp_status = check_qe_completed(self.gkkp_dir,self.prefix,self.out_gkkp,calc_type='gkkp')
            if self.gkkp_status: 
                self.yf.msg('gkkp calculation found!') 
            else:
                self.dvscf_status = check_qe_completed(self.gkkp_dir,self.prefix,self.out_dvscf,calc_type='ph')
                if self.dvscf_status: 
                    self.yf.msg('dvscf calculation found!')
                else:
                    self.scf_status = check_qe_completed(self.scf_dir,self.prefix,self.out_scf,calc_type='pw')
                    if self.scf_status: 
                        self.yf.msg('scf calculation found!')
                           
            self.nscf_status = check_qe_completed(self.nscf_dir,self.prefix,self.out_nscf,calc_type='pw')
            if self.nscf_status: self.yf.msg('nscf calculation found')
            
    def run_scf(self):
        """
        Run scf calculation
        """
        if self.scf_input.system['nbnd'] is None:
            raise ValueError('Please specify nbnd in the scf input in order to be able to compute the gkkp elements.')
        
        # Write down input
        inp_name = self.prefix + '.scf'
        self.scf_input.write('%s/%s'%(self.scf_dir,inp_name))
        
        # Submit calculation
        jname = 'scf_'+self.prefix
        self.scf_id = shell_qe_run(jname,inp_name,self.out_scf,self.scf_dir,exec=self.pw,scheduler=self.qejobrun)
        
    def run_dvscf(self):
        """
        Run dvscf calculation
        """
        # Write down input
        inp_name = self.prefix + '.dvscf'
        self.dvscf_input.write('%s/%s'%(self.gkkp_dir,inp_name))
        
        # Create symlink to qe save if needed
        commands = []
        if not os.path.isfile('%s/%s.save'%(self.gkkp_dir,self.prefix)):
            commands.append('ln -s %s/%s.save %s/'%(self.scf_dir,self.prefix,self.gkkp_dir)) 
        
        # Submit calculation
        jname = 'dvscf_'+self.prefix
        self.dvscf_id = shell_qe_run(jname,inp_name,self.out_dvscf,self.gkkp_dir,exec=self.ph,scheduler=self.qejobrun,\
                                     commands=commands,depend_on_JOBID=self.scf_id)   
                                     
    def run_gkkp(self):
        """
        Run gkkp calculation
        """
        # Write down input
        inp_name = self.prefix + '.gkkp'
        self.gkkp_input.write('%s/%s'%(self.gkkp_dir,inp_name))
        
        # Create symlink to qe save if needed
        commands = []
        if not os.path.isfile('%s/%s.save'%(self.gkkp_dir,self.prefix)):
            commands.append('ln -s %s/%s.save %s/'%(self.scf_dir,self.prefix,self.gkkp_dir)) 
        
        # Submit calculation
        jname = 'gkkp_'+self.prefix
        self.gkkp_id = shell_qe_run(jname,inp_name,self.out_gkkp,self.gkkp_dir,exec=self.ph,scheduler=self.qejobrun,\
                                    commands=commands,depend_on_JOBID=self.dvscf_id)    
                                    
    def run_nscf(self):
        """
        Run gkkp calculation
        """
        # Write down input
        nscf_input = self.generate_nscf_input()
        inp_name = self.prefix + '.nscf'
        nscf_input.write('%s/%s'%(self.nscf_dir,inp_name))
        
        # Create symlink to qe save if needed
        commands = []
        if not os.path.isfile('%s/%s.save'%(self.nscf_dir,self.prefix)):
            commands.append('ln -s %s/%s.save %s/'%(self.scf_dir,self.prefix,self.nscf_dir)) 
                
        # Dependency here may include gkkp job to ensure that this is the last job to be completed if SAVE is to be generated
        if self.wait_up: deps = '%d:%d'%(self.scf_id,self.gkkp_id) 
        else: deps = self.scf_id    
        
        # Submit calculation
        jname = 'gkkp_'+self.prefix
        self.gkkp_id = shell_qe_run(jname,inp_name,self.out_nscf,self.nscf_dir,exec=self.pw,scheduler=self.qejobrun,\
                                    commands=commands,depend_on_JOBID=deps,hang_python=self.wait_up)       

    def generate_nscf_input(self):
        """
        Create nscf input for yambo SAVE starting from scf input
        """
        
        nscf_input = deepcopy(self.scf_input)
        nscf_input.control['calculation']="'nscf'"
        nscf_input.electrons['diago_full_acc'] = ".true."
        nscf_input.electrons['conv_thr'] = 1e-8
        nscf_input.electrons['force_symmorphic'] = ".true."
        
        return nscf_input