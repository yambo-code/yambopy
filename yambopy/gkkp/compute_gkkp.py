#
# License-Identifier: GPL
#
# Copyright (C) 2024 The Yambo Team
#
# Authors: FP
#
# This file is part of the yambopy project
#
from yambopy.common.save_generation import CreateYamboSave
from yambopy.common.calculation_manager import check_qe_completed, shell_qe_run
from yambopy.io.iofile import YamboIO
from schedulerpy import Scheduler
import os
from copy import deepcopy

class YamboGkkpCompute():
    """
    Class to obtain qe s.dbph* and yambo ndb.elph* databases starting from scratch.

    [NOTE: this class uses the old el-ph yambo interface (via double ph.x calculation + ypp_ph)]
    [      A new interface making use of the LetzElphC code is available                       ]
    
    It runs the necessary pw.x and ph.x simulations, optionally followed by the yambo setup.
    
    Inputs needed:
        <required> 
        - scf_input: pw scf input file [NOTE: dvscf and gkkp inputs are automatically generated, check parameters if interested]

        <optional>
        - work_dir: directory where flow is run and yambo SAVE will appear
        - pw_exec_path: path to executables
        - qe_scheduler: optional scheduler for cluster submission
        - with_SAVE: if True, workflow will generate yambo SAVE at the end (the python master process will remain active).
                     The workflow can be called a second time switching with_SAVE to True to immediately generate the SAVE.

          [SUBOPTIONS for with_SAVE]

              -- expand: whether to expand the matrix elements
              -- yambo_exec_path: path to executables

    TODO: allow for random qpoints in ph calculations
    """
    
    def __init__(self,scf_input,work_dir='.',pw_exec_path='',qe_scheduler=None,with_SAVE=False,yambo_exec_path='',expand=False):

        if not os.path.isdir(work_dir): os.mkdir(work_dir)
        self.RUN_path = os.path.abspath(work_dir)
        self.wait_up = with_SAVE #Slightly restructure dependencies and waith for job completions if SAVE is to be created

        #Configuring schedulers
        if qe_scheduler is not None: self.qejobrun= qe_scheduler                          # Here we use, e.g., slurm 
        else:                        self.qejobrun = Scheduler.factory(scheduler="bash")  # Run without submission

        #Executables
        if yambo_exec_path != '': yambo_exec_path+='/'
        self.yambo     = yambo_exec_path + 'yambo'
        self.yambo_ph  = yambo_exec_path + 'yambo_ph' 
        self.p2y       = yambo_exec_path + 'p2y'
        self.ypp_ph    = yambo_exec_path + 'ypp_ph'
        if pw_exec_path != '': pw_exec_path+='/'
        self.pw = pw_exec_path + 'pw.x'
        self.ph = pw_exec_path + 'ph.x'
        self.dynmat = pw_exec_path + 'dynmat.x'

        # Inputs
        self.scf_input = scf_input
        self.prefix    = scf_input.prefix
        
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
            self.yf.msg('Running scf.')
            self.run_scf()
        
        if not self.dvscf_status:
            self.yf.msg('Running dvscf.')
            self.run_dvscf()

        if not self.gkkp_status:
            self.yf.msg('Running gkkp.')
            self.run_gkkp()
        
        if not self.nscf_status:
            self.yf.msg('Running nscf.')
            self.run_nscf()
            
        # [OPTIONAL] Create SAVE
        if with_SAVE:
            self.setup_SAVE()

            if not self.are_gkkp_there:
                if self.is_SAVE_there: 
                     import shutil
                     shutil.rmtree('%s/SAVE'%self.RUN_path)
                if expand: save_type = 'expanded_elph'
                else:      save_type = 'elph'
            
                self.yf.msg('---- Generating SAVE folder: ----') 
                CreateYamboSave(self.prefix,save_type=save_type,nscf=self.nscf_dir,elph_path=self.gkkp_dir,database=self.RUN_path,\
                                yambo_exec_path=yambo_exec_path,printIO=True)
                self.clean_rubbish()

        #End IO        
        self.yf.IO_close()
        
    def setup_calculations(self):
        """
        Generate workflow tree
        """
        
        # Directory names (hardcoded)
        dft_dir  = '%s/dft'%self.RUN_path
        self.scf_dir  = '%s/scf'%dft_dir
        self.gkkp_dir = '%s/gkkp'%dft_dir
        self.nscf_dir = '%s/nscf'%dft_dir
        
        # Logicals
        self.gkkp_status  = False
        self.dvscf_status = False
        self.scf_status   = False
        self.nscf_status  = False
        
        if not os.path.isdir(dft_dir):  os.mkdir(dft_dir)
        if not os.path.isdir(self.scf_dir):  os.mkdir(self.scf_dir)
        if not os.path.isdir(self.gkkp_dir): os.mkdir(self.gkkp_dir)
        if not os.path.isdir(self.nscf_dir): os.mkdir(self.nscf_dir)
        
        # Check if any gkkp->dvscf->scf calculations have been already done
        self.gkkp_status = check_qe_completed(self.gkkp_dir,self.prefix,self.out_gkkp,calc_type='gkkp')
        if self.gkkp_status: 
            self.yf.msg('gkkp calculation found!') 
            self.dvscf_status = True
        else:
            self.dvscf_status = check_qe_completed(self.gkkp_dir,self.prefix,self.out_dvscf,calc_type='ph')
            if self.dvscf_status: 
                self.yf.msg('dvscf calculation found!')
            else:
                self.scf_status = check_qe_completed(self.scf_dir,self.prefix,self.out_scf,calc_type='pw')
                if self.scf_status: 
                    self.yf.msg('scf calculation found!')
                       
        # Check if any nscf->scf calculations have been already done
        self.nscf_status = check_qe_completed(self.nscf_dir,self.prefix,self.out_nscf,calc_type='pw')
        if self.nscf_status: 
             self.yf.msg('nscf calculation found')
             self.scf_status = True
        else:
             if not self.scf_status:
                  self.scf_status = check_qe_completed(self.scf_dir,self.prefix,self.out_scf,calc_type='pw')
                  if self.scf_status: self.yf.msg('scf calculation found!')

    def setup_SAVE(self):
        """
        Expand the workflow tree to include yambo SAVE
        """

        # Check if SAVE and/or gkkp dbs are there already
        save_dir = '%s/SAVE'%self.RUN_path
        if os.path.isdir(save_dir):
            self.yf.msg('SAVE folder found!')
            self.is_SAVE_there = True
            if os.path.isfile('%s/ndb.elph_gkkp'%save_dir) or os.path.isfile('%s/ndb.elph_gkkp_expanded'%save_dir):
                self.yf.msg('ndb.elph databases already found!')
                self.are_gkkp_there = True
            else:
                self.are_gkkp_there = False
        else:
            self.is_SAVE_there  = False
            self.are_gkkp_there = False
            
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
        jname = 'scf'
        self.scf_id = shell_qe_run(jname,inp_name,self.out_scf,self.scf_dir,exec=self.pw,scheduler=self.qejobrun)
        
    def run_dvscf(self):
        """
        Run dvscf calculation
        """
        # Generate and write down input
        dvscf_input = self.generate_ph_input('dvscf')
        inp_name = self.prefix + '.dvscf'
        dvscf_input.write('%s/%s'%(self.gkkp_dir,inp_name))
        
        # Generate and write down dynmat input
        dynmat_input = self.generate_dynmat_input()
        dynp_name = self.prefix + '.dynmat'
        dynmat_input.write('%s/%s'%(self.gkkp_dir,dynp_name))

        # Set dynmat run after completion of main task
        dyn_run = ["mpirun -np 1 %s -inp %s > dynmat.out"%(self.dynmat,dynp_name)]

        # Create symlink to qe save if needed
        commands = []
        if not os.path.islink('%s/%s.save'%(self.gkkp_dir,self.prefix)):
            commands.append('ln -s %s/%s.save %s/'%(self.scf_dir,self.prefix,self.gkkp_dir)) 

        # Manage dependency
        if self.scf_status: depend = None # No dependency if scf was found
        else: depend = self.scf_id
        
        # Submit calculation
        jname = 'dvscf'
        self.dvscf_id = shell_qe_run(jname,inp_name,self.out_dvscf,self.gkkp_dir,exec=self.ph,shell_name='dvscf',\
                                     scheduler=self.qejobrun,pre_run=commands,pos_run=dyn_run,depend_on_JOBID=depend) 
                                     
    def run_gkkp(self):
        """
        Run gkkp calculation
        """
        # Generate and write down input
        gkkp_input = self.generate_ph_input('gkkp')
        inp_name = self.prefix + '.gkkp'
        gkkp_input.write('%s/%s'%(self.gkkp_dir,inp_name))
        
        # Create symlink to qe save if needed
        commands = []
        if not os.path.islink('%s/%s.save'%(self.gkkp_dir,self.prefix)):
            commands.append('ln -s %s/%s.save %s/'%(self.scf_dir,self.prefix,self.gkkp_dir)) 
        
        # Manage dependency
        if self.dvscf_status: depend = None # No dependency if dvscf was found
        else: depend = self.dvscf_id

        # Submit calculation
        jname = 'gkkp'
        self.gkkp_id = shell_qe_run(jname,inp_name,self.out_gkkp,self.gkkp_dir,exec=self.ph,shell_name='gkkp',\
                                    scheduler=self.qejobrun,pre_run=commands,depend_on_JOBID=depend)    
                                    
    def run_nscf(self):
        """
        Run nscf calculation
        """
        # Generate and write down input
        nscf_input = self.generate_nscf_input()
        inp_name = self.prefix + '.nscf'
        nscf_input.write('%s/%s'%(self.nscf_dir,inp_name))
        
        # Create symlink to qe save if needed
        commands = []
        if not os.path.isdir('%s/%s.save'%(self.nscf_dir,self.prefix)):
            commands.append('cp -r %s/%s.save %s/'%(self.scf_dir,self.prefix,self.nscf_dir)) 
                
        # Dependency here may include gkkp job to ensure that this is the last job to be completed if SAVE is to be generated
        if self.wait_up:
            if self.gkkp_status and self.scf_status:       depend = None # No dependency if scf and gkkp were found
            elif self.gkkp_status and not self.scf_status: depend = self.gkkp_id # scf was found, not gkkp
            elif not self.gkkp_status and self.scf_status: depend = self.scf_id  # gkkp was found, not scf
            else:                                          depend = '%s:%s'%(self.scf_id,self.gkkp_id) # double dependency
        else: 
            if self.scf_status: depend = None # No dependency if scf was found
            else:               depend = self.scf_id    
        
        # Submit calculation
        jname = 'nscf'
        self.nscf_id = shell_qe_run(jname,inp_name,self.out_nscf,self.nscf_dir,exec=self.pw,scheduler=self.qejobrun,\
                                    pre_run=commands,depend_on_JOBID=depend,hang_python=self.wait_up)       

    def generate_nscf_input(self):
        """
        Create nscf input for yambo SAVE starting from scf input
        """
        
        nscf_input = deepcopy(self.scf_input)
        nscf_input.control['calculation']="'nscf'"
        nscf_input.electrons['diago_full_acc'] = ".true."
        nscf_input.electrons['conv_thr'] = 1e-8
        nscf_input.system['force_symmorphic'] = ".true."
        
        return nscf_input

    def generate_ph_input(self,mode):
       """
       Create dvscf or gkkp input starting from scf input
 
       - mode: either 'dvscf' or 'gkkp'
       """

       from qepy import PhIn 
       ph_input = PhIn()
       # Common to dvscf and gkkp
       ph_input['prefix'] = "'%s'"%self.prefix
       ph_input['fildyn'] = "'%s'"%(self.prefix+'.dyn')
       nq1,nq2,nq3 = [ int(nk) for nk in self.scf_input.kpoints ]
       ph_input.set_nq(nq1,nq2,nq3)
       ph_input['tr2_ph'] = 1e-14
       ph_input['fildvscf']="'dvscf'"
       ph_input['ldisp']='.true.'
       ph_input['qplot']='.false.'
       # Only dvscf
       if mode=='dvscf':

           # Add effective charges if dealing with a non-metal
           is_insulator = 'occupations' not in self.scf_input.system or self.scf_input.system['occupations'] != "'smearing'"
           if is_insulator: ph_input['epsil']='.true.'
           self.is_insulator = is_insulator

           ph_input['electron_phonon']="'dvscf'"
           ph_input['recover']='.true.'
           ph_input['trans']='.true.'

       elif mode=='gkkp':
           ph_input['electron_phonon']="'yambo'"
           ph_input['trans']='.false.'
       else: raise ValueError("ph input mode not recognized (either 'dvscf' or 'gkkp')")

       return ph_input

    def generate_dynmat_input(self):
       """
       Create dynmat input file in order to treat issues with the frequencies at the Gamma point:
         -- Apply the acoustic sum rule
         -- Correct the LO mode (if non-metal)

         Outputs are saved in the gkkp folder as prefix.GAMMA_eigs_eivs and prefix.GAMMA_eigs_norm_eivs
       """
       from qepy import DynmatIn
       dm_input = DynmatIn()
       dm_input['asr']="'crystal'"
       dm_input['fildyn']="'%s.dyn1'"%self.prefix
       dm_input['fileig']="'%s.GAMMA_eigs_eivs'"%self.prefix
       dm_input['filout']="'%s.GAMMA_eigs_norm_eivs'"%self.prefix

       # Add LO correction along first cartesian axis if dealing with non-metal
       if self.is_insulator:
           dm_input['q(1)']=1
           dm_input['q(2)']=0
           dm_input['q(3)']=0

       return dm_input

    def clean_rubbish(self):
       """
       Remove logs, reports and inputs generated during SAVE creation
       """
       from glob import glob
       run_dir = self.RUN_path+'/'
       logs1 =    glob(run_dir+'l-*')
       logs2 =    glob(run_dir+'l_*')
       reports1 = glob(run_dir+'r-*')
       reports2 = glob(run_dir+'r_*')
       setups   = glob(run_dir+'setup.in*')
       for log in logs1:       os.remove(log)
       for log in logs2:       os.remove(log)
       for report in reports1: os.remove(report)
       for report in reports2: os.remove(report)
       for setup in setups:    os.remove(setup)
       os.remove(run_dir+'gkkp.in')
     
