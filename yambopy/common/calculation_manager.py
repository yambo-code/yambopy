from yambopy.common.workflow import wait_for_job
import subprocess
import os
from schedulerpy import *
from copy import deepcopy

"""
This file contains the basic functions needed to check on calculations for 
completion and more.

TODO: Include a shell_run function for all executables

"""
def shell_qe_run(job_name,inp_name,out_name,run_dir,exec='pw.x',shell_name='qe',scheduler=None,depend_on_JOBID=None,hang_python=False,pre_run=[],pos_run=[]):
    """ 
    Submit QUANTUM ESPRESSO job
    
        exec: executable to be run with full path.
              options: /path/to/pw.x, /path/to/ph.x
        
        job_name: job name
        shell_name: name of *.sh script which is generated
        depend_on_JOBID: job id of simulation that the present job has a dependency on
        run_dir: where job is run
        out_name: name of output file
        inp_name: name of input file
        scheduler: instance of scheduler class (if not present, bash is initialised)
        pre_run / pos_run: LISTS containing commands to be added before / after the mpirun command
        hang_python: if True, python process sleeps until job is completed
        
        returns id of present submitted job (-1 if scheduler is bash)         
    """
    # Check executable
    if not exec[-4:]=='pw.x' and not exec[-4:]=='ph.x':
        raise ValueError('Executable not recognised (pw.x and ph.x are the only options).')
        
    # Copy scheduler instance in order to safely edit it
    if scheduler is None: shell = Scheduler.factory(scheduler="bash")
    else:                 shell = deepcopy(scheduler)
    
    shell.name = '%s_%s'%(job_name,shell.name)
    
    # Add dependency if specified
    if depend_on_JOBID is not None and shell.schedulertype != 'bash':
        dependency='afterok:%s'%depend_on_JOBID
        shell.kwargs['dependency']=dependency
        
    # Add additional commands if present
    if len(pre_run) != 0:
        for command in pre_run: shell.pre_run.append(command)

    # Additional commands to be executed after the main mpirun command
    if len(pos_run) != 0:   
        for command in pos_run: shell.pos_run.append(command)
 
    # Main mpirun command
    shell.add_mpirun_command('%s -inp %s > %s'%(exec,inp_name,out_name))
    shell.run(filename='%s/%s.sh'%(run_dir,shell_name)) ### Specify run path
    
    # Manage submissions if specified
    if hang_python: wait_for_job(shell,run_dir)
    if shell.schedulertype!='bash': this_job_id = shell.jobid
    else:                           this_job_id = -1
    
    shell.clean()
    
    return this_job_id

def check_qe_completed(folder,prefix,output_file,calc_type='pw'):
    """ 
    Check if qe calculation has correctly completed.
    
    - folder: where the calculation has been run.
    - prefix: qe prefix
    - output_file: name of output file
    - calc_type: either 'pw' or 'ph' or 'gkkp'
    
    """
    status = True
    
    # If save folder does not exist, return False (= NOT completed) immediately
    if calc_type=='pw' and not os.path.isdir('%s/%s.save'%(folder,prefix)):
        status = False
        return status
    elif calc_type=='ph' and not os.path.isdir('%s/_ph0'%folder):
        status = False
        return status
    elif calc_type=='gkkp' and not os.path.isdir('%s/elph_dir'%folder):
        status = False
        return status
    if calc_type != 'pw' and calc_type != 'ph' and calc_type != 'gkkp':
        raise ValueError("calc_type not recognised: it has to be either 'pw' or 'ph' or 'gkkp'.")
        
    # Next, check if output is correctly completed
    try:
        check = subprocess.check_output("grep JOB %s/%s*"%(folder,output_file), shell=True, stderr=subprocess.STDOUT)
        check = check.decode('utf-8')
        check = check.strip().split()[-1]
    except subprocess.CalledProcessError as e:
        check = ""
    if check != "DONE.": status = False
    return status
