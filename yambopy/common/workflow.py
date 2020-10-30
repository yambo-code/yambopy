import time
from schedulerpy import *

"""
This file contains the basic functions needed to check and manage workflows

    - wait_for_job: Let the python execution sleep until job completion
    - TODO: submit_job, check_for_job_completion, ... 
"""

def wait_for_job(self,shell,run_dir,time_step=10.):
    """
    Let the python execution sleep until job completion.
    
    - shell: schedulerpy object relative to submitted job
    - run_dir: directory where the job is being run
    - time_step: checking period (seconds)
    
    Scheduler types supported:
    
        - bash
        - slurm
    """
    job_status = shell.check_job_status(run_dir)
    condition = job_status=='R' or job_status=='PD' or job_status=='CG'
    while condition:
        time.sleep(time_step)
        job_status = shell.check_job_status(run_dir) 
        condition = job_status=='R' or job_status=='PD' or job_status=='CG'

def wait_for_setup_operations(filename,run_dir,time_step=10.):
    """
    To be used inside submission scripts which (i) have depedencies, (ii) need setup.

    Example: yambo calculation can only run after nscf is finished (dependency), but
             needs the SAVE folder to be created by the python master process as well (setup)

    When the dependencies are fulfilled, the job is submitted by the scheduler. This function delays
    the execution of the actual code (i.e., yambo in the example above) until the setup is also completed
    by the master python process. 

    This usually lasts a much shorter time than the actual execution. By executing this function inside the
    submission script we avoid hanging the workflow.

    - filename: an empty file (i.e., named 'SETUP_DONE') touched by the setup routine to signal end of setup.
                when filenames 
    """
    stderr = True
    while stderr:
        time.sleep(time_step)
        p = subprocess.Popen(['ls','%s'%filename],stdout=subprocess.PIPE,stderr=subprocess.PIPE,cwd=run_dir)
        stdout,stderr = p.communicate() 
    
    
