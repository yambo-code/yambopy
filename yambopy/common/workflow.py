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
