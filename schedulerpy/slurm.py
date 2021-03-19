# Copyright (C) 2016 Henrique Pereira Coutada Miranda, Alejandro Molina-Sanchez
# All rights reserved.
#
# This file is part of yambopy
#
#
import os
import subprocess
from .scheduler import Scheduler

class Slurm(Scheduler):
    """
    Class to submit jobs through the Slurm scheduler

    ``_vardict`` states the default assignement of the nodes and cores variables
    from the schduler class to the variables needed in this class
    """
    _vardict = {"cores":"core",
                "nodes":"nodes",
                "cpus_per_task":"cpus_per_task"}
                          
                
    def initialize(self):
        self.get_vardict()

    def get_script(self):
        """
        get a .sh file to be submitted using sbatch
        sbatch <filename>.sh
        """
        lines = []; app = lines.append
        app('#!/bin/bash -l\n')
        
        partition = self.get_arg("partition",None)
        if self.name: app("#SBATCH -J %s"%self.name)
        if partition: app('#SBATCH --partition %s'%partition)

        qos = self.get_arg("qos",None)
        if qos: app('#SBATCH --qos=%s'%qos)

        dependency = self.get_arg("dependency",None)
        if dependency: app("#SBATCH --dependency=%s"%dependency)

        if self.nodes: app("#SBATCH -N %d" % self.nodes)
        if self.cores: app("#SBATCH --ntasks-per-node=%d" % self.cores)
        if self.cpus_per_task: app("#SBATCH --cpus-per-task=%d" % self.cpus_per_task )

        mem_per_cpu = self.get_arg("mem_per_cpu",None) 
        if mem_per_cpu: app("#SBATCH --mem-per-cpu=%s" % mem_per_cpu)

        app("#SBATCH --time=0-%s" % self.walltime)

        app(self.get_commands())
        return "\n".join(lines)

    def __str__(self):
        """
        create the string for this job
        """
        return self.get_script()

    def run(self,filename='run.sh',dry=False,command="sbatch",verbose=0):
        """
        Create the submission script and submit the job
        
        Arguments:
            dry - only print the commands to be run on the screen
        """
        if dry: 
            print(self)
            return 

        #create the submission script
        self.write(filename)        
        workdir  = os.path.dirname(filename)
        basename = os.path.basename(filename) 

        p = subprocess.Popen([command,basename],stdout=subprocess.PIPE,stderr=subprocess.PIPE,cwd=workdir)
        self.stdout,self.stderr = p.communicate()
        #Slurm-specific instruction to store jobid
        self.jobid = self.stdout.decode().split(' ')[-1].strip()
        
        #check if there is stderr
        if self.stderr: raise Exception(self.stderr)
        
        #check if there is stdout
        if verbose: print(self.stdout)
        
    def check_job_status(self,workdir):
        """
        Return status of slurm job (empty if job is not present)
        """
        p = subprocess.Popen(['squeue','-j %s'%self.jobid],stdout=subprocess.PIPE,stderr=subprocess.PIPE,cwd=workdir)
        stdout,stderr = p.communicate()
        if stdout:
            if stdout.decode()[-9:-1]=='(REASON)': job_status = 'NULL' 
            else: job_status = stdout.decode().split('\n')[1].split()[4]
        else: 
            job_status = 'NULL'
        return job_status
