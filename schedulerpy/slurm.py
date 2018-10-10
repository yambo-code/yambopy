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
                "nodes":"nodes"}
                
    def initialize(self):
        self.get_vardict()

    def get_script(self):
        """
        get a .pbs file to be submitted using qsub
        qsub <filename>.pbs
        """
        lines = []; app = lines.append
        app('#!/bin/bash -l\n')
        
        partition = self.get_arg("partition",None)
        if self.name: app("#SBATCH -J \"%s\""%self.name)
        if partition: app('#SBATCH --partition %s'%partition)

        qos = self.get_arg("qos",None)
        if qos: app('#SBATCH --qos=%s'%qos)

        if self.nodes: app("#SBATCH -N %d" % self.nodes)
        if self.cores: app("#SBATCH --ntasks-per-node=%d" % self.cores)

        mem_per_cpu = self.get_arg("mem-per-cpu",None) 
        if mem_per_cpu: app("#SBATCH --mem-per-cpu=%d" % mem_per_cpu)

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
        #create the submission script
        self.write(filename)        
        workdir  = os.path.dirname(filename)
        basename = os.path.basename(filename) 
 
        if dry:
            print(command)
        else:
            p = subprocess.Popen([command,basename],stdout=subprocess.PIPE,stderr=subprocess.PIPE,cwd=workdir)
            self.stdout,self.stderr = p.communicate()
            
            #check if there is stderr
            if self.stderr: raise Exception(self.stderr)
            
            #check if there is stdout
            if verbose: print(self.stdout)
