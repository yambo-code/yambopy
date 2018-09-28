# Copyright (C) 2016 Henrique Pereira Coutada Miranda, Alejandro Molina-Sanchez
# All rights reserved.
#
# This file is part of yambopy
#
#
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
        args = self.arguments
         
        if self.name:  args.append("#SBATCH -J \"%s\"\n"%self.name)
        if self.get_arg("name"): args.append("")
        resources_line = self.get_resources_line()
        if resources_line:
            args.append(resources_line)
        
    def get_resources_line(self):
        """
        get the the line with the resources
        """
        s = " "
        if self.nodes: s = "#SBATCH -n %d\n"                % (self.nodes)
        #if self.ntask: s += "#SBATCH --ntasks-per-node=%d\n" % (self.ntask)
        #if self.memory: s = "#SBATCH --memory %d" %(self.memory)
        #if self.memory_per_cpu: s = "#SBATCH --memory-per-cpu %d" %(self.memory_per_cpu)
        s += "#SBATCH --time=0-%s" % self.walltime
        return s
    
    def get_script(self):
        """
        get a .pbs file to be submitted using qsub
        qsub <filename>.pbs
        """
        s = '#!/bin/bash -l \n'
        s += '#SBATCH -p batch \n'
        s += '#SBATCH --qos=qos-batch \n'
        s += self.get_commands()
        return s

    def get_bash(self):
        """
        get a bash command to submit the job
        """
        command  = "sbatch "
        command += " \\\n".join(self.arguments)+" "
        command += "\"%s\""%self.get_commands().replace("'","\"").replace("\"","\\\"")
        return command
        
    def __str__(self):
        """
        create the string for this job
        """
        return self.get_script()

    def run(self,silent=True,dry=False):
        """
        run the command
        arguments:
        dry - only print the commands to be run on the screen
        """
        command = self.get_bash()
        
        if dry:
            print(command)
        else:
            p = subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
            self.stdout,self.stderr = p.communicate()
            
            #check if there is stderr
            if self.stderr: raise Exception(self.stderr)
            
            #check if there is stdout
            if not silent: print(self.stdout)
                
            #get jobid
            for line in self.stdout.split('\n'):
                if 'OAR_JOB_ID' in line:
                    self.jobid = int(line.strip().split('=')[1])
            print("jobid:",self.jobid)
