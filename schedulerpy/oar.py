# Copyright (C) 2018 Henrique Pereira Coutada Miranda, Alejandro Molina-Sanchez
# All rights reserved.
#
# This file is part of yambopy
#
#
from __future__ import print_function
import subprocess
from .scheduler import Scheduler

class Oar(Scheduler):
    """
    Class to submit jobs through the OAR scheduler

    ``_vardict`` states the default assignement of the nodes and cores variables
    from the schduler class to the variables needed in this class
    """
    _vardict = {"cores":"core",
                "nodes":"nodes"}
                
    def initialize(self):
        self.get_vardict()
        args = self.arguments
        if self.get_arg("besteffort"): args.append("-t besteffort")
        if self.get_arg("idempotent"): args.append("-t idempotent")
        if self.get_arg("bigmem"): args.append("-t bigmem")
        
        if self.name:  args.append("-n \"%s\""%self.name)
        dependent = self.get_arg("dependent")
        if dependent: args.append("-a %d"%dependent)
        
        resources_line = self.get_resources_line()
        if resources_line:
            args.append(resources_line)
        
    def get_resources_line(self):
        """
        get the the line with the resources
        """
        s = "-l "
        if self.nodes: s += "%s=%d"%(self.vardict['nodes'],self.nodes)
        if self.nodes and self.cores: s += "/"
        if self.cores: s += "%s=%d"%(self.vardict['cores'],self.cores)
        if self.nodes or self.cores: s+= ","
        s += "walltime=%s"%self.walltime
        return s
    
    def get_script(self):
        """
        get a .pbs file to be submitted using qsub
        qsub <filename>.pbs
        """
        s = '#!/bin/bash\n'
        s += "\n".join(["#OAR %s"%s for s in self.arguments])+'\n'
        s += self.get_commands()
        return s
        
    def get_bash(self):
        """
        get a bash command to submit the job
        """
        command  = "oarsub "
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
