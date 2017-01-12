# Copyright (C) 2016 Henrique Pereira Coutada Miranda, Alejandro Molina-Sanchez
# All rights reserved.
#
# This file is part of yambopy
#
#
import subprocess
from schedulerpy import *

class Bash(Scheduler):
    """
    Class to submit jobs using bash
    """
    def initialize(self):
        pass
        
    def __str__(self):
        return self.get_commands()

    def get_bash(self):
        return str(self)

    def get_script(self):
        return str(self)

    def add_mpirun_command(self, cmd):
        threads = 1
        if self.cores: threads*=self.cores
        if self.nodes: threads*=self.nodes
        self.add_command("mpirun -np %d %s"%(threads,cmd))
        
    def run(self,dry=False):
        if dry:
            print str(self)
        else:
            p = subprocess.Popen(str(self),stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,executable='/bin/bash')
            self.stdout, self.stderr = p.communicate()
            if self.stderr != "":
                raise ValueError("ERROR:\n%s"%self.stderr)
            print self.stdout
