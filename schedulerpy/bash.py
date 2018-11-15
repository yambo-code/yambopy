# Copyright (C) 2018 Henrique Pereira Coutada Miranda, Alejandro Molina-Sanchez
# All rights reserved.
#
# This file is part of yambopy
#
from __future__ import print_function
import os
from builtins import str
import subprocess
import sys
from .scheduler import Scheduler

class Bash(Scheduler):
    """
    Class to submit jobs using BASH
    """
    _vardict = {"cores":"core",
                "nodes":"nodes"}

    def initialize(self):
        self.get_vardict()

    def __str__(self):
        return self.get_commands()

    def get_bash(self):
        return str(self)

    def get_script(self):
        return str(self)

    def add_mpirun_command(self, cmd):
        threads = 1
        if self.cores: threads*=self.cores
        mpirun = self.get_arg("mpirun","mpirun")
        np = self.get_arg("np","-np")
        self.add_command("%s %s %d %s"%(mpirun,np,threads,cmd))

    def run(self,filename='./run.sh',command='sh',dry=False):
        #create the submission script
        self.write(filename)
        workdir  = os.path.dirname(filename)
        basename = os.path.basename(filename)

        if dry:
            print(command)
        else:
            p = subprocess.Popen([command,basename],stdout=subprocess.PIPE,stderr=subprocess.PIPE,cwd=workdir)
            self.stdout, self.stderr = p.communicate()
            # In Python 3, Popen.communicate() returns bytes
            try:
                self.stdout = self.stdout.decode()
                self.stderr = self.stderr.decode()
            # If Python 2, <str>.decode() will raise an error that we ignore
            except AttributeError:
                pass
