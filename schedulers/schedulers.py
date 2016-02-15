# Copyright (C) 2015 Henrique Pereira Coutada Miranda, Alejandro Molina-Sanchez
# All rights reserved.
#
# This file is part of yambopy
#
#
import subprocess
#
# Scheduler of HPC cluster of the University of Luxembourg (Gaia and Chaos)
#
class oarsub():
    def __init__(self, nodes=None, core=None, be=False, dependent=0, idempotent=False, walltime="1:00:00", name=None, option=None ):
        self.nodes = nodes
        self.name = name
        self.core = core
        self.commands = []
        self.be = be
        self.dependent = dependent
        self.idempotent = idempotent
        self.walltime = walltime
        self.option = option

    def __str__(self):
        s = 'oarsub '
        if self.be:         s += " -t besteffort "
        if self.idempotent: s += " -t idempotent "
        if self.dependent:  s += " -a %d "%self.dependent
	if self.option:     s += " -p %s "%self.option
        if self.name:  s += "-n \"%s\""%self.name
        s += " -l \""        
        if self.nodes: s += "nodes=%d"%self.nodes
        if self.nodes and self.core: s += "/"
        if self.core:  s += "core=%d"%self.core
        if self.nodes or self.core: s+= ","
        s += "walltime=%s\" "%self.walltime
        s += """\"if [ -f  /etc/profile ]; then . /etc/profile; fi ; module use /work/projects/tss-physics/modules/ \n"""
        s += "\n".join(self.commands)
        s += "\"\n"
        return s

    def add_command(self,cmd):
        self.commands.append(cmd)

    def run(self):
        p = subprocess.Popen(str(self),stdout=subprocess.PIPE,shell=True)
        out = p.communicate()[0].split('\n')
        for line in out:
            print line
            if 'OAR_JOB_ID' in line:
                self.jobid = int(line.strip().split('=')[1])
        print "OAR_JOB_ID=%d"%self.jobid

    def clean(self):
        self.commands = []

    def write(self,filename=None):
        if not filename: filename = self.name
        f = open(filename,"w")
        f.write(str(self))
        f.close()
