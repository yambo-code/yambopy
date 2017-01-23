# Copyright (C) 2015 Henrique Pereira Coutada Miranda, Alejandro Molina-Sanchez
# All rights reserved.
#
# This file is part of yambopy
#
#
import subprocess
from textwrap import dedent
#
# Scheduler of HPC cluster of CENAERO
#
class pbs():
    def __init__(self, name='pbs', nodes=1, core=1, dependent=0, mem=2624, queue='main', group_list=None, walltime="1:00:00", option=None, modules=[] ):
        self.nodes = nodes
        self.core = core
        self.name = name
        self.dependent = dependent
        self.queue = queue
        self.mem = mem
        self.group_list = group_list
        self.commands = []
        self.walltime = walltime
        self.option = option
        self.modules = modules
        self.header = dedent('''
                             exec > ${PBS_O_WORKDIR}/${PBS_JOBNAME}_${PBS_JOBID}.log 
                             echo "------------------ Work dir --------------------" 
                             cd ${PBS_O_WORKDIR} && echo ${PBS_O_WORKDIR} 
                             echo "------------------ Job Info --------------------" 
                             echo "jobid : $PBS_JOBID" 
                             echo "jobname : $PBS_JOBNAME" 
                             echo "job type : $PBS_ENVIRONMENT" 
                             echo "submit dir : $PBS_O_WORKDIR" 
                             echo "queue : $PBS_O_QUEUE" 
                             echo "user : $PBS_O_LOGNAME" 
                             echo "threads : $OMP_NUM_THREADS"\n
                             ''')

    def __str__(self):
        s = '#!/bin/bash\n'
        if self.name:  s += "#PBS -N %s\n"%self.name
        s += "#PBS -l select=%d:ncpus=%d:mpiprocs=%d:vmem=%dmb:ompthreads=1\n"%(self.nodes,self.core,self.core,self.mem*self.core)
        s += "#PBS -l pvmem=%dmb\n"%(self.mem*self.core)
        s += "#PBS -q %s\n"%(self.queue)
        s += "#PBS -r y\n"
        if self.group_list: s += "#PBS -W group_list=%s\n"%self.group_list
        if self.dependent:  s += "#PBS -W depend=afterok:%s\n"%self.dependent
        s += "#PBS -l walltime=%s\n"%self.walltime
        s += dedent(self.header)
        for names in self.modules:
          s += "module load %s\n" % (names)
        s += "\n".join(self.commands)
        return s

    def add_command(self,cmd):
        self.commands.append(cmd)

    def run(self):
	self.write('%s.sh'%self.name)
        p = subprocess.Popen('qsub %s.sh'%self.name,stdout=subprocess.PIPE,shell=True)
	self.jobid = p.communicate()[0].split('\n')[0]
	print self.jobid
	return self.jobid

    def clean(self):
        self.commands = []

    def write(self,filename=None):
        if not filename: filename = self.name
        f = open(filename,"w")
        f.write(str(self))
        f.close()
