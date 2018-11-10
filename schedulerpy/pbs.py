# Copyright (C) 2018 Henrique Pereira Coutada Miranda, Alejandro Molina-Sanchez
# All rights reserved.
#
# This file is part of yambopy
#
#
from __future__ import print_function
import os
import subprocess
from textwrap import dedent
from copy import deepcopy
from collections import OrderedDict
from .scheduler import Scheduler

class Pbs(Scheduler):
    """
    Class to submit jobs through the PBS scheduler.

    ``_vardict`` states the default assignement of the nodes and cores variables
    from the schduler class to the variables needed in this class
    """
    _vardict = {"cores":"core",
                "nodes":"select"}
                   
    def initialize(self):
        self.get_vardict()
        args = self.arguments
        queue = self.get_arg("queue")
        rerunable= self.get_arg("rerunable")
        mem = self.get_mem()
        if self.name: args.append("-N %s"%self.name)
        
        if queue: args.append("-q %s"%(queue))
        group_list = self.get_arg("group_list")
        
        if group_list: args.append("-W group_list=%s"%group_list)
        dependent = self.get_arg("dependent")
        
        if dependent: args.append("-W depend=afterok:%s"%dependent)
        args.append("-l walltime=%s"%self.walltime)

        if rerunable: args.append("-r y")

        if self.get_arg("pvmem"): args.append("-l pvmem=%dMB"%mem)
        
        resources_line = self.get_resources_line()
        if resources_line:
            args.append("-l %s"%resources_line)

    def get_mem(self):
        """
        get the memory for this job
        """

        #block to evaluate expressions from
        #http://stackoverflow.com/questions/2371436/evaluating-a-mathematical-expression-in-a-string
        import ast
        import operator as op

        # supported operators
        operators = {ast.Add: op.add, ast.Sub: op.sub, ast.Mult: op.mul,
                            ast.Div: op.truediv, ast.Pow: op.pow, ast.BitXor: op.xor,
                                         ast.USub: op.neg}

        def eval_expr(expr):
            return eval_(ast.parse(expr, mode='eval').body)

        def eval_(node):
            if isinstance(node, ast.Num): # <number>
                return node.n
            elif isinstance(node, ast.BinOp): # <left> <operator> <right>
                return operators[type(node.op)](eval_(node.left), eval_(node.right))
            elif isinstance(node, ast.UnaryOp): # <operator> <operand> e.g., -1
                return operators[type(node.op)](eval_(node.operand))
            else:
                raise TypeError(node)
        ######

        mem = self.get_arg("mem")
        if mem:
            if self.cores: cores = self.cores
            else: cores = 1
            if self.nodes: nodes = self.nodes
            else: nodes = 1
            mem = mem.replace("nodes",str(nodes))
            mem = mem.replace("cores",str(cores))
            mem = eval_expr(mem) 
        return mem 

    def get_resources_line(self):
        """
        get the the line with the resources
        """
        tags = ['select','nodes','core','ppn','ncpus','mpiprocs','ompthreads']
        args = [self.get_arg(tag) for tag in tags]
        resources = []
        if self.nodes: resources.append((self.vardict['nodes'],self.nodes))
        if self.cores: resources.append((self.vardict['cores'],self.cores))
        resources += [(tag,value) for tag,value in zip(tags,args) if value is not None]
        resources = OrderedDict(resources)

        # memory stuff
        mem = self.get_mem()
        if mem: resources["mem"]  = "%dMB"%mem
        
        resources_line = ":".join(["%s=%s"%(item,value) for item,value in list(resources.items())])
       
        return resources_line
    
    def get_script(self):
        """
        get a .pbs file to be submitted using qsub
        qsub <filename>.pbs
        """
        s = '#!/bin/bash\n'
        s += "\n".join(["#PBS %s"%s for s in self.arguments])+'\n'
        s += self.get_commands()
        return s
        
    def get_bash(self):
        """
        get a bash command to submit the job
        """
        s = "echo \"%s\" | "%self.get_commands().replace("'","\"").replace("\"","\\\"").replace("$","\$")
        s += "qsub \\\n"
        s += " \\\n".join(self.arguments)
        return s
        
    def __str__(self):
        return self.get_script()
        
    def run(self,filename='run.sh',dry=False,command="qsub",verbose=0):
        """
        run the command
        arguments:
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
        
        #check if there is stderr
        if self.stderr: raise Exception(self.stderr)
        
        #check if there is stdout
        if verbose: print(self.stdout)
        
        #get jobid
        self.jobid = str(self.stdout).split("\n")[0]
        if verbose: print("jobid:",self.jobid)

    
