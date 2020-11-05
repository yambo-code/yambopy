# Copyright (C) 2018 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
#
from __future__ import print_function, absolute_import
import os
import subprocess
import json
from copy import deepcopy

class Scheduler(object):
    """
    Generic scheduler class

    Initialize the class, by default we look for a config file in `~/.yambopy/config.json`

    This file has the following information:

        .. code-block:: javascript

            {
                "default": "oar <or> bash <or> pbs",
                "oar <or> bash <or> pbs":
                {
                    "mpirun": "<mpirun command>",
                    "modules": {"<module tag>":"<module name>"},
                    "pre_run": ["line1 to execute before the run",
                                "line2 to execute before the run"],
                    "pos_run": "file:<file_in_config_folder>",
                    "<tag>":"<value>"
                }
            }

    The "default" tag chooses the scheduler to use when no scheduler is
    specified by the user.

    "<tag>" and "<value>" in are additional variables that are handled
    differently for each scheduler. They are stored as arguments in the variable kwargs.
    """

    #TODO: Add a way to get the home directory
    _config_path     = os.path.expanduser("~") + "/.yambopy"
    _config_filename = "%s/config.json"%_config_path

    def __init__(self, name=None, nodes=None, cores=None, cpus_per_task=None, walltime="1:00:00", **kwargs ):
        self.name = name
        self.nodes = nodes
        self.cores = cores
        self.cpus_per_task = cpus_per_task
        self.walltime = walltime
        self.kwargs = kwargs

        self.pre_run = self.get_arg("pre_run",[])
        self.pos_run = self.get_arg("pos_run",[])

        self.modules_dict = self.get_arg("modules_dict",{})
        self.modules_list = []
        for mod in self.get_arg("modules_list",[]):
            self.add_module(mod)
        self.arguments = []
        self.commands  = []

        self.initialize()


    @classmethod    
    def factory(cls,scheduler=None,cores=None,nodes=None,walltime="1:00:00",**kwargs):
        """
        Initialize a scheduler instance.

        Default arguments:
        cores    - Number of cores to use
        nodes    - Number of nodes to use
        walltime - Walltime limit
        **kwargs - Additional tags specific for each scheduler
        """
        from .oar import Oar
        from .pbs import Pbs
        from .bash import Bash
        from .slurm import Slurm

        schedulers = { "oar":   Oar,
                       "pbs":   Pbs,
                       "slurm": Slurm,
                       "bash":  Bash }

        #load configurations file
        config = cls.load_config()

        #set from scheduler
        schedulername = scheduler
        if not scheduler: schedulername = config.get('default',None)

        if schedulername is None:
            raise ValueError('scheduler not specified and "default" not set in configuration file %s.'%self._config_filename)

        #load the configurations
        schedulerconfig = config.get(schedulername,{})

        #determine the scheduler type
        if "type" in schedulerconfig:
            schedulertype = schedulerconfig["type"]
        elif "oar"   in schedulername: schedulertype = "oar"
        elif "pbs"   in schedulername: schedulertype = "pbs"
        elif "bash"  in schedulername: schedulertype = "bash"
        elif "slurm" in schedulername: schedulertype = "slurm"
        else:
            raise ValueError("Could not determine the scheduler type. "
                             "Please specify it in the name of the scheduler or using the 'type' tag.")

        #sanity check
        if schedulername not in list(schedulers.keys()):
            raise ValueError("Scheduler name %s is invalid"%schedulername)

        #check type from outside class
        cls.schedulertype = schedulertype

        if "nodes" in schedulerconfig and nodes is None:
            nodes = int(schedulerconfig["nodes"])
            del schedulerconfig["nodes"]
        if "cores" in schedulerconfig and cores is None:
            cores = int(schedulerconfig["cores"])
            del schedulerconfig["cores"]

        walltime = schedulerconfig.pop("walltime",walltime)

        #add scheduler config arguments
        kwargs.update(schedulerconfig)

        #create an instance of the scheduler to use
        return schedulers[schedulertype](cores=cores,nodes=nodes,walltime=walltime,**kwargs)

    @staticmethod
    def load_config():
        """
        load the configuration file and check its sanity
        """
        try:
            with open(Scheduler._config_filename) as f:
                config = json.load(f)
        except IOError:
            config = {  "default":"bash",
                        "bash":
                        {
                            "mpirun": "mpirun",
                            "modules": "None",
                            "pre_run": [],
                            "pos_run": []
                        }
                     }

        #TODO: put sanity checks here
        schedulername = "bash"
        if 'default' in config:
            schedulername = config['default']

        return config

    def initialize(self):
        raise NotImplementedError("Initialize not implemented in this class!")

    def clean(self):
        """
        clean the command list
        """
        self.commands = []
        self.modules_list = []

    def get_vardict(self):
        """
        initialize the vardict from the default values and check in the arguments if
        a new assignment of the variable was specified.
        For example to assing the cores to ppn the user can set in the configurations file the
        "var_cores":"ppn" or to the arguments of the scheduler.factory() function as var_cores="ppn".
        """
        #load default vardict
        self.vardict = deepcopy(self._vardict)
        cores = self.get_arg('var_cores')
        if cores: self.vardict['cores'] = cores
        nodes = self.get_arg('var_nodes')
        if nodes: self.vardict['nodes'] = nodes
    
    def copy(self):
        """
        return a copy of this instance
        """
        return deepcopy(self)

    def write(self,filename):
        """
        write the bash file to a file in the disk
        """
        if not filename: filename = self.name
        with open(filename,"w") as f:
            f.write(str(self))

    def get_arg(self,argument,default=None):
        """
        get an argument from the list of optional variables kwargs

        return:
        None - if the variable does not exist
        value - if the variable exists return the value of the variable
        """
        if argument in self.kwargs:
            arg = self.kwargs[argument]
            if arg == "true":
                return True
            elif arg == "false":
                return False
            elif type(arg) in [str,str] and "file:" in arg:
                #load the text from a textfile
                filename = arg.split(":")[-1]
                f = open("%s/%s"%(self._config_path,filename),'r')
                return f.read().strip().split('\n')
                f.close()
            else:
                return self.kwargs[argument]
        else:
            return default

    def add_command(self,cmd):
        """
        add commands to be run by the scheduler
        """
        self.commands.append(cmd)

    def add_arguments(self,arguments):
        """
        add commands to be run by the scheduler
        """
        self.arguments.append(arguments)

    def add_mpirun_command(self,cmd):
        """
        add commands to be run by the scheduler
        """
        threads = 1
        if self.cores: threads*=self.cores
        if self.nodes: threads*=self.nodes
        mpirun = self.get_arg("mpirun","mpirun")
        np = self.get_arg("np","-np")
        self.add_command("%s %s %d %s"%(mpirun,np,threads,cmd))

    def add_module(self,mod):
        """
        Add module to be loaded by the scheduler
        
        Arguments:
            mod: if the mod string exists in the "modules_dict" dictionary get the key value
                 if it does not exist simply add the string assuming the user knows what he is doing
        """
        if self.modules_dict:
          if mod in self.modules_dict:
            self.modules_list.append(self.modules_dict[mod])
            return
        self.modules_list.append(mod)

    def run(self):
        """
        run the command in the bash
        """
        raise NotImplementedError('Run not implemented')

    def set_posrun(self,posrun):
        self.pos_run = posrun

    def set_prerun(self,prerun):
        self.pre_run = prerun
       
    @property
    def modulelist(self):
        return ["module load %s"%mod for mod in self.modules_list]

    def get_commands(self):
        """
        get commands to execute
        """
        lines = []; app = lines.append
        if self.pre_run:
            app("\n#pre_run")
            lines += self.pre_run
        if self.modules_list:
            app("\n#modules")
            lines += self.modulelist
        app("\n#commands")
        lines += self.commands
        if self.pos_run:
            app("\n#pos_run")
            lines += self.pos_run
        return "\n".join(lines)

if __name__ == "__main__":
    #after
    s = Scheduler(nodes=1,cores=4)
    s.add_command('abinit < run.sh')
    print(s)
