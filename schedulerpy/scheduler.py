# Copyright (C) 2018 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
#
from __future__ import print_function, absolute_import
import subprocess
import json
from abc import ABCMeta, abstractmethod
from textwrap import dedent
from copy import deepcopy
import os


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

    def __init__(self, name=None, nodes=None, cores=None, walltime="1:00:00", **kwargs ):
        self.name = name
        self.nodes = nodes
        self.cores = cores
        self.walltime = walltime
        self.kwargs = kwargs

        self.pre_run = self.get_arg("pre_run")
        self.pos_run = self.get_arg("pos_run")

        #modules = self.get_arg("modules")
        #if modules is None: self.modules = []
        #else:               self.modules = modules
        self.modules   = []
        self.arguments = []
        self.commands  = []

        self.initialize()


    @classmethod    
    def factory(cls,scheduler=None,cores=None,nodes=None,walltime="1:00:00",**kwargs):
        """
        Initialize a schduler instance.

        Default arguments:
        cores    - Number of cores to use
        nodes    - Number of nodes to use
        walltime - Walltime limit
        **kwargs - Additional tags specific for each scheduler
        """

        from . import oar
        from . import pbs
        from . import bash

        schedulers = { "oar":  oar.Oar,
                       "pbs":  pbs.Pbs,
                       "bash": bash.Bash }

        #load configurations file
        config = cls.load_config()
        if scheduler is None:
            schedulername = config['default']
        else:
            schedulername = scheduler

        #load the configurations
        if schedulername in config:
            schedulerconfig = config[schedulername]
        else:
            schedulerconfig = {}

        #determine the scheduler type
        if "type" in schedulerconfig:
            schedulertype = schedulerconfig["type"]
        elif "oar"  in schedulername: schedulertype = "oar"
        elif "pbs"  in schedulername: schedulertype = "pbs"
        elif "bash" in schedulername: schedulertype = "bash"
        else:
            raise ValueError("Could not determine the scheduler type. "
                             "Please specify it in the name of the scheduler or using the 'type' tag.")

        #sanity check
        if schedulername not in list(schedulers.keys()):
            raise ValueError("Scheduler name %s is invalid"%schedulername)

        if "nodes" in schedulerconfig and nodes is None:
            nodes = int(schedulerconfig["nodes"])
            del schedulerconfig["nodes"]
        if "cores" in schedulerconfig and cores is None:
            cores = int(schedulerconfig["cores"])
            del schedulerconfig["cores"]

        #add scheduler config arguments
        kwargs.update(schedulerconfig)

        #create an instance of the scheduler to use
        return schedulers[schedulertype](cores=cores,nodes=nodes,walltime=walltime,**kwargs)

    def load_config():
        """
        load the configuration file and check its sanity
        """
        try:
            f = open(Scheduler._config_filename)
            config = json.load(f)
            f.close()
        except IOError:
            config = {  "default":"bash",
                        "bash":
                        {
                            "mpirun": "mpirun",
                            "modules": "None",
                            "pre_run": ["echo 'running job...'"],
                            "pos_run": ["echo 'done!'"]
                        }
                     }

        #put sanity checks here
        schedulername = "bash"
        if 'default' in config:
            schedulername = config['default']
        #print json.dumps(config,indent=2)

        return config
    load_config = staticmethod(load_config)

    def initialize(self):
        raise NotImplementedError("Initialize not implemented in this class!")

    def clean(self):
        """
        clean the command list
        """
        self.commands = []
        self.modules  = []

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

    def write(self,filename):
        """
        write the bash file to a file in the disk
        """
        if not filename: filename = self.name
        f = open(filename,"w")
        f.write(str(self))
        f.close()

    def get_arg(self,argument):
        """
        get an argument from the list of optional variables kwargs

        return:
        None - if the variable does not exist
        value - if the variable existe return the value of the variable
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
            return None

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
        mpirun = self.get_arg("mpirun")
        if mpirun is None: mpirun = "mpirun"
        self.commands.append("%s %s"%(mpirun,cmd))

    def add_module(self,mod):
        """
        add module to be loaded by the scheduler
        """
        if "modules" in self.kwargs:
            config_modules = self.kwargs["modules"]
            if not config_modules == "None":
              if mod in config_modules:
                  self.modules.append(config_modules[mod])
              else:
                  raise ValueError("Module \"%s\" is now known. "
                                   "Add it to the config file in %s. "
                                   "Known modules are: %s"%(mod,self._config_filename,config_modules))
        else:
            raise ValueError("Option 'modules' not found in the specified type of scheduler. "
                             "Add the option to the config file in %s. "
                             "If you do not use modules, specify 'modules':'None' "%(self._config_filename))

    def run(self):
        """
        run the command in the bash
        """
        raise NotImplementedError('Run not implemented')

    def get_commands(self):
        """
        get commands to execute
        """
        s="\n"
        if self.pre_run: s += "#pre_run\n" + "\n".join(self.pre_run)+"\n\n"
        if self.modules: s += "#modules\n" + "\n".join(["module load %s"%mod for mod in self.modules])+"\n\n"
        s += "#commands\n"+ "\n".join(self.commands)+"\n\n"
        if self.pos_run: s += "#pos_run\n" + "\n".join(self.pos_run)+"\n\n"
        return s

if __name__ == "__main__":
    #after
    s = Scheduler(nodes=1,cores=4)
    s.add_command('abinit < run.sh')
    print(s)
