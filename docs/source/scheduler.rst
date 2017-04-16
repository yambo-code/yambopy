Schedulerpy
==========================
Schedulerpy is a simple set of python modules to run applications on different schedulers (PBS, OAR) and bash using the same python scripts.

Basic concept
--------------------------
We initialize the python scheduler with:

.. code-block:: python

    s = Scheduler.factory(nodes=1,cpus=4)
    s.add_module('yambo/4.0')
    s.add_command('yambo -F yambo.in')
    s.run()

The class will read the environment to run the command from a `config.json` file present in `~/.yambopy/config.json`. If the job is to be run through a scheduler (OAR or PBS), the script will submit a job to it with the indicated requirements. If the job is to be run locally using bash, then the code will do that.

The handling of the command is made using different python classes that define the different interfaces.
There is one class per interface, and the currently implemented interfaces are OAR, PBS and bash.
Different interfaces can easily be added. We can for example create a class to run jobs remotely through ssh.

The main goal is to create only one python script that says which code to execute and to be able to run it on different computers, schedulers and environments. For that we define a local configuration file that instructs schedulerpy how to run the jobs, be it through a scheduler or bash.

Configurations file
----------------------------
Currently available options are:

* bash - Execute the job in the bash
* oar  - Use the OAR scheduler
* pbs  - Use the PBS scheduler

Initialize the class, by default we look for a config file in `~/.yambopy/configure.json`

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

The "default" tag chooses the scheduler to use when no scheduler is specified by the user.

`"<tag>"` and `"<value>"` in are additional variables that are handled differently for each scheduler. They are stored as arguments in the variable kwargs.

`"modules"` is a dictionary that matches a tag to a local module. A certain module might
have different names across different clusters or environments.
Using this simple scheme the user can define one name that is translated differently in each platform 
to the name of the module.

Examples of configurations file
--------------------------------------------

1. Zenobe


2. Mac computer

.. code-block:: bash

   { "default":"bash"
   }

3. Gaia

.. code-block:: bash

   { "default":"oar",
     "bash" :
     {
      "modules": "None"
     },
     "oar" :
     {
      "mpirun": "mpirun",
      "modules": {"abinit"  :"abinit/8.0",
                  "espresso":"espresso/5.4.0-gcc",
                  "yambo":"yambo/master-intel"},
                  "pre_run": "file:pre_run_oar.sh",
                  "pos_run": ["echo 'done!'"]
     }
   }


