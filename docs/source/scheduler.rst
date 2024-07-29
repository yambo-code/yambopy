schedulerpy
==========================
Schedulerpy is a simple set of python modules to run applications on different schedulers (PBS, OAR) and BASH using the same python scripts.

Basic concept
--------------------------
The idea is to allow the user to create a python script describing a worflow that can be used in different environments
with minor modifications. For that, there is an abstract ``scheduler`` class that takes all the commands to be executed.
Once all the commands are specified, the scheduler class reads from a ``config.json`` file that should be
created in ``~/.yambopy/config.json``. In this file the user can choose which scheduler should be used by default and
other configurations.
We initialize the python scheduler with:

.. code-block:: python

    s = Scheduler.factory(nodes=1,cpus=4)
    s.add_module('yambo/4.0')
    s.add_command('yambo -F yambo.in')
    s.run()

If the job is to be run through a scheduler (OAR or PBS), the script will submit a job to it with the indicated requirements.
If the job is to be run locally using BASH, then the code use ``os.system`` from the standard python distribution to do that

The handling of the command is made using different python classes that define the different interfaces.
There is one class per interface, and the currently implemented interfaces are OAR, PBS and BASH.
Different interfaces can easily be added. We can for example create a class to run jobs remotely through ssh.

The main goal is to have only one python script that says which code to execute and to be able to
run it on different computers, schedulers and environments.
For that we specify the configuration in the local ``~/.yambopy/config.json`` file that instructs 
schedulerpy how to run the jobs, be it through a scheduler or BASH.

Configurations file
----------------------------
Currently available options for schedulersare:

* ``bash`` - Execute the job in the bash
* ``oar``  - Use the OAR scheduler
* ``pbs``  - Use the PBS scheduler

Initialize the class, by default we look for a config file in ``~/.yambopy/config.json``

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

`"<tag>"` and `"<value>"` in are additional variables that are handled differently for each scheduler.
They are stored as arguments in the variable kwargs.

`"modules"` is a dictionary that matches a tag to a local module. A certain module might
have different names across different clusters or environments.
Using this simple scheme the user can define one name that is translated differently in each platform 
to the name of the module.

Examples of configurations file
--------------------------------------------
Here are some example configuration files that we are currently using.

1. Local computer (Linux, Mac)
~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   { "default":"bash" }


2. Zenobe
~~~~~~~~~~~~

Zenobe is a cluster in Belgium part of Céci
http://www.ceci-hpc.be/clusters.html#zenobe

.. code-block:: bash

    {
        "default":"pbs",
        "bash":
        {
            "mpirun": "mpirun",
            "pre_run": ["echo 'running job...'"],
            "pos_run": ["echo 'done!'"]
        },
        "pbs":
        {
            "modules": {"yambo":"yambo/git-slepc-intel"},
            "mpirun": "mpirun",
            "mem": 2600,
            "var_nodes":"select",
            "var_cores":"ncpus",
            "group_list": "<group_list_name>",
            "pre_run": "file:pre_run_pbs.sh",
            "pos_run": ["echo 'done!'"]
        }
    }


3. Gaia
~~~~~~~~~~~~~~~~~~~

Gaia is a cluster in Luxembourg part of the University of Luxembourg
https://hpc.uni.lu/systems/gaia/

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


Also add the ``pre_run_oar.sh`` file in  ``~/.yambopy/`` with the following

.. code-block:: bash

    if [ -f  /etc/profile ]; then
        . /etc/profile;
    fi;
    module use /work/projects/tss-physics/modules/

