#
# Author: Henrique Pereira Coutada Miranda
# Tests for the yambopy library
# scheduler
#
from __future__ import print_function
import unittest
import sys
import os
import argparse
import subprocess
import filecmp
import json
from schedulerpy import *
from qepy import *
from textwrap import dedent

#name of the configuration file to be created when running the tests
_pre_run_pbs  = dedent( """
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
                        echo "threads : $OMP_NUM_THREADS"
                        """)
_pre_run_pbs_file = "pre_run_pbs.sh"
_test_config_file = "testconfig.json"
_test_config =  { "default":"bash",
                  "bash":
                    {
                        "mpirun": "mpirun",
                        "modules": "None",
                        "pos_run": ["echo \"done!\""],
                        "modules": {"abinit"  :"abinit/8.0",
                                    "espresso":"espresso/5.4.0"}
                    },
                  "oar":
                    {
                        "mpirun": "mpirun",
                        "pos_run": ["echo \"done!\""],
                        "idempotent": "true",
                        "modules": {"abinit"  :"abinit/8.0",
                                    "espresso":"espresso/5.4.0"}
                    },
                  "slurm":
                    {
                        "mpirun": "mpirun",
                        "pos_run": ["echo \"done!\""],
                        "modules": {"abinit"  :"abinit/8.0",
                                    "espresso":"espresso/5.4.0"}
                    },
                  "pbs":
                    {
                        "mpirun": "mpirun",
                        "pre_run": "file:%s"%_pre_run_pbs_file,
                        "pos_run": ["qstat -f $PBS_JOBID"],
                        "queue": "myqueue",
                        "group_list": "mygroup",
                        "mem": "2624",
                        "var_cores": "ppn",
                        "var_nodes": "select",
                        "modules": {"abinit"  :"abinit/8.0",
                                    "espresso":"espresso/5.4.0"}
                    }
                }
_test_file = "test.sh"

def header(text):
    return "\n%s\n%s\n%s\n"%("="*10,text,"="*10)

def init_config():
    """
    create basic configuration file
    """
    #write pbs pre_run example
    f = open(_pre_run_pbs_file,'w')
    f.write(_pre_run_pbs)
    f.close()
    #write test configuration file
    f = open(_test_config_file,'w')
    json.dump(_test_config,f)
    f.close()
    #scheduler
    Scheduler._config_path     = os.getcwd()
    Scheduler._config_filename = _test_config_file

def clean_config():
    os.remove(_pre_run_pbs_file)
    os.remove(_test_config_file)

class TestScheduler(unittest.TestCase):
    """
    This class tests the scheduler from yambnopy.

    1. Deploy a basic configuration file in ~/.yambopy
    compare the produced bash scripts to be run reference files
    """
    def test_default(self):
        """
        run the echo command on the scheduler using the default scheduler
        """
        s = Scheduler.factory(cores=1)
        s.add_module("abinit")
        s.add_command("echo 'hello'")
        print(header("default: %s"%s.__class__))
        print(s)
        s.write(_test_file)

        #remove files
        os.remove(_test_file)

    def test_print(self):
        """
        run the echo command on the different schedulers
        """
        #create basic configuration file
        init_config()

        #run using that configuration file
        for schedulername in ["oar","pbs","bash"]:
            print(header(schedulername))
            s = Scheduler.factory(scheduler=schedulername,cores=1,nodes=2)
            s.add_module("abinit")
            s.add_command("echo 'hello'")
            print(s)

        #remove files
        clean_config()

class TestSchedulerRun(unittest.TestCase):
    """
    This class tests the scheduler from yambnopy.

    2. Execute the echo command and check if the local configuration is working
    """
    def test_default(self):
        """
        run a echo hello world using the default scheduler on this machine
        """
        #run using that configuration file
        s = Scheduler.factory(cores=1)
        print(header("default: %s"%s.__class__))
        s.add_module("abinit")
        s.add_command("echo 'hello world'")
        s.run()

    def test_run(self):
        """
        run a echo hello world on the different possible configurations
        """
        init_config()

        for schedulername in ["oar","pbs","bash","slurm"]:
            print(header(schedulername))
            s = Scheduler.factory(scheduler=schedulername,cores=1)
            s._config = _test_config_file
            s.add_module("abinit")
            s.add_command("echo 'hello world'")
            s.run(dry=True)

        #remove files
        clean_config()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Test the yambopy script.')
    parser.add_argument('-i', '--input', action="store_true",
                        help='Generate the bash files and compare with the reference ones')
    parser.add_argument('-f', '--full',  action="store_true",
                        help='Generate the bash files, run them and compare the results')
    parser.add_argument('-c', '--clean',  action="store_true",
                        help='Clean all the data from a previous run')
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    # Count the number of errors
    nerrors = 0

    ul = unittest.TestLoader()
    tr = unittest.TextTestRunner(verbosity=2)

    # Run the test
    if args.input:
        suite = ul.loadTestsFromTestCase(TestScheduler)
        nerrors += not tr.run(suite).wasSuccessful()

    if args.full:
        suite = ul.loadTestsFromTestCase(TestSchedulerRun)
        nerrors += not tr.run(suite).wasSuccessful()

    #clean tests
    if args.clean:
        print("cleaning...")
        os.system('rm -rf scf')
        print("done!")

    sys.exit(nerrors)
