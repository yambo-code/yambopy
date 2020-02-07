# Choose the scheduler of your cluster

from schedulerpy import Scheduler
import sys
import argparse
import os

parser = argparse.ArgumentParser(description='Test schedulerpy.')

parser.add_argument('-b','--bash', action="store_true", help='Run bash')
parser.add_argument('-s','--slurm', action="store_true", help='Run slurm')
parser.add_argument('-o','--oar', action="store_true", help='Run Oar')
parser.add_argument('-p','--pbs', action="store_true", help='Run pbs')

args = parser.parse_args()

# set current directory
cwd = os.getcwd()

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

if args.bash:
   print('check the run.sh file')
   sch = Scheduler.factory(scheduler="bash",ntasks=1,walltime="01:00:00")
   sch.add_module('quantumespresso/6.1')
   sch.add_command('echo $PWD')
   sch.add_command('echo $? > __yambopystatus__')
   sch.run(filename=cwd+'/run.sh')#  '/Users/alejandro/Software/yambopy/tutorial/bn/run.sh')

if args.slurm:
   sch = Scheduler.factory(scheduler="slurm",ntasks=2,walltime="01:00:00")
   sch.add_module('quantumespresso/6.1')
   sch.add_command('echo $PWD')
   sch.add_command('echo $? > __yambopystatus__')
   sch.run(filename=cwd+'/run.sh',dry='dry')

if args.oar:
   sch = Scheduler.factory(scheduler="oar",ntasks=2,walltime="01:00:00")
   sch.add_module('quantumespresso/6.1')
   sch.add_command('echo $PWD')
   sch.add_command('echo $? > __yambopystatus__')
   sch.run(dry=True)

if args.pbs:
   sch = Scheduler.factory(scheduler="pbs",ntasks=2,walltime="01:00:00")
   sch.add_module('quantumespresso/6.1')
   sch.add_command('echo $PWD')
   sch.add_command('echo $? > __yambopystatus__')
   sch.run(filename=cwd+'/run.sh',dry='dry')
