from schedulerpy import *
s = Scheduler.factory(nodes=1,cores=5)
s.add_module('abinit')
s.add_command('echo \'hello\'')
s.add_mpirun_command('echo \'hello\'')
print s
