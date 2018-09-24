# Copyright (C) 2018 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
"""
This file contains classes to handle tasks.

The rules to implment anything here are:
  1. Minimalistic
  2. The state is stored in human readable files written in the disk to easily correct problems during the run
  3. The python code defines or is aware of where the main code (pw.x, ph.x, yambo, p2y, a2y, ypp) writes all the files
  4. Each task has a list of dependencies, only the last task is executed by the flow which will imply that 
     all the dependenecies need to be executed before. Only once the results of all the dependencies are obtained
     we execute the last task.
  5. Parallelism is possible when a task depends on a list of other tasks, meaning those tasks can be executed in parallel.
  6. The flow is the container of all the tasks and has the information about absolute paths and the overall scheduler.

The hierarchy is:

->YambopyFlow
  ->YambopyTask
    ->YamboTask
    ->P2yTask
    ->YppTask
    ->PwTask
"""

import os
import shutil
from schedulerpy import Scheduler
from qepy.pw import PwIn
from qepy import qepyenv
from yambopy import yambopyenv
from yambopy.tools.duck import isiter
from yambopy.io.inputfile import YamboIn

__all__ = [
    'YambopyFlow',
    'YambopyTask',
    'YamboTask',
    'P2yTask',
    'YppTask',
    'PwTask',
]

def write_fake(filename):
    with open(filename,'w') as f:
        f.write('')

def merge_two_dicts(x, y):
    """ taken from: https://stackoverflow.com/questions/38987/how-to-merge-two-dictionaries-in-a-single-expression"""
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z

class YambopyFlow():
    """
    Handle multiple tasks and their interdependencies
    Monitor the progress
    """
    def __init__(self,foldername,tasks):
        if not isiter(tasks): 
            self._tasks = [tasks]
        self._tasks = tasks 
        self.foldername = foldername

    @classmethod
    def from_tasks(cls,foldername,tasks):
        return cls(foldername,tasks)

    @property
    def dependencies(self):
        """Get a list of all the tasks"""
        dependencies = []
        for task in self._tasks:
            dependencies.append(task.dependencies)
        return dependencies
   
    @property
    def tasks(self):
        return self._tasks

    @property
    def ntasks(self):
        return len(self.tasks)

    def create(self):
        """Create a folder to run the flow"""
        if os.path.isdir(self.foldername):
            raise ValueError('A folder with name %s already exists. Please remove it of change the name of the flow')
        os.mkdir(self.foldername)

        #initialize each task
        for it,task in enumerate(self.tasks):
            #create folder   
            os.mkdir(os.path.join(self.foldername,'t%d'%it))
            #initialize each task
            task.initialize(os.path.join(self.foldername,'t%d'%it))

        #create a general run script
        lines = ['cd t%s; sh run.sh; cd ..'%n for n in range(self.ntasks)] 
        with open(os.path.join(self.foldername,'run.sh'),'w') as f:
            f.write('\n'.join(lines))

    def clean(self):
        shutil.rmtree(self.foldername)

    def get_status(self):
        lines = []; app = lines.append
        return "\n".join(lines)

    def __str__(self):
        lines = []; app = lines.append
        for it,task in enumerate(self.tasks):
            app(str(task))
        return "\n".join(lines)


class YambopyTask():
    """
    Generic task for yambopy

    Arguments:
        inputs: a list of input files
        commands: a list of commands to de executed with this executable
        executable: the name of the executable
        outputs: the output files produced when running the executable
        scheduler: the scheduler to run the task
    """
    def __init__(self,inputs,executable,scheduler,dependencies=None):
        if not isiter(inputs): inputs = [inputs]
        self.inputs = inputs
        self.executable = executable
        self.nlog = None
        self.nerr = None
        self.scheduler = scheduler
        if dependencies is None:
            self._dependencies = None
        elif not isiter(dependencies): 
            self._dependencies = [dependencies]
        else:
            ValueError('Unknown dependency type')
        self.outputs = None

    @property
    def log(self):
        if self.nlog: return 'run%d.log'%self.nlog
        self.nlog = 1
        return 'run.log' 

    @property
    def err(self):
        if self.nerr: return 'err%d.log'%self.nerr
        self.nerr = 1
        return 'err.log' 

    @property
    def dependencies(self):
        """ Return all the dependencies of this task"""
        if self._dependencies is None: return (self, None)
        dependencies = []
        for dependency in self._dependencies:
            dependencies.append(dependency.dependencies)
        return (self, dependencies)

    def get_instances_from_inputs(self,instance):
        """get all the instances of class from inputs"""
        return [inp for inp in self.inputs if isinstance(inp,instance)]

    def run(self):
        """
        Run this task using the specified scheduler
        """
        pass
#
# code specific tasks
#
class YamboTask(YambopyTask):
    """
    Routines specific to yambo task
    """
    yamboin = 'yambo.in'

    @classmethod
    def from_runlevel(cls,interface_task,runlevel,yamboin_dict={},
                      executable=yambopyenv.YAMBO,scheduler=yambopyenv.SCHEDULER,dependencies=None):
        """ run yambo with the runlevel string to generate the inputfile """
        instance = cls(inputs=interface_task,executable=executable,scheduler=scheduler,dependencies=dependencies)
        instance.runlevel = runlevel
        instance.yamboin_dict = yamboin_dict
        return instance

    @classmethod
    def from_folder(cls,folder):
        """
        create a task container from a pre-existing folder
        will use the YamboFolder and YamboFile class to detect inputs and outputs
        """
        raise NotImplementedError('TODO')
        #search input file
        #search output files
        return cls(inputs=yamboinput,executable=executable,scheduler=scheduler,dependencies=dependencies)

    @classmethod
    def from_inputfile(cls,yamboinput,executable=yambopyenv.YAMBO,scheduler=yambopyenv.SCHEDULER,dependencies=None):
        """Create a yambotask from an existing input file"""
        raise NotImplementedError('TODO')
        return cls(inputs=yamboinput,executable=execulable,scheduler=scheduler,dependencies=dependencies)

    def initialize(self,path):
        """ Initialize an yambo task """
        #get output from p2y task
        p2ys = self.get_instances_from_inputs(P2yTask)
        if len(p2ys) > 1:  raise ValueError('More than one P2yTask instance in input, cannot continue')
        if len(p2ys) == 0: raise ValueError('No P2yTask found in input, cannot continue')
        p2y = p2ys[0]
        src = os.path.abspath(p2y.output)
        dst = os.path.abspath(os.path.join(path,'SAVE'))
        os.symlink(src,dst)

        db1_path = os.path.join(dst,'ns.db1')
        if os.path.isfile(db1_path):
            #create inputfile
            yamboin = YamboIn.from_runlevel(self.runlevel,folder=path)
            merge_two_dicts(yamboin.arguments,self.yamboin_dict)
            yamboin.write(os.path.join(path,'yambo.in'))
        #else:
        #    raise FileNotFoundError('SAVE/ns.db1 not available in %s'%db1_path)

        #create running script
        abs_path = os.path.abspath(path)
        self._run = os.path.join(path,'run.sh')
        self.scheduler = Scheduler.factory(self.scheduler)
        self.scheduler.add_command('%s -F yambo.in -J run > %s 2> %s'%(self.executable,self.log,self.err))
        self.scheduler.write(self._run)

        self.path = path

class P2yTask(YambopyTask):
    """
    Run a P2Y calculation
    """
    @classmethod
    def from_nscf_task(cls,nscf_task,executable=yambopyenv.P2Y,scheduler=yambopyenv.SCHEDULER,setup=yambopyenv.YAMBO):
        """
        specify a nscf .save folder to start the calculation from
        """
        instance = cls(nscf_task,executable,scheduler,dependencies=nscf_task)
        instance.setup = setup
        return instance

    def initialize(self,path):
        #get output from nscf task
        nscfs = self.get_instances_from_inputs(PwTask)
        if len(nscfs) > 1: raise ValueError('More than one PwTask instance in input, cannot continue')
        if len(nscfs):
            nscf = nscfs[0]
            src = os.path.abspath(nscf.output)
            dst = os.path.abspath(os.path.join(path,"%s.save"%nscf.pwinput.prefix))
            os.symlink(src,dst)

        #create running script
        abs_path = os.path.abspath(path)
        self._run = os.path.join(path,'run.sh')
        self.scheduler = Scheduler.factory(self.scheduler)
        ac = self.scheduler.add_command
        ac('cd %s; %s > %s 2> %s; cd %s'%(dst,self.executable,self.log,self.err,abs_path))
        ac('mv %s/SAVE .'%(dst))
        if self.setup: ac('%s > %s 2> %s'%(self.setup,self.log,self.err))
        self.scheduler.write(self._run)

        self.path = path

    @property
    def output(self):
        return os.path.join(self.path,'SAVE')

class YppTask(YambopyTask):
    """
    Run a ypp calculation
    """
    @classmethod
    def from_runlevel(cls,runlevel,executable=yambopyenv.YPP,dependencies=None):
        raise NotImplementedError('TODO')
        #search input file
        #search output files
        return cls(inputs=yppinput,executable=executable,scheduler=scheduler,dependencies=dependencies)

    def initialize(self,path):
        write_fake(os.path.join(path,'ypp'))

class PwTask(YambopyTask):
    """
    Routines specific to quantum espresso task
    """
    @classmethod
    def from_input(cls,pwinputs,executable=qepyenv.PW,scheduler=yambopyenv.SCHEDULER,dependencies=None):
        if not isiter(pwinputs): pwinputs = [pwinputs]
        if not all([isinstance(pwi,(PwIn,cls)) for pwi in pwinputs]):
            raise ValueError('The input is not an instance of PwIn or PwTask but %s'%(pwinput))
        return cls(inputs=pwinputs,executable=executable,scheduler=scheduler,dependencies=dependencies)

    @property
    def pwinput(self):
        return self.get_instances_from_inputs(PwIn)[0]

    def initialize(self,path):
        """write inputs and get pseudopotential"""
        self.pwinput.write(os.path.join(path,'pw.in'))
        self.pwinput.get_pseudos(destpath=path)

        #in case there is another PwTask task in inputs link it
        self.link_pwtask(path)

        #create running script
        self._run = os.path.join(path,'run.sh')
        self.scheduler = Scheduler.factory(self.scheduler)
        self.scheduler.add_command('%s < pw.in > %s'%(self.executable,self.log))
        self.scheduler.write(self._run)

        self.path = path

    def link_pwtask(self,path):
        """Link pwtask from inputs with the current one"""
        pwinstances = self.get_instances_from_inputs(PwTask)
        if len(pwinstances) > 1: raise ValueError('More than one PwTask instance in input, cannot continue')
        if len(pwinstances):
            src = os.path.abspath(pwinstances[0].output)
            dst = os.path.abspath(os.path.join(path,"%s.save"%self.pwinput.prefix))
            os.symlink(src,dst)

    @property
    def output(self):
        """normally the output of a pw task is a .save folder"""
        return os.path.join(self.path,"%s.save"%self.pwinput.prefix)
