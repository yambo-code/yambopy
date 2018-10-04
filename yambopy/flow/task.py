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
import time
from schedulerpy import Scheduler
from qepy.pw import PwIn
from qepy.ph import PhIn
from qepy.dynmat import DynmatIn
from qepy import qepyenv
from yambopy import yambopyenv
from yambopy.tools.string import marquee
from yambopy.tools.duck import isiter, isstring
from yambopy.io.inputfile import YamboIn

__all__ = [
'YambopyFlow',
'YambopyTask',
'YamboTask',
'P2yTask',
'YppTask',
'PwTask',
'PhTask',
'DynmatTask',
]

def write_fake(filename):
    with open(filename,'w') as f:
        f.write('')

class YambopyFlow(object):
    """
    Handle multiple tasks and their interdependencies
    Monitor the progress
    """
    _picklename = "__yambopyflow__"

    def __init__(self,path,tasks):
        if not isiter(tasks): 
            self._tasks = [tasks]
        self._tasks = tasks 
        self.path = path
        self.initialized = False

    @classmethod
    def from_tasks(cls,path,tasks):
        return cls(path,tasks)

    @classmethod
    def from_pickle(cls,pickle_filename):
        import pickle
        with open(pickle_filename,'rb') as f:
            pickle = pickle.load(f)
        return pickle

    @classmethod
    def from_folder(cls,flow_folder):
        pickle_filename = os.path.join(flow_folder,cls._picklename)
        if not os.path.isfile(pickle_filename):
            raise FileNotFoundError('pickle not found in %s'%pickle_filename)
        return cls.from_pickle(pickle_filename)

    @property
    def dependencies(self):
        """Get a list of all the tasks"""
        dependencies = []
        for task in self._tasks:
            dependencies.append(task.dependencies)
        return dependencies

    @property  
    def name(self):
        return self.__class__.__name__
 
    @property
    def tasks(self):
        return self._tasks

    @property
    def ntasks(self):
        return len(self.tasks)

    @property
    def readytasks(self):
        return [(it,task) for it,task in enumerate(self.tasks) if task.status == "ready"]

    @property
    def alldone(self):
        return all([task.status == "done" for task in self.tasks])

    def create(self):
        """Create a folder to run the flow"""
        if os.path.isdir(self.path):
            raise ValueError('A folder with name %s already exists.\n'
                             'Please remove it of change the name of the flow'%self.path)
        os.mkdir(self.path)

        #initialize each task
        for it,task in enumerate(self.tasks):
            #get path
            path = os.path.join(self.path,'t%d'%it)
            #create folder   
            os.mkdir(path)

            #initialize each task
            if not task.initialized: task.initialize(path)

        self.initialized = True
        self.pickle()

    def dump_run(self):
        """Create a bash script to run the whole flow"""
        #create a general run script
        lines = ['cd t%s; sh run.sh; cd ..'%n for n in range(self.ntasks)] 
        with open(os.path.join(self.path,'run.sh'),'w') as f:
            f.write('\n'.join(lines))

    def run(self,maxexecs=1,sleep=1,verbose=0):
        """Run all the tasks"""
        if not self.initialized: self.create()

        while not self.alldone:
            #exeute maxexecs ready tasks
            for it,task in self.readytasks[:maxexecs]:
                if not task.initialized: 
                    path = os.path.join(self.path,'t%d'%it)
                    task.initialize(path)
                task.run()

            #wait some seconds
            time.sleep(sleep)
 
    def pickle(self):
        """store the flow in a pickle"""
        import pickle
        pickle_filename = os.path.join(self.path,self._picklename)
        with open(pickle_filename,'wb') as f:
            pickle.dump(self,f)

    def clean(self):
        shutil.rmtree(self.path)

    def get_status(self):
        lines = []; app = lines.append
        return "\n".join(lines)

    def __getitem__(self,idx):
        return self.tasks[idx]

    def __str__(self):
        lines = []; app = lines.append
        app(marquee(self.name))
        for it,task in enumerate(self.tasks):
            app("%5s %10s  %s"%("t%d"%it,str(task.name),task.status))
        return "\n".join(lines)

def task_init(initialize):
    def new_initialize(self,path):
        initialize(self,path)
        self.path = path
        self.initialized = True
    return new_initialize
    
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
    _yambopystatus = "__yambopystatus__"

    def __init__(self,inputs,executable,scheduler,dependencies=None,initialized=False):
        if not isiter(inputs): inputs = [inputs]
        self.inputs = inputs
        self.executable = executable
        self.nlog = 0
        self.nerr = 0
        self.initialized = initialized
        
        #code_injection
        self.code_dict = {}
        self.vars_dict = {}

        #scheduler
        self.scheduler_setup(scheduler)
 
        #dependencies
        if dependencies is None: self._dependencies = None
        elif not isiter(dependencies): self._dependencies = [dependencies]
        else: ValueError('Unknown dependency type')

    def scheduler_setup(self,scheduler):
        """Setup the scheduler for this task"""
        if isstring(scheduler): self.scheduler = Scheduler.factory(scheduler)
        elif isinstance(scheduler,Scheduler): self.scheduler = scheduler.copy()
        else: raise ValueError('Invalid scheduler')

        #add a check for exitcode
        new_posrun = ['echo $? > %s'%self._yambopystatus] + self.scheduler.pos_run
        self.scheduler.set_posrun(new_posrun)

    def get_vars(self,key):
        return self.vars_dict[key]

    def set_vars(self,key,value):
        self.vars_dict[key] = value 

    def get_code(self,key):
        """get a funnction to execute"""
        def null_func(self): pass
        return self.code_dict.get(key,null_func)
 
    def set_code(self,key,code):
        """set a function to execute"""
        #TODO check consistency of the function signature
        self.code_dict[key] = code

    @property
    def ready(self):
        if self._dependencies == None: return True
        if all([dep.status == "done" for dep in self._dependencies]): return True
        return False

    @property
    def status(self):
        """
        status:
                        meaning
            - waiting - waiting for dependencies
            - ready   - ready to run
            - failed  - run and exitcode 1
            - done    - run and exitcode 0
        """
        if self.ready and self.exitcode == None: return "ready"
        if   self.exitcode == "success": return "done"
        elif self.exitcode == None: return "waiting"
        else: return "waiting"

    @property
    def name(self):
        return self.__class__.__name__

    @property
    def exitcode(self):
        """Check the exit code of the application"""
        if not self.initialized: return None
        exitcode_file = os.path.join(self.path,self._yambopystatus)
        if not os.path.isfile(exitcode_file):  return None
        with open(exitcode_file,'r') as f:
            exitcode = 'failure' if int(f.read()) else 'success'
        return exitcode

    @property
    def done(self):
        return self.exitcode == 0

    @property
    def log(self):
        if self.nlog: log = 'run%d.log'%self.nlog
        else:         log = 'run.log'
        self.nlog += 1
        return log

    @property
    def err(self):
        if self.nerr: err = 'err%d.log'%self.nerr
        else:         err = 'err.log'
        self.nerr += 1
        return err

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

    def run(self,dry=False,verbose=1):
        """
        Run this task using the specified scheduler
        """
        #initialize the task
        if not self.initialized: self.initialize()

        #if initialized run it
        if self.initialized:
            os.system('cd %s; sh run.sh'%self.path)
            #self.scheduler.run(dry=dry)
        else:
            raise ValueError('could not initialize task')

    def copy(self):
        import copy
        return copy.deepcopy(self)

    def __str__(self):
        lines = []; app = lines.append
        app(marquee(self.name))
        app('initialized: {}'.format(self.initialized))
        app('status:      {}'.format(self.status))
        app('exitcode:    {}'.format(self.exitcode))
        if self.initialized:
            app('path: {}'.format(self.path))
        return "\n".join(lines)
#
# code specific tasks
#
class YamboTask(YambopyTask):
    """
    Routines specific to yambo task
    """
    yamboin = 'yambo.in'

    @classmethod
    def from_runlevel(cls,interface_task,runlevel,yamboin_dict={},dependencies=None,**kwargs):
        """ Run yambo with the runlevel string to generate the inputfile """
        scheduler = kwargs.pop('scheduler',yambopyenv.SCHEDULER)
        executable = kwargs.pop('executable',yambopyenv.YAMBO)
        instance = cls(inputs=interface_task,executable=executable,
                       scheduler=scheduler,dependencies=dependencies)
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
        return cls(inputs=yamboinput,executable=executable,
                   scheduler=scheduler,dependencies=dependencies)

    @classmethod
    def from_inputfile(cls,yamboinput,dependencies=None,**kwargs):
        """Create a yambotask from an existing input file"""
        raise NotImplementedError('TODO')
        scheduler = kwargs.pop('scheduler',yambopyenv.SCHEDULER)
        executable = kwargs.pop('executable',yambopyenv.YAMBO)
        return cls(inputs=yamboinput,executable=execulable,
                   scheduler=scheduler,dependencies=dependencies)

    def initialize(self,path,verbose=0):
        """ Initialize an yambo task """
        #get output from p2y task
        if self.status != "ready": return

        #link p2y tasks
        p2ytasks = self.get_instances_from_inputs(P2yTask)
        if len(p2ytasks) == 1:
            src = os.path.abspath(p2ytasks[0].output)
            dst = os.path.abspath(os.path.join(path,'SAVE'))
            os.symlink(src,dst)

        #link yambo tasks
        yambotasks = self.get_instances_from_inputs(YamboTask)
        if len(yambotasks) == 1:
            run_path = os.path.join(path,'run')
            os.mkdir(run_path)
            for src in yambotasks[0].output:
                dst = os.path.abspath(os.path.join(run_path,os.path.basename(src)))
                os.symlink(src,dst)

        #set to initiailized
        self.initialized = True
        self.path = path

        #code injector
        code = self.get_code('initialize')
        code(self)

        #create inputfile
        db1_path = os.path.join(path,'SAVE','ns.db1')
        if os.path.isfile(db1_path):
            if verbose: print("Creating inputfile in %s"%path)
            yamboin = YamboIn.from_runlevel(self.runlevel,executable=self.executable,folder=path)
            yamboin.set_fromdict(self.yamboin_dict)
            yamboin.write(os.path.join(path,'run.in'))
        else:
            raise FileNotFoundError('SAVE/ns.db1 not available in %s'%db1_path)

        #create running script
        abs_path = os.path.abspath(path)
        self._run = os.path.join(path,'run.sh')
        self.scheduler.add_mpirun_command('%s -F run.in -J run > %s 2> %s'%(self.executable,self.log,self.err))
        self.scheduler.write(self._run)


    @property
    def output(self):
        """ 
        There are many different outputs from an Yambo task depending on the runlevel
        In general the different database outputs will be insite the run folder
        """
        outputs = []
        run_path = os.path.join(self.path,'run')
        for filename in os.listdir(run_path):
            outputfile = os.path.abspath(os.path.join(run_path,filename))
            outputs.append(outputfile)
        return outputs

class P2yTask(YambopyTask):
    """
    Run a P2Y calculation
    """
    @classmethod
    def from_nscf_task(cls,nscf_task,**kwargs):
        """
        specify a nscf .save folder to start the calculation from
        """
        scheduler = kwargs.pop('scheduler',yambopyenv.SCHEDULER)
        setup = kwargs.pop('setup',yambopyenv.YAMBO)
        executable = kwargs.pop('executable',yambopyenv.P2Y)
        instance = cls(nscf_task,executable,scheduler,dependencies=nscf_task)
        instance.setup = setup
        return instance

    @classmethod
    def from_folder(cls,folder,**kwargs):
        """
        Initialize a P2Y task from a folder.
        Useful to start from a SAVE folder calculation
        """
        scheduler = kwargs.pop('scheduler',yambopyenv.SCHEDULER)
        executable = kwargs.pop('executable',yambopyenv.P2Y)
        instance = cls(inputs=None,executable=executable,scheduler=scheduler,dependencies=None)
        instance.path = folder
        instance.initialized = True
        return instance

    @task_init
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
        ac = self.scheduler.add_command
        ac('cd %s; %s > %s 2> %s; cd %s'%(dst,self.executable,self.log,self.err,abs_path))
        ac('mv %s/SAVE .'%(dst))
        if self.setup: ac('%s > %s 2> %s'%(self.setup,self.log,self.err))
        self.scheduler.write(self._run)

        self.path = path
        self.initialized = True

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
        return cls(inputs=yppinput,executable=executable,
                   scheduler=scheduler,dependencies=dependencies)

    def initialize(self,path):
        write_fake(os.path.join(path,'ypp'))

class PwTask(YambopyTask):
    """
    Routines specific to quantum espresso task
    """
    @classmethod
    def from_input(cls,pwinputs,dependencies=None,**kwargs):
        scheduler = kwargs.pop('scheduler',yambopyenv.SCHEDULER)
        executable = kwargs.pop('executable',qepyenv.PW)
        if not isiter(pwinputs): pwinputs = [pwinputs]
        if not all([isinstance(pwi,(PwIn,cls)) for pwi in pwinputs]):
            raise ValueError('The input is not an instance of PwIn or PwTask but %s'%(pwinput))

        return cls(inputs=pwinputs,executable=executable,
                   scheduler=scheduler,dependencies=dependencies)

    @property
    def pwinput(self):
        return self.get_instances_from_inputs(PwIn)[0]

    @task_init
    def initialize(self,path):
        """write inputs and get pseudopotential"""
        self.pwinput.write(os.path.join(path,'pw.in'))
        self.pwinput.get_pseudos(destpath=path)

        #in case there is another PwTask task in inputs link it
        self.link_pwtask(path)

        #create running script
        self._run = os.path.join(path,'run.sh')
        self.scheduler.add_mpirun_command('%s -inp pw.in > %s'%(self.executable,self.log))
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


class PhTask(YambopyTask):
    """
    Routines specific to quantum espresso task
    """
    @classmethod
    def from_scf_task(cls,pwinputs,dependencies=None,**kwargs):
        scheduler = kwargs.pop('scheduler',yambopyenv.SCHEDULER)
        executable = kwargs.pop('executable',qepyenv.PH)
        if not isiter(pwinputs): pwinputs = [pwinputs]
        if not any([isinstance(pwi,(PhIn,cls)) for pwi in pwinputs]):
            raise ValueError('The input is not an instance of PwTask but %s'%(pwinputs))
        instance = cls(inputs=pwinputs,executable=executable,
                   scheduler=scheduler,dependencies=dependencies)
        #set PhIn instance prefix
        instance.phinput.prefix = instance.pwinput.prefix
        instance.phinput.fildyn = "%s.dyn"%instance.pwinput.prefix
        return instance

    @property
    def phinput(self):
        return self.get_instances_from_inputs(PhIn)[0]

    @property
    def pwtask(self):
        return self.get_instances_from_inputs(PwTask)[0]

    @property
    def pwinput(self):
        return self.pwtask.pwinput

    @task_init
    def initialize(self,path):
        """write inputs"""
        self.phinput.write(os.path.join(path,'ph.in'))
        self.pwinput.get_pseudos(destpath=path)

        #in case there is another PwTask task in inputs link it
        self.link_pwtask(path)

        #create running script
        self._run = os.path.join(path,'run.sh')
        self.scheduler.add_mpirun_command('%s -inp ph.in > %s'%(self.executable,self.log))
        self.scheduler.write(self._run)

        self.path = path

    def link_pwtask(self,path):
        """Link pwtask from inputs with the current one"""
        pwinstances = self.get_instances_from_inputs(PwTask)
        if len(pwinstances) > 1: raise ValueError('More than one PwTask instance in input, cannot continue')
        if len(pwinstances):
            src = os.path.abspath(pwinstances[0].output)
            dst = os.path.abspath(os.path.join(path,"%s.save"%self.phinput.prefix))
            os.symlink(src,dst)

    @property
    def output(self):
        """normally the output of a ph task is a set of dyn files"""
        outputs = []
        for filename in os.listdir(self.path):
            if "dyn" in filename:
                outputs.append(os.path.abspath(os.path.join(self.path,filename)))
        return outputs 

class DynmatTask(YambopyTask):
    """
    Routines specific to quantum espresso task
    """
    @classmethod
    def from_phonon_task(cls,phtask,dependencies=None,**kwargs):
        scheduler = kwargs.pop('scheduler',yambopyenv.SCHEDULER)
        executable = kwargs.pop('executable',qepyenv.DYNMAT)
        if not isinstance(phtask,PhTask):
            raise ValueError('The input is not an instance of PhTask but %s'%(pwinputs))
        return cls(inputs=phtask,executable=executable,
                   scheduler=scheduler,dependencies=dependencies)

    @property
    def phtask(self):
        return self.get_instances_from_inputs(PhTask)[0]

    @property
    def phinput(self):
        return self.phtask.pwinput

    def initialize(self,path):
        """write inputs and get pseudopotential"""
        if self.status != "ready": return

        #crate dynmat input file
        dynmat = DynmatIn.from_prefix(self.phinput.prefix) 
        dynmat.write(os.path.join(path,'dynmat.in'))

        #copy dyn files in this foldeer
        phtasks = self.get_instances_from_inputs(PhTask)
        for src in phtasks[0].output:
            dst = os.path.abspath(os.path.join(path,os.path.basename(src)))
            shutil.copy(src,dst)

        #create running script
        self._run = os.path.join(path,'run.sh')
        #dynamt is always serial
        self.scheduler.add_command('%s < dynmat.in > %s'%(self.executable,self.log))
        self.scheduler.write(self._run)

        self.path = path
        self.initialized = True
    
    @property
    def output(self):
        """normally the output of a pw task is a .save folder"""
        return os.path.join(self.path,"%s.save"%self.phinput.prefix)



