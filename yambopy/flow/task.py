# Copyright (C) 2018 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
"""
This file contains classes to handle tasks.

The rules to implement anything here are:
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
from yambopy.env import yambopyenv
from yambopy.tools.string import marquee
from yambopy.tools.duck import isiter, isstring
from yambopy.io.inputfile import YamboIn

__all__ = [
'YambopyFlow',
'YambopyTask',
'YamboTask',
'YamboChiTask',
'P2yTask',
'YppTask',
'AbinitTask',
'E2yTask',
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
            tasks = [tasks]
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
        return [(it,task) for it,task in enumerate(self.tasks) if task.status == "ready" and not task.launched]

    @property
    def alldone(self):
        return all([task.status == "done" for task in self.tasks])

    def create(self,agressive=False):
        """Create a folder to run the flow"""
        if agressive:
            if os.path.isdir(self.path): shutil.rmtree(self.path)

        if os.path.isdir(self.path):
            raise ValueError('A folder with name %s already exists.\n'
                             'Please remove it of change the name of the flow'%self.path)
        os.mkdir(self.path)

        #initialize each task
        for it,task in enumerate(self.tasks):
            self.initialize_task(it,verbose=False)

        self.initialized = True
        self.pickle()

    def dump_run(self):
        """Create a bash script to run the whole flow"""
        #create a general run script
        lines = ['cd t%s; sh run.sh; cd ..'%n for n in range(self.ntasks)] 
        with open(os.path.join(self.path,'run.sh'),'w') as f:
            f.write('\n'.join(lines))

    def run(self,maxexecs=1,sleep=5,dry=False,verbose=0):
        """Run all the tasks"""
        if not self.initialized: self.create()

        print(marquee("YambopyFlow.run"))
        while not self.alldone:
            #exeute maxexecs ready tasks
            for it,task in self.readytasks[:maxexecs]:
                print("%5s %10s  %s"%("t%d"%it,str(task.name),task.status))
                self.initialize_task(it,verbose=False)
                task.run(dry=dry)

            #wait some seconds
            time.sleep(sleep)
        
            #transform into a pickle
            self.pickle()

    def initialize_task(self,itask,verbose=True):
        task = self[itask]
        if not task.initialized:
            path = os.path.join(self.path,'t%d'%itask)
            if not os.path.isdir(path): os.mkdir(path)
            task.initialize(path)
            return 0
        if verbose: print('This task was already initialized')
        return 1

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
        if self.status != "ready": return
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
        if not isinstance(inputs,list): inputs = [inputs]
        self.inputs = inputs
        self.executable = executable
        self.nlog = 0
        self.nerr = 0
        self.initialized = initialized
        self.launched = False
        
        #code_injection
        self.code_dict = {}
        self.vars_dict = {}

        #scheduler
        self.scheduler_setup(scheduler)
 
        #dependencies
        if dependencies is None: self._dependencies = None
        elif not isiter(dependencies): dependencies = [dependencies]
        else: ValueError('Unknown dependency type')
        self._dependencies = dependencies

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

    def get_instances_from_inputs(self,instance,n=None):
        """get all the instances of class from inputs"""
        instances = [inp for inp in self.inputs if isinstance(inp,instance)]
        #if len(instances) == 0:
        #    inputs_str = ",".join([inp.__class__.__name__ for inp in self.inputs])
        #    raise RuntimeError('Did not find an instance of %s among the inputs: [%s]'%(instance.__name__,inputs_str))
        if n==1: return instances[0]
        return instances[:n]

    def run(self,dry=False,verbose=1):
        """
        Run this task using the specified scheduler
        """
        #initialize the task
        if not self.initialized: self.initialize()

        if self.launched: return 0

        #if initialized run it
        if self.initialized:
            self.scheduler.run(self._run,dry=dry)
        else:
            raise ValueError('could not initialize task')

        #set it as launched
        self.launched = True
        return 1

    def copy(self):
        import copy
        return copy.deepcopy(self)

    def to_string(self,mark='='):
        lines = []; app = lines.append
        app(marquee(self.name,mark=mark))
        app('initialized: {}'.format(self.initialized))
        app('status:      {}'.format(self.status))
        app('exitcode:    {}'.format(self.exitcode))
        if self.initialized:
            app('path: {}'.format(self.path))
        return "\n".join(lines)

    def __str__(self):
        return self.to_string()
#
# code specific tasks
#
class YamboTask(YambopyTask):
    """
    Routines specific to yambo task
    """
    yamboin = 'yambo.in'

    @classmethod
    def from_runlevel(cls,interface_task,runlevel,yamboin_dict={},yamboin_args=[],dependencies=None,**kwargs):
        """ Run yambo with the runlevel string to generate the inputfile """
        scheduler = kwargs.pop('scheduler',yambopyenv.SCHEDULER)
        executable = kwargs.pop('executable',yambopyenv.YAMBO)
        instance = cls(inputs=interface_task,executable=executable,
                       scheduler=scheduler,dependencies=dependencies)
        instance.runlevel = runlevel
        instance.yamboin_dict = yamboin_dict
        instance.yamboin_args = yamboin_args
        return instance

    @classmethod
    def from_folder(cls,folder,filename='run.in',dependencies=None,**kwargs):
        """
        create a task container from a pre-existing folder
        will use YamboFile class to read the inputs
        """
        scheduler = kwargs.pop('scheduler',yambopyenv.SCHEDULER)
        executable = kwargs.pop('executable',yambopyenv.YAMBO)
        #search input file
        yamboinput = YamboIn.from_file(filename=filename,folder=folder)
        instance = cls(inputs=yamboinput,executable=executable,
                       scheduler=scheduler,dependencies=dependencies)
        instance.initialized = True
        instance.yamboinput = yamboinput
        instance.path = folder
        return instance

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
        #get output from interface task
        if self.status != "ready": return

        #link interface tasks
        x2ytasks = self.get_instances_from_inputs((P2yTask,E2yTask))
        for x2ytask in x2ytasks:
            src = os.path.abspath(x2ytask.output)
            dst = os.path.abspath(os.path.join(path,'SAVE'))
            if not os.path.isdir(dst): os.symlink(src,dst)

        #link yambo tasks
        yambotasks = self.get_instances_from_inputs(YamboTask)
        run_path = os.path.join(path,'run')
        if not os.path.isdir(run_path): os.mkdir(run_path)
        for yambotask in yambotasks:
            for src in yambotask.output:
                dst = os.path.abspath(os.path.join(run_path,os.path.basename(src)))
                if not os.path.isdir(dst): os.symlink(src,dst)

        #set to initiailized
        self.initialized = True
        self.path = path

        #code injector
        code = self.get_code('initialize')
        code(self)

        #create inputfile
        db1_path = os.path.join(path,'SAVE','ns.db1')
        if not os.path.isfile(db1_path):
            raise FileNotFoundError('SAVE/ns.db1 not available in %s'%db1_path)

        if verbose: print("Creating inputfile in %s"%path)
        self.yamboinput = YamboIn.from_runlevel(self.runlevel,executable=self.executable,folder=path)
        self.yamboinput.set_fromdict(self.yamboin_dict)
        self.yamboinput.set_fromargs(self.yamboin_args)
        self.yamboinput.write(os.path.join(path,'run.in'))

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

class YamboChiTask(YamboTask):
    """
    Task containing a yambo chi calculation
    The purpose of this class is to paralelize the calculation of chi.
    This is done by overloading the run method with one similar to the
    YambopyFlow.
    The input files for the missing q are created in initialize
    One job per q-point is launched at a time.
    Once the database is ready the input file for that q-point is removed
    """
    def run(self,maxexecs=None,sleep=5,dry=False):
        #initialize the task
        if not self.initialized: self.initialize()

        if not self.initialized:
            raise ValueError('could not initialize task')

        #if initialized run it
        while not self.alldone:
            #execute maxexecs ready tasks
            for iq,(done,launched) in enumerate(zip(self.done_chi,self.launched_chi)):
                if done or launched: continue
                iq = iq + 1 #pad iq
                #create input for this q-point
                inpq = self.yamboinput.copy()
                inpq.set_q(iq)
                if dry: print(inpq)
                else:   inpq.write(os.path.join(self.path,'runchiq%d.in'%iq))
                #create submission script for this q-point
                this_scheduler = self.scheduler.copy()
                run = os.path.join(self.path,'runchiq%d.sh'%iq)
                cmd = '%s -F runchiq%d.in -J run -C runchiq%d > runchiq%d.log 2> runchiq%d.err'%(self.executable,iq,iq,iq,iq)
                this_scheduler.add_mpirun_command(cmd)
                #launch job
                print("%10s"%("q%d"%iq))
                this_scheduler.run(run,dry=dry)
                if maxexecs: 
                    maxexecs = maxexecs-1
                    if maxexecs == 0: return

                #wait some seconds
                time.sleep(sleep)

            #check for deadlocks
            all_done_or_launched = all([ done or launched for (done,launched) in zip(self.done_chi,self.launched_chi)])
            if all_done_or_launched and not self.alldone:
                print('waiting for the jobs to finish')
                time.sleep(25)


    @property
    def alldone(self):
        """Check if all the chi databases are present"""
        return all(self.done_chi)

    @property
    def exitcode(self):
        if not self.initialized: return None
        if self.alldone: return "success"
        return None

    @property
    def nqpoints(self):
        """get how many q points for chi need to be calculated"""
        if not self.initialized: return None
        for var in ['QpntsRXp','QpntsRXd','QpntsRXs']:
            qpts = self.yamboinput.variables.get(var,None)
            if qpts is not None: return qpts[0][1]

    @property
    def launched_chi(self):
        """
        A job is known to be launched when the input file is present in the folder
        We use run.in for the job with all the q points and runchiqn.in for each n q point
        """ 
        import re
        import glob
        runchis = glob.glob(os.path.join(self.path,'runchiq*.in'))
        # get a list of launched/not launched
        launched = [False]*self.nqpoints
        for runchi in runchis:
            qidx = int(re.findall('([0-9]+)',os.path.basename(runchi))[-1])-1
            launched[qidx] = True
        return launched

    @property
    def done_chi(self):
        """Get list of all the chi objects"""
        import glob
        chis = glob.glob(os.path.join(self.path,'run','ndb.pp*'))
        # get a list of done/not done chis
        done = [False]*self.nqpoints
        for chi in chis:
            if 'fragment' not in chi: continue
            qidx = int(chi.split('_')[-1])-1
            done[qidx] = True
        return done

    def __str__(self):
        lines = []; app = lines.append
        app(self.to_string())
        app('initialized: %d'%self.initialized)
        if self.initialized: 
            app('nqpoints: %d'%self.nqpoints)
            for iq,(done,launched) in enumerate(zip(self.done_chi,self.launched_chi)):
                app("%5s %5s %5s"%("q%d"%(iq+1),done,launched))
        return '\n'.join(lines)

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

class AbinitTask(YambopyTask):
    """
    Simple class to run the Abinit code in a similar way as QE
    For a more complete framework to run Abinit please use Abipy
    """
    @classmethod
    def from_input(cls,abinitinput,executable="abinit",dependencies=None,**kwargs):
        from abipy.abio.inputs import AbinitInput
        if not isinstance(abinitinput,list): abinitinput = [abinitinput]
        scheduler = kwargs.pop('scheduler',yambopyenv.SCHEDULER)
        return cls(inputs=abinitinput,executable=executable,
                   scheduler=scheduler,dependencies=dependencies)

    def initialize(self,path):
        #get output from previous tasks
        if self.status != "ready": return

        #create an Abipy AbinitTask to prepare the workdir
        from abipy.flowtk import AbinitTask as AbipyAbinitTask
        from abipy.abio.inputs import AbinitInput

        class FakeManager():
            def write_jobfile(self,filename):
                pass

        #append mandatory flags for yambo
        self.abinitinput.set_spell_check(False)
        self.abinitinput.set_vars(prtkbff=1,istwfk='*1',iomode=3)

        #use abipy to create files and folders
        abinit_task = AbipyAbinitTask.from_input(self.abinitinput,workdir=path)
        abinit_task.manager = FakeManager()
        abinit_task.build()
 
        #create links to files
        self.link_abinittask(path)
       
        #create submission script
        self._run = os.path.join(path,'run.sh')
        self.scheduler.add_mpirun_command('%s < run.files > %s'%(self.executable,self.log))
        self.scheduler.write(self._run)

        #set to initiailized
        self.initialized = True
        self.path = path

    def link_abinittask(self,path):
        """ Basic linker for Abinit tasks
        """
        #if getden is in the input file, link the DEN from outdir with the DEN from indir
        if 'getden' in self.abinitinput:
            abinittasks = self.get_instances_from_inputs(AbinitTask)
            if len(abinittasks) > 1: raise ValueError('More than one PwTask instance in input, cannot continue')
            if len(abinittasks):
                src = os.path.abspath(os.path.join(abinittasks[0].output,'out_DEN'))
                dst = os.path.abspath(os.path.join(path,"indata",'in_DEN'))
                if not os.path.isfile(src):
                    src += '.nc'
                    dst += '.nc'
                os.symlink(src,dst)

    @property
    def abinitinput(self):
        from abipy.abio.inputs import AbinitInput
        return self.get_instances_from_inputs(AbinitInput)[0]

    @property
    def output(self):
        """The output of an abinit task is contained in the outdata dir"""
        return os.path.join(self.path,"outdata")

class E2yTask(YambopyTask):
    """ Run the e2y driver to convert Abinit databases
    """
    @classmethod
    def from_wfk_file(cls,wfk_file):
        raise NotImplementedError("TODO")
        return cls(inputs=pwinputs,executable=executable,
                   scheduler=scheduler,dependencies=dependencies)

    @classmethod
    def from_folder(cls,folder,**kwargs):
        """
        Initialize a E2Y task from a folder.
        Useful to start from a SAVE folder calculation
        """
        scheduler = kwargs.pop('scheduler',yambopyenv.SCHEDULER)
        executable = kwargs.pop('executable',yambopyenv.E2Y)
        instance = cls(inputs=None,executable=executable,scheduler=scheduler,dependencies=None)
        instance.path = folder
        instance.initialized = True
        return instance

    @classmethod
    def from_nscf_task(cls,abinittask,**kwargs):
        setup = kwargs.pop('setup',yambopyenv.YAMBO)
        executable = kwargs.pop('executable',yambopyenv.E2Y)
        scheduler = kwargs.pop('scheduler',yambopyenv.SCHEDULER)
        instance = cls(inputs=abinittask,executable=executable,
                       scheduler=scheduler,dependencies=abinittask)
        instance.setup = setup
        return instance 

    def initialize(self,path):
        if self.status != "ready": return
        #get output from nscf task
        nscfs = self.get_instances_from_inputs(AbinitTask)
        if len(nscfs) > 1: raise ValueError('More than one AbinitTask instance in input, cannot continue')
        if len(nscfs):
            nscf = nscfs[0]
            src = os.path.abspath(os.path.join(nscf.output,'out_WFK.nc'))
            dst = os.path.abspath(os.path.join(path,'in_WFK.nc'))
            os.symlink(src,dst)

        #create running script
        abs_path = os.path.abspath(path)
        self._run = os.path.join(path,'run.sh')
        ac = self.scheduler.add_command
        ac('%s -F %s > %s 2> %s'%(self.executable,dst,self.log,self.err))
        if self.setup: ac('%s > %s 2> %s'%(self.setup,self.log,self.err))
        self.scheduler.write(self._run)

        #set to initiailized
        self.initialized = True
        self.path = path

    @property
    def output(self):
        return os.path.join(self.path,'SAVE')

class PwTask(YambopyTask):
    """
    Routines specific to quantum espresso task
    """
    @classmethod
    def from_input(cls,pwinputs,dependencies=None,**kwargs):
        scheduler = kwargs.pop('scheduler',yambopyenv.SCHEDULER)
        executable = kwargs.pop('executable',qepyenv.PW)
        paralelization = kwargs.pop('paralelization','')
        if not isiter(pwinputs): pwinputs = [pwinputs]
        if not all([isinstance(pwi,(PwIn,cls)) for pwi in pwinputs]):
            raise ValueError('The input is not an instance of PwIn or PwTask but %s'%(pwinput))

        instance = cls(inputs=pwinputs,executable=executable,
                   scheduler=scheduler,dependencies=dependencies)
        instance.paralelization = paralelization
        return instance

    @property
    def pwinput(self):
        return self.get_instances_from_inputs(PwIn)[0]

    def initialize(self,path):
        """write inputs and get pseudopotential"""
        #get output from interface task
        if self.status != "ready": return

        #code injector
        code = self.get_code('initialize')
        code(self)

        #set to initiailized
        self.initialized = True
        self.path = path

        #write input
        self.pwinput.write(os.path.join(path,'pw.in'))
        self.pwinput.get_pseudos(destpath=path)

        #in case there is another PwTask task in inputs link it
        self.link_pwtask(path)

        paralelization = self.paralelization if hasattr(self,'paralelization') else ""

        #create running script
        self._run = os.path.join(path,'run.sh')
        self.scheduler.add_mpirun_command('%s %s -inp pw.in > %s'%(self.executable,paralelization,self.log))
        self.scheduler.write(self._run)

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



