# Copyright (C) 2018 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
"""
This file contains classes and methods to generate input files

The FiniteDifferencesPhononFlow
"""

from qepy.pw import PwIn
from qepy.ph import PhIn
from qepy.matdyn import Matdyn
from qepy import qepyenv
from yambopy.io.inputfile import YamboIn
from yambopy.tools.duck import isiter
from yambopy.flow import PwTask, PhTask, P2yTask, YamboTask, DynmatTask, YambopyFlow

__all__ = [
"FiniteDifferencesPhononFlow",
"KpointsConvergenceFlow",
"PwNscfYamboIPChiTasks",
"YamboQPBSETasks",
"PwNscfTasks",
"PhPhononTasks"
]

CONV_THR = qepyenv.CONV_THR

class FiniteDifferencesPhononFlow():
    """
    This class takes as an input one structure and a phonon calculation.
    It produces a flow with the QE input files displaced along the phonon modes

    Author: Henrique Miranda
    """
    def __init__(self,structure,phonon_modes):
        self.structure = structure
        if not isinstance(phonon_modes,Matdyn):
            raise ValueError('phonon_modes must be an instance of Matdyn')
        self.phonon_modes = phonon_modes

    def get_tasks(self,path,kpoints,ecut,nscf_bands,nscf_kpoints=None,
                  imodes_list=None,displacements=[0.01,-0.01],iqpoint=0,**kwargs):
        """
        Create a flow with all the tasks to perform the calculation
        
        Arguments:
            generator: a builder function that takes structure, kpoints, ecut, 
                       nscf_bands, nscf_kpoints and kwargs as argument and returns the
                       tasks to be performed at each displacement
        """

        if imodes_list is None: imodes_list = list(range(self.phonon_modes.nmodes))
        if not isiter(displacements): displacements = [displacements]
        generator = kwargs.pop("generator",PwNscfYamboIPChiTasks)
             
        #create qe input from structure
        pwin = PwIn.from_structure_dict(self.structure,kpoints=kpoints,ecut=ecut)

        tasks = []
        #
        # Undisplaced structure
        #
        tasks.extend(generator(pwin.get_structure(),kpoints,ecut,nscf_bands,nscf_kpoints,**kwargs))

        #
        # Displaced structures
        #
        reference = [dict(imode=None,displacement=0)]

        for imode in imodes_list:
            #get phonon mode
            cart_mode = self.phonon_modes.modes[iqpoint,imode]
            #iterate over displacements
            for displacement in displacements:
                #displace structure
                input_mock = pwin.get_displaced(cart_mode, displacement=displacement)
                displaced_structure = input_mock.get_structure()

                #generate tasks
                new_tasks = generator(displaced_structure,kpoints,ecut,
                                      nscf_bands,nscf_kpoints,**kwargs)
                #save task
                tasks.extend(new_tasks)
                reference.append([dict(imode=None,displacement=displacement)])

        #store a reference for each of the tasks
        self.reference = reference
        return tasks

    
    def get_flow(self,path,kpoints,ecut,nscf_bands,nscf_kpoints=None,imodes_list=None,**kwargs):

        tasks = self.get_tasks(path=path,kpoints=kpoints,ecut=ecut,nscf_bands=nscf_bands,
                               nscf_kpoints=nscf_kpoints,imodes_list=imodes_list,**kwargs)
       
        #put all the tasks in a flow
        self.yambo_flow = YambopyFlow.from_tasks(path,tasks)
        return self.yambo_flow

    def get_dchi(self):
        """Collect all tasks producing eps files"""
        raise NotImplementedError('TODO')
        #find all the optical tasks
        #calculate finite differences
        pass

    def plot_ax(self,what):
        """ Collect and plot the results
        
        Arguments:
            what: what to plot can be eps for optical response
        """
        raise NotImplementedError('TODO')
        #plot finite differences
        pass        

    @property
    def path(self):
        return self.yambo_flow.path

class KpointsConvergenceFlow():
    """
    This class takes as an input one structure.
    It produces a flow with the QE input files with different number of kpoints
    and then runs a task staring from the QE input

    Author: Henrique Miranda
    """
    def __init__(self,structure):
        self.structure = structure
    
    def get_tasks(self,scf_kpoints,nscf_kpoints_list,ecut,nscf_bands,**kwargs):
        """
        Create a flow with all the tasks to perform the calculation
        
        Arguments:
            generator: a builder function that takes structure, kpoints, ecut, 
                       nscf_bands, nscf_kpoints and kwargs as argument and returns the
                       tasks to be performed at each displacement
        """
        generator = kwargs.pop("generator",PwNscfYamboIPChiTasks)

        #create a QE scf task and run
        qe_input = PwIn.from_structure_dict(self.structure,kpoints=scf_kpoints,ecut=ecut)
        qe_scf_task = PwTask.from_input(qe_input)
        tasks = [qe_scf_task]

        for nscf_kpoints in nscf_kpoints_list:

            #generate tasks
            new_tasks = generator(structure=self.structure,kpoints=scf_kpoints,
                                  nscf_kpoints=nscf_kpoints,
                                  ecut=ecut,nscf_bands=nscf_bands,**kwargs)
            tasks.extend(new_tasks)

        return tasks

    def get_flow(self,path,scf_kpoints,nscf_kpoints_list,ecut,nscf_bands,**kwargs):
        tasks = self.get_tasks(scf_kpoints=scf_kpoints,ecut=ecut,nscf_bands=nscf_bands,
                               nscf_kpoints_list=nscf_kpoints_list,**kwargs)
       
        #put all the tasks in a flow
        self.yambo_flow = YambopyFlow.from_tasks(path,tasks,**kwargs)
        return self.yambo_flow

class BandsConvergenceFlow():
    """
    This class takes as an input one structure.
    It produces a flow with the QE input files for nscf calculation with the max number of bands
    Then adds one yambo task with different number of bands

    Author: Henrique Miranda
    """
    def __init__(self,structure):
        self.structure = structure
    
    def get_tasks(self,scf_kpoints,ecut,nscf_kpoints,bands_list,**kwargs):
        """
        Create a flow with all the tasks to perform the calculation
        
        Arguments:
            generator: a builder function that takes structure, kpoints, ecut, 
                       nscf_bands, nscf_kpoints and kwargs as argument and returns the
                       tasks to be performed at each displacement
        """
        generator = kwargs.pop("generator",YamboIPChiTask)
        conv_thr = kwargs.pop("conv_thr",qepyenv.CONV_THR)

        #create a QE scf task and run
        nscf_bands = max(bands_list)
        tasks = PwNscfTasks(self.structure,kpoints=scf_kpoints,ecut=ecut,
                            nscf_bands=nscf_bands,nscf_kpoints=nscf_kpoints,conv_thr=conv_thr)
        qe_scf_task, qe_nscf_task, p2y_task = tasks
        tasks = list(tasks)

        link_task = p2y_task
        new_task  = p2y_task
        #we invert the list and compute from higher bands to lower
        #this avoids recomputing the dipoles
        for bands in sorted(bands_list,reverse=True):
            #generate tasks
            new_task = generator(link_task,bands=bands,dependencies=new_task,**kwargs)
            tasks.append(new_task)
            link_task = [p2y_task,new_task]

        return tasks

    def get_flow(self,path,scf_kpoints,ecut,nscf_kpoints,bands_list,**kwargs):
        tasks = self.get_tasks(scf_kpoints=scf_kpoints,ecut=ecut,nscf_kpoints=nscf_kpoints,
                               bands_list=bands_list,**kwargs)
       
        #put all the tasks in a flow
        self.yambo_flow = YambopyFlow.from_tasks(path,tasks,**kwargs)
        return self.yambo_flow


#inject some code in the initialization
def copy_pp(self):
    """
    This code is to be executed during the initialization of the SOC calculation
    It uses the dielectric function calculated without spin orbit coupling to
    calculate the QP corrections with spin-orbit coupling
    """
    import os
    import glob
    import shutil
    from netCDF4 import Dataset
    nospin_task = self.get_vars('nospin_task')
    spin_bands = self.get_vars('spin_bands')
    src = os.path.join(nospin_task.path,'run')
    dst = os.path.join(self.path,'run')
    if not os.path.isdir(dst): os.mkdir(dst)

    #link all the pp files
    ppfiles = glob.glob(os.path.join(src,'ndb.pp*'))
    for ppfile_path in ppfiles:
        ppfile = os.path.basename(ppfile_path)
        shutil.copy(ppfile_path,os.path.join(dst,ppfile))

    #change the spin value in the main file
    #also need to change the number of bands so yambo will read it
    main_pp = os.path.join(dst,'ndb.pp')
    with Dataset(main_pp,'r+') as dst_db:
        dst_db.variables['SPIN_VARS'][1]=2
        dst_db.variables['X_PARS_1'][2]=spin_bands


def get_scissor(self):
    """
    This code is to be executed during the intialization of a BSE calculation
    with scissor operator. It computes the scissor shift from a QP database
    and sets the BSE input accordingly
    """
    import os
    from yambopy.dbs.qpdb import YamboQPDB
    qp_task = self.get_vars('qp_task')
    valence = self.get_vars('valence')
    qp_path = os.path.join(qp_task.path,'run')
    #read qp database
    qpdb = YamboQPDB.from_db(folder=qp_path)
    scissor = qpdb.get_scissor(valence,verbose=0)[:3]
    #set scissor in the current input file
    self.yamboin_dict['KfnQP_E'] = list(scissor)

class SpinOrbitFlow():
    """
    Perform a yambo calculation including spin-orbit coupling.
    The dielectric function is calculated from the system without spin-orbit
    This has been tested for GW and BSE

    Author: Henrique Miranda
    """
    def __init__(self,structure,structure_spin=None):
        self.structure_nospin = structure
        if structure_spin is None: self.structure_spin = structure
  
    def get_tasks(self,scf_kpoints,ecut,nscf_kpoints,chi_bands,spin_bands,**kwargs):
        """
        Get a list of tasks executing this flow
        """
        tasks=[]
        conv_thr = kwargs.pop("conv_thr",qepyenv.CONV_THR)

        #create a yambo qp run
        yamboin_default_dict = dict(BndsRnXp=[1,chi_bands],
                                    NGsBlkXp=[1,'Ry'],
                                    EXXRLvcs=[10,'Ry'],
                                    QPkrange=[1,1,1,spin_bands],
                                    GbndRnge=[1,spin_bands])

        yamboin_dict = kwargs.pop("yamboin_dict",yamboin_default_dict)
        generator = kwargs.pop("generator",YamboQPTask)
        pp_runlevel = kwargs.pop("pp_runlevel",'-p p -V all')
        qp_runlevel = kwargs.pop("qp_runlevel",'-p p -g n -V all')
        spin_runlevel = kwargs.pop("spin_runlevel",qp_runlevel) 

        #without spin
        new_tasks = PwNscfTasks(self.structure_nospin,kpoints=scf_kpoints,ecut=ecut,
                            nscf_bands=chi_bands,nscf_kpoints=nscf_kpoints,conv_thr=conv_thr)
        qe_scf_task, qe_nscf_task, p2y_task = new_tasks

        nospin_task = generator(p2y_task,runlevel=pp_runlevel,yamboin_dict=yamboin_dict,dependencies=p2y_task,**kwargs)

        tasks.extend( [qe_scf_task, qe_nscf_task, p2y_task, nospin_task] )

        #with spin
        new_tasks = PwNscfTasks(self.structure_spin,kpoints=scf_kpoints,ecut=ecut,
                            nscf_bands=spin_bands,nscf_kpoints=nscf_kpoints,conv_thr=conv_thr)
        qe_scf_task, qe_nscf_task, p2y_task = new_tasks

        #patch the input files to include spin-orbit coupling
        qe_scf_task.get_instances_from_inputs(PwIn)[0].set_spinorbit()
        qe_nscf_task.get_instances_from_inputs(PwIn)[0].set_spinorbit()

        spin_task = generator(p2y_task,runlevel=spin_runlevel,yamboin_dict=yamboin_dict,dependencies=p2y_task,**kwargs)

        tasks.extend( [qe_scf_task, qe_nscf_task, p2y_task, spin_task] )

        spin_task.set_vars('spin_bands',spin_bands)
        spin_task.set_vars('nospin_task',nospin_task)
        spin_task.set_code('initialize',copy_pp)

        return tasks

    def get_flow(self,path,scf_kpoints,ecut,nscf_kpoints,chi_bands,spin_bands,**kwargs):
        tasks = self.get_tasks(scf_kpoints=scf_kpoints,ecut=ecut,nscf_kpoints=nscf_kpoints,
                               chi_bands=chi_bands,spin_bands=spin_bands,**kwargs)
       
        #put all the tasks in a flow
        self.yambo_flow = YambopyFlow.from_tasks(path,tasks,**kwargs)
        return self.yambo_flow

def PwNscfYamboIPChiTasks(structure,kpoints,ecut,nscf_bands,
                          yambo_runlevel='-o c -V all',nscf_kpoints=None,**kwargs):
    """
    Return the PwNscfTasks and a YamboTask for and IP calculation
    """ 
    conv_thr = kwargs.pop("conv_thr",qepyenv.CONV_THR)

    #create scf, nscf and p2y task
    tmp_tasks = PwNscfTasks(structure,kpoints,ecut,nscf_bands,nscf_kpoints,conv_thr=conv_thr)
    qe_scf_task,qe_nscf_task,p2y_task = tmp_tasks

    yambo_task = YamboIPChiTask(p2y_task,yambo_runlevel=yambo_runlevel,**kwargs)

    return qe_scf_task,qe_nscf_task,p2y_task,yambo_task

def YamboIPChiTask(p2y_task,**kwargs):
    """
    Return a yambo IP task
    """
    yambo_ip_default_dict = dict(QpntsRXd=[1,1],
                                 ETStpsXd=1000)
    yambo_ip_dict = kwargs.pop('yambo_ip_dict',yambo_ip_default_dict)
    yambo_runlevel = kwargs.pop('yambo_runlevel','-o c -V all')
    bands = kwargs.pop('bands',None)
    dependencies = kwargs.pop('dependencies',p2y_task)
    if bands: yambo_ip_dict['BndsRnXd'] = [1,bands]
    
    #add yambo_task
    yambo_task = YamboTask.from_runlevel(p2y_task,yambo_runlevel,yambo_ip_dict,
                                         dependencies=dependencies)
    return yambo_task

def YamboQPBSETasks(p2y_task,qp_dict,bse_dict,qp_runlevel='-p p -g n -V all',
                    bse_runlevel='-p p -k sex -y d -V all',valence=False,**kwargs):
    """
    Return a QP and BSE calculation
    """
    import copy
    dependencies = kwargs.pop('dependencies',p2y_task)
    inputs = [p2y_task]
    additional_inputs = kwargs.pop('additional_inputs',None)
    if additional_inputs: 
        if not isiter(additional_inputs): additional_inputs = [additional_inputs]
        inputs.extend(additional_inputs)

    #create a yambo qp run
    qp_task = YamboTask.from_runlevel(inputs,qp_runlevel,qp_dict,dependencies=dependencies)

    #create a yambo qp+bse run
    bse_dict['KfnQPdb']="E < run/ndb.QP"
    qp_bse_task = YamboTask.from_runlevel([p2y_task,qp_task],bse_runlevel,bse_dict,
                                          dependencies=qp_task)
    
    #create a yambo scissor+bse run
    if valence:
        bse_dict = copy.deepcopy(bse_dict)
        bse_dict.pop('KfnQPdb')
        scissor_bse_task = YamboTask.from_runlevel([p2y_task,qp_task],bse_runlevel,bse_dict,
                                                   dependencies=qp_task)
        scissor_bse_task.set_vars('qp_task',qp_task)
        scissor_bse_task.set_vars('valence',valence)
        scissor_bse_task.set_code('initialize',get_scissor)
        return qp_task, qp_bse_task, scissor_bse_task
    
    return qp_task, qp_bse_task


def YamboQPTask(p2y_task,yamboin_dict,runlevel='-p p -g n -V all',**kwargs):
    """
    Return a QP calculation
    """
    #create a yambo qp run
    dependencies = kwargs.pop('dependencies',p2y_task)
    qp_task = YamboTask.from_runlevel(p2y_task,runlevel,yamboin_dict,dependencies=p2y_task)

    return qp_task

def PhPhononTasks(structure,kpoints,ecut,qpoints=None):
    """
    Return a ScfTask, a PhTask and Matdyn task
    """

    #create a QE scf task and run
    qe_input = PwIn.from_structure_dict(structure,kpoints=kpoints,ecut=ecut)
    qe_scf_task = PwTask.from_input(qe_input)

    #create phonon tasks
    if qpoints is None: qpoints = qe_input.kpoints
    ph_input = PhIn.from_qpoints(qpoints)
    ph_task = PhTask.from_scf_task([ph_input,qe_scf_task],dependencies=qe_scf_task)

    #create matdyn task
    matdyn_task = DynmatTask.from_phonon_task(ph_task,dependencies=ph_task)
 
    return qe_scf_task, ph_task, matdyn_task

def PwNscfTasks(structure,kpoints,ecut,nscf_bands,nscf_kpoints=None,**kwargs):
    """
    Return a ScfTask, NscfTask and P2yTask preparing for a Yambo calculation
    """
    scf_conv_thr = kwargs.pop("conv_thr",qepyenv.CONV_THR)
    scf_conv_thr = kwargs.pop("scf_conv_thr",scf_conv_thr)
    nscf_conv_thr = kwargs.pop("nscf_conv_thr",scf_conv_thr*10)

    #create a QE scf task and run
    qe_input = PwIn.from_structure_dict(structure,kpoints=kpoints,ecut=ecut,conv_thr=scf_conv_thr)
    qe_scf_task = PwTask.from_input(qe_input)

    #create a QE nscf task and run
    qe_input = qe_input.copy().set_nscf(nscf_bands,nscf_kpoints,conv_thr=nscf_conv_thr)
    qe_nscf_task = PwTask.from_input([qe_input,qe_scf_task],dependencies=qe_scf_task)

    #create a p2y nscf task and run
    p2y_task = P2yTask.from_nscf_task(qe_nscf_task)

    return qe_scf_task, qe_nscf_task, p2y_task


def AbinitNscfTasks(inp,kpoints,ecut,nscf_bands,nscf_kpoints=None,**kwargs):
    tolwfr = kwargs.pop("tolwfr",1e-22)
    nbdbuf = kwargs.pop("nbdbuf",int(nscf_bands*0.15))

    #scf
    inp = inp.deepcopy()
    inp.set_kmesh(ngkpt=kpoints,
                  kptopt=1,
                  shiftk=(0,0,0))
    abinit_scf = AbinitTask.from_input(inp)

    #nscf
    inp = inp.deepcopy()
    inp.set_kmesh(ngkpt=nscf_kpoints,
                  kptopt=1,
                  shiftk=(0,0,0))
    inp.set_vars(iscf=-2,
                 getden=1,
                 nband=nscf_bands+nbdbuf,
                 nbdbuf=nbdbuf)
    abinit_nscf = AbinitTask.from_input([inp,abinit_scf],dependencies=abinit_scf)

    #e2y task
    e2y_task = E2yTask.from_nscf_task(abinit_nscf)

    return abinit_scf_task, abinit_nscf_task, e2y_task
