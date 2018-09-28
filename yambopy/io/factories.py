# Copyright (C) 2018 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
"""
This file contains classes and methods to generate input files
"""

from qepy.pw import PwIn
from qepy.ph import PhIn
from yambopy.flow import PwTask, PhTask, P2yTask

__all__ = [
"FiniteDifferencesPhonon",
"PwNscfTask",
"PhPhononTask"
]

class FiniteDifferencesPhonon():
    """
    This class takes as an input one structure and a phonon calculation.
    It produces a flow with the QE input files displaced along the phonon modes
    """
    def __init__(self,structure,phonons):
        self.structure = structure
        self.phonons = phonons
        self.spectra = None

    def write_flow(self,path,modes_list=None):
        """
        Create a flow with all the tasks to perform the calculation
        """
        raise NotImplementedError("Not implemented yet")
        #apply the displacement in the structure

        #create scf, nscf and p2y task
        PwNscfTasks(displaced_structure,kpoints,ecut,nscf_bands)

        #put all the tasks in a flow
        yambo_flow = YambopyFlow.from_tasks(path,[qe_scf_task,qe_nscf_task,p2y_task,yambo_task])
        return yambo_flow


def PhPhononTask(structure,kpoints,ecut,qpoints=None):
    """
    Return a ScfTask and a series of phonon tasks
    """

    #create a QE scf task and run
    qe_input = PwIn.from_structure_dict(structure,kpoints=kpoints,ecut=ecut)
    qe_scf_task = PwTask.from_input(qe_input)

    #create phonon tasks
    if qpoints is None: qpoints = qe_input.kpoints
    ph_input = PhIn.from_qpoints(qpoints)
    ph_task = PhTask.from_scf_task([ph_input,qe_scf_task],dependencies=qe_scf_task)

    return qe_scf_task, ph_task 

def PwNscfTask(structure,kpoints,ecut,nscf_bands):
    """
    Return a ScfTask, NscfTask and P2yTask preparing for a Yambo calculation
    """

    #create a QE scf task and run
    qe_input = PwIn.from_structure_dict(structure,kpoints=kpoints,ecut=ecut)
    qe_scf_task = PwTask.from_input(qe_input)

    #create a QE nscf task and run
    qe_input = qe_input.copy().set_nscf(nscf_bands)
    qe_nscf_task = PwTask.from_input([qe_input,qe_scf_task],dependencies=qe_scf_task)

    #create a p2y nscf task and run
    p2y_task = P2yTask.from_nscf_task(qe_nscf_task)

    return qe_scf_task, qe_nscf_task, p2y_task
