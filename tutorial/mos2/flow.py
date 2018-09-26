# Copyright (C) 2018 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
"""
 This is an example of a BSE calculation for MoS2 
 using the new Flow/Task methods of yambopy.
 The approach is the same as the one implemented in Abipy.
"""
from yambopy.data.structures import MoS2
from qepy.pw import PwIn
from yambopy.flow import YambopyFlow, PwTask, P2yTask, YamboTask

#create a QE scf task and run
qe_input = PwIn.from_structure_dict(MoS2,kpoints=[12,12,1],ecut=30)
qe_scf_task = PwTask.from_input(qe_input)

#create a QE nscf task and run
qe_input = qe_input.copy().set_nscf(20)
qe_nscf_task = PwTask.from_input([qe_input,qe_scf_task],dependencies=qe_scf_task)

#create a p2y nscf task and run
p2y_task = P2yTask.from_nscf_task(qe_nscf_task)

#create a yambo optics task and run
yamboin_dict = dict(NGsBlkXs=[1,'Ry'],
                    BndsRnXs=[[1,20],''],
                    BSEBands=[[8,11],''],
                    BEnRange=[[0.0,6.0],'eV'],
                    BEnSteps=[1000,''])

yambo_task = YamboTask.from_runlevel(p2y_task,'-b -o b -k sex -y d',
                                     yamboin_dict,dependencies=p2y_task)
print(yambo_task)

#create yamboflow
yambo_flow = YambopyFlow.from_tasks('flow',[qe_scf_task,qe_nscf_task,p2y_task,yambo_task])
print(yambo_flow)
yambo_flow.create()
yambo_flow.dump_run()
print(yambo_flow)
yambo_flow.run()
