# Copyright (C) 2018 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
import unittest
import os
import shutil
from qepy.pw import PwIn
from yambopy.data.structures import BN
from yambopy.flow import YambopyFlow, PwTask, P2yTask, YamboTask 

test_path = os.path.join(os.path.dirname(__file__),'..','..','data','refs','bse')

class TestTask(unittest.TestCase):

    def setUp(self):
        if os.path.isdir('flow'): shutil.rmtree('flow') 

    def test_flow(self):
        #create a QE scf task and run
        qe_input = PwIn.from_structure_dict(BN,kpoints=[9,9,1],ecut=20)
        qe_scf_task = PwTask.from_input(qe_input)
        
        #create a QE nscf task and run
        qe_input = qe_input.copy().set_nscf(10)
        qe_nscf_task = PwTask.from_input([qe_input,qe_scf_task],dependencies=qe_scf_task)

        #create a p2y nscf task and run
        p2y_task = P2yTask.from_nscf_task(qe_nscf_task)

        #create a yambo optics task and run
        yamboin_dict = dict(FFTGvecs=[30,'Ry'],
                            BndsRnXs=[1,30],
                            QpntsRXd=[1,1],
                            ETStpsXd=500)

        yambo_task = YamboTask.from_runlevel(p2y_task,'-o c',yamboin_dict,dependencies=p2y_task)
        print(yambo_task)
        yambo_task.run()

        #create yamboflow
        yambo_flow = YambopyFlow.from_tasks('flow',[qe_scf_task,qe_nscf_task,p2y_task,yambo_task])
        print(yambo_flow)
        yambo_flow.create()
        #yambo_flow.run()
    
        #store for cleanup
        self.yambo_flow = yambo_flow

    def tearDown(self):
        #self.yambo_flow.clean()
        pass

if __name__ == '__main__':
    unittest.main()
