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

class TestFlow(unittest.TestCase):

    def setUp(self):
        flows_to_clean = ['flow','save_flow','bse_flow']
        for flow in flows_to_clean:
            if os.path.isdir(flow): shutil.rmtree(flow) 

    def test_full_flow(self):
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
                            BndsRnXs=[[1,8],''],
                            QpntsRXd=[[1,1],''],
                            ETStpsXd=[500,''])

        yambo_task = YamboTask.from_runlevel(p2y_task,'-o c',yamboin_dict,dependencies=p2y_task)
        print(yambo_task)

        #create yamboflow
        yambo_flow = YambopyFlow.from_tasks('flow',[qe_scf_task,qe_nscf_task,p2y_task,yambo_task])
        print(yambo_flow)
        yambo_flow.create()
        yambo_flow.dump_run()
        print(yambo_flow)
        yambo_flow.run()
    
        #store for cleanup
        self.yambo_flow = yambo_flow

    def test_flow_2_parts(self):
        """
        This is a two parts flow
        In the first flow we prepare the SAVE folder
        In the second flow we run a BSE convergence test with number of bands
        """
        #
        # SAVE Flow
        #

        #create a QE scf task and run
        qe_input = PwIn.from_structure_dict(BN,kpoints=[9,9,1],ecut=20)
        qe_scf_task = PwTask.from_input(qe_input)
        
        #create a QE nscf task and run
        qe_input = qe_input.copy().set_nscf(10)
        qe_nscf_task = PwTask.from_input([qe_input,qe_scf_task],dependencies=qe_scf_task)

        #create a p2y nscf task and run
        p2y_task = P2yTask.from_nscf_task(qe_nscf_task)

        #create yamboflow
        yambo_flow = YambopyFlow.from_tasks('save_flow',[qe_scf_task,qe_nscf_task,p2y_task])
        yambo_flow.create()
        yambo_flow.dump_run()
        print(yambo_flow)
        yambo_flow.run()

        #
        # BSE convergence Flow
        #
        #create a yambo optics task and run
        tasks = [p2y_task]
        for BSEEhEny in [5,6,7,8]:
            #create a list of yambo optics task and run
            yamboin_dict = dict(NGsBlkXs=[1,'Ry'],
                                BndsRnXs=[1,10],
                                BSEBands=[1,10],
                                BEnRange=[[0.0,6.0],'eV'],
                                BSEEhEny=[[0,BSEEhEny],'eV'],
                                BEnSteps=1000)

            yambo_task = YamboTask.from_runlevel(p2y_task,'-b -o b -k sex -y d -V all',
                                             yamboin_dict,dependencies=p2y_task)

            tasks.append(yambo_task)

        yambo_flow = YambopyFlow.from_tasks('bse_flow',[p2y_task,yambo_task])
        print(yambo_flow)
        yambo_flow.run()

    def tearDown(self):
        #self.yambo_flow.clean()
        pass

if __name__ == '__main__':
    unittest.main()
