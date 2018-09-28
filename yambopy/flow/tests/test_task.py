# Copyright (C) 2018 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
import unittest
import os
import shutil
from qepy.pw import PwIn
from yambopy.data.structures import BN, Si
from yambopy.io.factories import PhPhononTask, PwNscfTask
from yambopy.flow import YambopyFlow, PwTask, P2yTask, YamboTask 

test_path = os.path.join(os.path.dirname(__file__),'..','..','data','refs','bse')

class TestFlow(unittest.TestCase):

    def clean(self,flows):
        if not isinstance(flows,list): flows = [flows]
        for flow in flows:
            if os.path.isdir(flow): shutil.rmtree(flow) 

    def test_full_flow(self):
        self.clean('flow')

        qe_scf_task, qe_nscf_task, p2y_task = PwNscfTask(BN,kpoints=[3,3,1],ecut=20,nscf_bands=10)

        #create a yambo optics task and run
        yamboin_dict = dict(FFTGvecs=[20,'Ry'],
                            BndsRnXs=[[1,10],''],
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

    def test_flow_2_parts(self):
        """
        This is a two parts flow
        In the first flow we prepare the SAVE folder
        In the second flow we run a BSE convergence test with number of bands
        """
        #
        # SAVE Flow
        #
        self.clean('save_flow')

        qe_scf_task, qe_nscf_task, p2y_task = PwNscfTask(BN,kpoints=[3,3,1],ecut=20,nscf_bands=10)

        #create yamboflow
        yambo_flow = YambopyFlow.from_tasks('save_flow',[qe_scf_task,qe_nscf_task,p2y_task])
        yambo_flow.create()
        yambo_flow.dump_run()
        print(yambo_flow)
        yambo_flow.run()

        #
        # BSE convergence Flow
        #
        self.clean('bse_flow')

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

        yambo_flow = YambopyFlow.from_tasks('bse_flow',tasks)
        print(yambo_flow)
        yambo_flow.run()

    def test_phonon_flow(self):
        self.clean('phonon_flow')        

        tasks = PhPhononTask(Si,kpoints=[3,3,3],qpoints=[1,1,1],ecut=30)
        yambo_flow = YambopyFlow.from_tasks('phonon_flow',tasks)
        print(yambo_flow)
        yambo_flow.create()
        yambo_flow.run()
        
    def tearDown(self):
        self.clean(['flow','bse_flow','save_flow'])

if __name__ == '__main__':
    unittest.main()
