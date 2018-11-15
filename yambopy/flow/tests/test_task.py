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
from yambopy.io.factories import PhPhononTasks, PwNscfTasks, YamboQPBSETasks, KpointsConvergenceFlow
from yambopy.flow import YambopyFlow, PwTask, P2yTask, YamboTask 

test_path = os.path.join(os.path.dirname(__file__),'..','..','data','refs','bse')

class TestFlow(unittest.TestCase):

    def clean(self,flows):
        if not isinstance(flows,list): flows = [flows]
        for flow in flows:
            if os.path.isdir(flow): shutil.rmtree(flow) 

    def test_full_flow(self):
        self.clean('flow')

        qe_scf_task, qe_nscf_task, p2y_task = PwNscfTasks(BN,kpoints=[3,3,1],ecut=20,nscf_bands=10)

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

        qe_scf_task, qe_nscf_task, p2y_task = PwNscfTasks(BN,kpoints=[3,3,1],ecut=20,nscf_bands=10)

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

        tasks = PhPhononTasks(Si,kpoints=[3,3,3],qpoints=[1,1,1],ecut=30)
        yambo_flow = YambopyFlow.from_tasks('phonon_flow',tasks)
        print(yambo_flow)
        yambo_flow.create()
        yambo_flow.run()

    def fd_flow(self):    
        self.clean('fd_flow')

        phonon_modes = Matdyn.from_modes_file(folder='phonon_flow/t2',filename='pw.modes')
        print(phonon_modes)
        fd = FiniteDifferencesPhononFlow(Si,phonon_modes)
        yambo_input=dict(ETStpsXd=1000,
                         EnRngeXd=[[0,5],'eV'])
        fd_flow = fd.get_flow('fd_flow',imodes_list=[3,4,5],kpoints=[1,1,1],ecut=20,nscf_bands=10,
                               yambo_input=yambo_input,yambo_runlevel='-o c')
        fd_flow.create()
        fd_flow.run()
        print(fd_flow)
   
    def test_kpoint_flow(self):
        """ Run kpoint convergence flow for the number of kpoints """
        self.clean('kpoint_flow')
        kp = KpointsConvergenceFlow(Si)
        kp_flow = kp.get_flow('kpoint_flow',scf_kpoints=[4,4,4],
                              nscf_kpoints_list=[[i,i,i] for i in [2,4,6]],
                              ecut=20, nscf_bands=10)
        kp_flow.create()
        kp_flow.run()
        print(kp_flow) 

    def test_qpbse_flow(self):
        self.clean('qpbse_flow')

        tasks = []
        tmp_tasks = PwNscfTasks(Si,kpoints=[4,4,4],ecut=20,nscf_bands=10)
        qe_scf_task, qe_nscf_task, p2y_task = tmp_tasks
        tasks.extend(tmp_tasks)

        #create a yambo qp run
        qp_dict = dict(FFTGvecs=[20,'Ry'],
                       BndsRnXp=[1,20],
                       NGsBlkXp=[1,'Ry'],
                       EXXRLvcs=[10,'Ry'],
                       QPkrange=[1,13,3,6],
                       GbndRnge=[1,10])
       
        #create a yambo bse run
        bse_dict = dict(BEnSteps=1000,
                        BEnRange=[[0,5],'eV'],
                        BndsRnXp=[1,20],
                        NGsBlkXp=[1,'Ry'],
                        FFTGvecs=[20,'Ry'],
                        BSENGexx=[10,'Ry'],
                        BSENGBlk=[1,'Ry'],
                        KfnQPdb="E < run/ndb.QP",
                        BSEBands=[3,6])

        qp_task, bse_task = YamboQPBSETasks(p2y_task,qp_dict,bse_dict) 
        tasks.extend([qp_task, bse_task])

        yambo_flow = YambopyFlow.from_tasks('qpbse_flow',tasks)
        yambo_flow.create()
        print(yambo_flow)
        yambo_flow.run()


    def tearDown(self):
        #self.clean(['flow','bse_flow','save_flow','phonon_flow','fd_flow','qpbse_flow'])
        pass

if __name__ == '__main__':
    unittest.main()
