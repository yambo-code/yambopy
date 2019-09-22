# Copyright (C) 2019 Alejandro Molina Sanchez - Henrique PC Miranda 
# All rights reserved.
#
# This file is part of yambopy
#
# Tutorial File of Yambopy Tasks. BSE flow
# 

import argparse
import os
import shutil
#from yambopy.data.structures import Si
#from qepy import PwXML
#from qepy.lattice import Path
#from qepy.matdyn import Matdyn
#from yambopy.io.factories import PwNscfTasks, PwBandsTasks, PwRelaxTasks
from yambopy.flow import YambopyFlow, P2yTask, YamboTask
from schedulerpy import Scheduler
from yambopy import yambopyenv

# Set list of task and dictionary of yambo variables
tasks = []
yamboin_dict = dict()


# Set origin of SAVE folder
p2y_task = P2yTask.from_folder('nscf_flow/t2')

print(p2y_task)

# Coulomb-cutoff and RIM dictionary
cutoffdict = dict(RandQpts=1000000,RandGvec=[1,'RL'])

# Parallel Environment dictionary
paradict = dict(X_all_q_ROLEs="q",X_all_q_CPU="2")

# BSE variables dictionary
bse_dict = dict(BEnSteps=1000,
                FFTGvecs=[10,'Ry'],
                BEnRange=[[0,5],'eV'],
                BndsRnXp=[1,20],
                NGsBlkXp=[1,'Ry'],
                BSENGexx=[10,'Ry'],
                BSENGBlk=[1,'Ry'],
                BSEBands=[2,7])

yamboin_dict = {**yamboin_dict,**cutoffdict,**paradict}

bse_task = YamboTask.from_runlevel([p2y_task],'-r -o b -b -k sex -y h -V all',yamboin_dict)

tasks.append(bse_task)

yambo_flow = YambopyFlow.from_tasks('bse_flow',tasks)
print(yambo_flow)

yambo_flow.create(agressive=True)
yambo_flow.run()
print(yambo_flow)
