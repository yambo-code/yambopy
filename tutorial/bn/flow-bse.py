# Copyright (C) 2019 Alejandro Molina Sanchez - Henrique PC Miranda 
# All rights reserved.
#
# This file is part of yambopy
#
# Tutorial File of Yambopy Tasks. BSE flow
# 

#import argparse
#import os
#import shutil
from yambopy.flow import YambopyFlow, P2yTask, YamboTask
#from schedulerpy import Scheduler
from yambopy import yambopyenv

# Set list of task and dictionary of yambo variables
tasks = []
yamboin_dict = dict()

# Set origin of SAVE folder
p2y_task = P2yTask.from_folder('nscf_flow/t2')

print(p2y_task)

# Coulomb-cutoff and RIM dictionary
cutoffdict = dict(CUTBox = [0,0,10],CUTGeo='box z',RandQpts=1000000,RandGvec=[1,'RL'])

# Parallel Environment dictionary
paradict = dict(X_all_q_ROLEs="q",X_all_q_CPU="2")

# BSE variables dictionary
bsedict = dict(BEnSteps=1000,
                FFTGvecs=[10,'Ry'],
                BEnRange=[[0,5],'eV'],
                BndsRnXs=[1,60],
                NGsBlkXs=[1,'Ry'],
                BSENGexx=[10,'Ry'],
                BSENGBlk=[1,'Ry'],
                BSEBands=[7,10])

# Merge all dict variables
yamboin_dict = {**yamboin_dict,**cutoffdict,**paradict,**bsedict}

# Set Yambo task (BSE in this case)
# yamboin_args >> Add arguments (ExtendOut, WRbsWF, EvalKerr, etc.)

bse_task = YamboTask.from_runlevel([p2y_task],'-r -o b -b -k sex -y d -V all',yamboin_dict,yamboin_args=['WRbsWF'])

# Introduce each task in the list of task
tasks.append(bse_task)

# Set the Yambo flow
yambo_flow = YambopyFlow.from_tasks('bse_flow',tasks)
print(yambo_flow)

# Create the Yambo flow
yambo_flow.create(agressive=True)
# Run the Yambo flow
yambo_flow.run()
print(yambo_flow)
