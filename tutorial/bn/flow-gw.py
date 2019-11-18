# Copyright (C) 2019 Alejandro Molina Sanchez - Henrique PC Miranda 
# All rights reserved.
#
# This file is part of yambopy
#
# Tutorial File of Yambopy Tasks. GW flow
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

# Parallel Environment dictionary (serial in this example)
paradict = dict(X_all_q_ROLEs="",X_all_q_CPU="",SE_CPU= "",SE_ROLEs= "")

# GW variables dictionary (standard variables, more advanced in Yambo Website)
gwdict = dict(FFTGvecs=[10,'Ry'],
              BndsRnXp=[1,60],
              NGsBlkXp=[1,'Ry'],
              GbndRnge=[1,60],
              EXXRLvcs=[10,'Ry'],
              VXCRLvcs=[10,'Ry'],
              QPkrange=[1,19,7,10])

# Merge all dict variables
yamboin_dict = {**yamboin_dict,**cutoffdict,**paradict,**gwdict}

# Set Yambo task (GW in this case)
# yamboin_args >> Add arguments (ExtendOut, WRbsWF, EvalKerr, etc.)

gw_task = YamboTask.from_runlevel([p2y_task],'-r -g n -p p -V all',yamboin_dict,yamboin_args=['ExtendOut'])

# Introduce each task in the list of task
tasks.append(gw_task)

# Set the Yambo flow
yambo_flow = YambopyFlow.from_tasks('gw_flow',tasks)
print(yambo_flow)

# Create the Yambo flow
yambo_flow.create(agressive=True)
# Run the Yambo flow
yambo_flow.run()
print(yambo_flow)
