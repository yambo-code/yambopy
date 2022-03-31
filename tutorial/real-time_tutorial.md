Tutorial
========

How to converge the time step of a Yambo RT simulation.

- Examples are in: tutorial/real-time
- Yambopy code is at: yambopy/rt, yambopy/io, yambopy/dbs

## Ground state calculation

0. Calculate the ground state properties of your system using Quantum espresso (scf and nscf runs).
    - Examples: python gs_bn.py -sn

## RT Convergence

1. Setup the RT yambo database automatically
    - call function YamboRTSetup(field_direction,prefix,[OPTIONAL VARIABLES])
    - Optional variables include setting nscf, SAVE, and yambo executable paths
    - Info on YamboRTSetup in yambopy/rt/rt_setup.py
    - Examples: python prepare_rt.py -f E_x E_y E_z -p qe_prefix

2. Run convergence tests for time steps (optimize_time_step.py)
    - call function YamboRTStep_Optimize(input_path,SAVE_path,TStep_MAX,TStep_increase,NSimulations,[OPTIONAL VARIABLES])
    - Optional variables include setting max time step, time step increase, max number of runs, run duration, tolerance for convergence tests
    - Info on YamboRTStep_Optimize in yambopy/rt/rt_timestep_optimize.py
    - Examples: python optimize_time_step.py -F input_file_path -D RUN_path

3. Minimal python script to run the bn tutorial:

 .. code-block:: python

    from yambopy import *

    YamboRTSetup([1,0,0],'bn') #Field direction and QE prefix

    YamboRTStep_Optimize('TD_inputs/td_ip.in','database/FixSymm/SAVE') #RT input and SAVE paths
 ..
