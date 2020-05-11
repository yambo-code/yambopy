Tutorial
========

UNDER DEVELOPMENT

Boron Nitride REAL TIME
==============

0. Calculate the Ground state properties of boron nitride using Quantum espresso (gs_bn.py)
    - TO RUN SCF AND NSCF: python gs_bn.py -sn


1. Setup the RT yambo database automatically (prepare_rt.py)
    - TO RUN RT setup: python prepare_rt.py -f E_x E_y E_z -p qe_prefix

2. Run convergence tests for time steps (optimize_time_step.py)
    - TO RUN RT time step optimization: python optimize_time_step.py -F input_file_path



