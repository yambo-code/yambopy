Tutorial
========

UNDER DEVELOPMENT

Boron Nitride REAL TIME
==============

0. Calculate the Ground state properties of boron nitride using Quantum espresso (gs_bn.py)
    - TO RUN SCF AND NSCF: python gs_bn.py -sn


1. Setup the RT yambo database automatically (prepare_rt.py)
    - TO RUN RT setup: call function YamboRTSetup(field_direction,prefix[,OPTIONAL VARIABLES])

    Info on YamboRTSetup (../yambopy/rt/rt_setup.py):

    Class to run the setup for RT calculations.
    Must be run outside the folder where the nscf calculation took place.
    Example of use:
    Generate a SAVE file with reduced symmetries:
        .. code-block:: python
            YamboRTSetup(FIELD_direction,QE_prefix,nscf=nscf_path,database=save_path)


    (file already prepared: prepare_rt.py -f E_x E_y E_z -p qe_prefix)

2. Run convergence tests for time steps (optimize_time_step.py)
    - TO RUN RT time step optimization: call function YamboRTStep_Optimize(input_path,SAVE_path[,OPTIONAL VARIABLES])

    Info on YamboRTStep_Optimize (../yambopy/rt/rt_timestep_optimize.py):

    Class to run convergence tests for the RT time step.
    Note: time steps must be given in as units.
    Example of use:
        .. code-block:: python
            YamboRTStep_Optimize(input_path,SAVE_path,RUN_path)

    (file already prepared: python optimize_time_step.py -F input_file_path)



