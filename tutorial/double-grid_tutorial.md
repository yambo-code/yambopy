Tutorial
========

How to perform a double grid workflow and optionally converge it against optical absorption

- Examples are in: tutorial/double-grid
- Yambopy code is at: yambopy/double_grid, yambopy/io, ./materials

## Input data
- Path to executables pw, p2y, yambo, ypp [OPTIONAL]
- Material prefix for quantum espresso
- List of coarse grids CG_1, CG_2, ... CG_N
- List of fine grids FG_1i, FG_2i, .. FG_Mi for each CG_i
- Energy of laser impinging on the sample [OPTIONAL: if converging]
- Quantum espresso save folder of an scf calculation (previously computed) and path to it
- Path to pseudopotentials
- Path to work directory
- Scheduler objects for quantum espresso and yambo [OPTIONAL: if submitting on HPC facility]
- Base input file for nscf (qe) and independent-particle (yambo) calculations [Can be taken from ./materials]
- Prefix of the output file name for qe and yambo (the report file).

## Workflow

1. Use the function YamboDG_Optimize
    - Optional variables:
        - Run individual steps of the workflow (see below) 
        - Generate folder trees and inputs without running qe and yambo
    - NSCF calculations are in nscf_grids
    - IP calculations are in [yambo_calc_type]\_grids
    - Plots of the results - if converging - are in plots
    - Ready made example: dg_test.py script (using monolayer hBN)
    
2. Minimal python script to run the bn tutorial (assuming all executables are in the PATH):

 .. code-block:: python

    from yambopy import *

 .. code-block:: All the non-optional inputs

    # Standard workflow
    YamboDG_Optimize(cg_grids,fg_grids,prefix,qe_input,yambo_input,STEPS='all',\
                     scf_save_path=scf_save_path,pseudo_path=pseudo_path,RUN_path=work_dir,\
                     nscf_out=nscf_out,y_out_dir=y_out_dir,yambo_calc_type='ip')

 ..
    
    # Double grid convergence
    YamboDG_Optimize(cg_grids,fg_grids,prefix,qe_input,yambo_input,E_laser=E_laser,STEPS='all',converge_DG=True,\
                     scf_save_path=scf_save_path,pseudo_path=pseudo_path,RUN_path=work_dir,\
                     nscf_out=nscf_out,y_out_dir=y_out_dir)
 ..

3. Scheme of the workflow

    - The workflow is divided in FOUR STEPS that can be executed separately or together:
        1. nscf CG [STEPS='1']
        2. nscf FG and ip CG [STEPS='2']
        3. yambo FG [STEPS='3']
        4. plot results [STEPS='4' when converging]
        
    - Scheme of the workflow:
    
            NSCF                                       IP
            |                                          |
    step 1  CG_1              CG_2 ... CG_N            |
            |                 |                        |
    step 2  FG_11 ... FG_M1   FG_12 ... FG_M2 ...      CG_1              CG_2 ... CG_N
                                                       |                 |
    step 3                                             FG_11 ... FG_M1   FG_12 ... FG_M2 ...
                                                        \         |      |         /
                                                         \        \      |        /
    step 4                                                 _________ PLOTS ______
