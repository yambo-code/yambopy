Tutorial
========

Here you find a basic tutorial to get you started on using yambopy.
Run the files, see what happens.
The idea of yambopy is to call qe and yambo to generate the base input files, read them into python classes allowing the user to manipulate them and then run series of calculations (i.e., convergence tests). Finally, yambopy also manages data analysis and plotting.

Boron Nitride GW+BSE
==============

0. Calculate the Ground state properties of silicon using Quantum espresso (gs_bn.py)
    - Relax unit cell
    - Self-consistent cycle
    - Non self-consistent cycle
    - Phonon dispersion (DFPT)

1. Generate the yambo databases automatically (included in the scripts)
    - Run p2y
    - Run Yambo

2. GW calculation for boron nitride (gw_bn.py, plot-qp.py)
    - Set the variables for a yambo input file using python
    - Run the calculation

3. Convergence of GW calculation for boron nitride (gw_conv_bn.py, plot-gw-conv.py)
    - Set a python dictionary with different values for the variables to converge
    - Run multiple calculations
    - Plot the results

4. IP calculation for boron nitride (ip_bn.py)
    - Set the variables for a yambo input file using python
    - Run the calculation

5. BSE calculation for boron nitride (bse_bn.py, plot-bse.py, plot-excitondb.py)
    - Set the variables for a yambo input file using python
    - Run the calculation

6. Convergence of BSE calculation for boron nitride (bse_conv_bn.py, bse_cutoff.py, plot-bse-conv.py)
    - Set a python dictionary with different values for the variables to converge
    - Run multiple calculations
    - Plot the results

7. GW+BSE calculation for boron nitride (gw_bse_bn.py)
    - Run GW calculation using yambo
    - Run BSE calculation using yambo using the dielectric function from the previous calculation

8. BSE parallelisation and job submission for boron nitride ()
    - TO BE DONE
