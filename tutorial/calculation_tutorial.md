Tutorial on managing calculations
========

Here you find a basic tutorial to get you started on using yambopy.
Run the files, see what happens.
The idea of yambopy is to call qe and yambo to generate the base input files, read them into python classes allowing the user to manipulate them and then run series of calculations (i.e., convergence tests). Finally, yambopy also manages data analysis and plotting.

Boron Nitride GW+BSE
==============

0. Calculate the Ground state properties of boron nitride using Quantum espresso (gs\_bn.py)
    - Relax unit cell
    - Self-consistent cycle
    - Non self-consistent cycle
    - Band structure
    - Phonon dispersion (DFPT)

1. Generate the yambo databases automatically (included in the scripts)
    - Run p2y
    - Run Yambo

2. GW calculation for boron nitride (gw\_bn.py, plot-qp.py)
    - Set the variables for a yambo input file using python
    - Run the calculation

3. Convergence of GW calculation for boron nitride (gw\_conv\_bn.py, plot-gw-conv.py)
    - Set a python dictionary with different values for the variables to converge
    - Run multiple calculations
    - Plot the results

4. IP calculation for boron nitride (ip\_bn.py)
    - Set the variables for a yambo input file using python
    - Run the calculation

5. BSE calculation for boron nitride (bse\_bn.py, plot-bse.py, plot-excitondb.py)
    - Set the variables for a yambo input file using python
    - Run the calculation

6. Convergence of BSE calculation for boron nitride (bse\_conv\_bn.py, bse\_cutoff.py, plot-bse-conv.py)
    - Set a python dictionary with different values for the variables to converge
    - Run multiple calculations
    - Plot the results

7. GW+BSE calculation for boron nitride (gw\_bse\_bn.py)
    - Run full GW+BSE calculation using yambo

8. TO BE DONE: BSE parallelisation and job submission for boron nitride (gw\_par\_bn.py, plot-par.py)
    - Submit job to HPC cluster
    - Manage parallelisation
    - Plot CPU and memory usage
