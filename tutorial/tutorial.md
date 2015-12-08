Tutorial
========

Here you find a basic tutorial to get you started on using yambopy.
Run the files, see what happens.

0. Calculate the Ground state properties of silicon using Quantum espresso (gs_si.py)
    - Relax unit cell
    - Self-consistent cycle
    - Non self-consistent cycle

1. Generate the yambo databases
    - Run p2y
    - Run Yambo

2. GW calculation for silicon (gw_si.py)
    - Set the variables for a yambo input file using python
    - Run the calculation

3. Convergence of GW calculation for silicon (gw_conv_si.py)
    - Set a python dictionary with different values for the variables to converge
    - Run multiple calculations
    - Plot the results

4. BSE calculation for silicon (bse_si.py)
    - Set the variables for a yambo input file using python
    - Run the calculation

5. BSE calculation for silicon (bse_conv_si.py)
    - Set a python dictionary with different values for the variables to converge
    - Run multiple calculations
    - Plot the results

6. GW+BSE calculation for silicon (gw_bse_si.py)
    - Run GW calculation using yambo
    - Run BSE calculation using yambo using the dielectric function from the previous calculation
