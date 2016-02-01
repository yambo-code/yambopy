Introduction
=============

A typical yambo calculation proceeds as follows:
    - Obtain the ground state proprieties from a DFT code (`pw.x` or `abinit`)
    - Create the yambo netCDF databases using the correspooonding interface: (`p2y` for `pw.x` or `a2y` for `abinit`)
    - Run `yambo` once to complete the database
    - Run `yambo` specifying the `runlevels <http://www.yambo-code.org/input_file/yambo_3.4.0/index.php>`_

    - Edit the `yambo` input file
    - Run `yambo`
    - Plot the data results

Since many of the parameters of the calculation have to be converged the user might end up running the last three steps many times.
This is rather time consuming without an automatization script.

The `yambo-py` project aims to provide a simple set of python scripts to read and edit yambo input files. The primary objective is to make the convergence tests easier.

In the future `yambopy` might be used to run `yambo` automatically on large sets of different materials.
The facilities to read and stor output files can also be used for consistency checks between different versions of the `yambo` code.
