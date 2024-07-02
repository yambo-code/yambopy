Introduction
=============

A typical `yambo` calculation proceeds as follows:

    - Obtain the ground state proprieties from a DFT code (`pw.x` or `abinit`)
    - Create the `yambo` netCDF databases using the corresponding interface: (`p2y` for `pw.x` or `a2y` for `abinit`)
    - Run `yambo` once to complete the database
    - Run `yambo` specifying the `runlevels <http://www.yambo-code.org/input_file/yambo_3.4.0/index.php>`_
    - Edit the `yambo` input file
    - Run `yambo`
    - Plot the data results

Since many of the parameters of the calculation have to be converged the user might end up running the last three steps many times.
This is rather time consuming without an automatization script.

The `yambopy` project aims to provide a simple set of python scripts to read and
edit `yambo` input files. The primary objective is to make the convergence tests easier.

In the future `yambopy` might be used to run `yambo` automatically on large sets
of different materials.
The facilities to read and store output files can also be used for consistency
checks of the code.

Keep in mind that this code is still in **beta** version.
Any bug reports and feedback are welcome!
You can report them at:
https://github.com/henriquemiranda/yambopy/issues

The code is hosted in Github:
https://github.com/henriquemiranda/yambopy
