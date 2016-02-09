yambopy
=======

Create automatic workflows for yambo using python

Documentation
-------------
You can read the documentation in:
http://yambopy.readthedocs.org/en/latest/

Features
--------
- Create Yambo input files using a transparent python script
- Collect the output data in .json files for posterior analysis
- Plot the results using matplotlib
- Create Quantum Espresso input files using python
- Test suite
- Tutorial

Tests
------
- Generate input files and compare to reference
- Generate input files, run the calculations using Quantum Espresso and Yambo
- Plot the results using matplotlib

Requirements
------------
- numpy: http://www.numpy.org/
- matplotlib: http://matplotlib.org/
- netCDF4: http://unidata.github.io/netcdf4-python/

TODO
----
- Automatic convergence tests (increase a certain variable until the final result changes less than a certain threshold)

Authors
------
- Henrique Pereira Coutada Miranda
- Alejandro Molina Sanchez

The code is at an early stage of development, help us by sending bug reports, patches and suggestions!
