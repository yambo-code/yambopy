yambopy
=======

Create automatic workflows for yambo using python.

Yambo official website: http://www.yambo-code.org/ 

Yambo Github page: https://github.com/yambo-code/yambo

Documentation
-------------
You can read the documentation in:
http://yambopy.readthedocs.org/en/latest/

Features
--------
- Create Yambo and Quantum Espresso input files from python
- Collect the output data in .json files for analysis
- Plot the results using matplotlib
- Test suite
- Tutorial

Tests
------
- Generate input files and compare to reference
- Generate input files, run the calculations using Quantum Espresso and Yambo
- Plot the results using matplotlib

Requirements
------------
- yambo (>4.0.0): http://www.yambo-code.org/
- numpy: http://www.numpy.org/
- matplotlib: http://matplotlib.org/
- netCDF4: http://unidata.github.io/netcdf4-python/
- Quantum Espresso (optional): http://www.quantum-espresso.org/

TODO
----
- Automatic convergence tests (increase a certain variable until the final result changes less than a certain threshold)

Authors
------
- [Henrique Pereira Coutada Miranda](http://henriquemiranda.github.io/)
- [Alejandro Molina Sanchez](http://alexmoratalla.github.io/)
- Fulvio Paleari
- Alexandre Morlet

The code is at an early stage of development, help us by sending bug reports, patches and suggestions!
