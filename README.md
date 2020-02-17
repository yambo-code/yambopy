yambopy
=======

Create automatic workflows for yambo and quantum espresso using python.
Do pre/post-processing, data analysis and plotting for yambo and quantum espresso.

Yambo official website: http://www.yambo-code.org/ 

Yambo Github page: https://github.com/yambo-code/yambo

Documentation
-------------
You can read the documentation in:
http://yambopy.readthedocs.org/en/latest/

Features
--------
- Create Yambo and Quantum Espresso input files from python
- Collect and manipulate the output data for analysis
- Automatic submissions of calculations (e.g., convergence workflows)
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
- scipy: https://www.scipy.org/
- matplotlib: http://matplotlib.org/
- netCDF4: http://unidata.github.io/netcdf4-python/
- Quantum Espresso (optional): http://www.quantum-espresso.org/
- Abipy (optional): https://abinit.github.io/abipy/

TODO
----
- Enhance modularisation (task-oriented instead of goal-oriented) to support wider developments.

Authors
------
- [Henrique Pereira Coutada Miranda](http://henriquemiranda.github.io/)
- [Alejandro Molina Sanchez](http://alexmoratalla.github.io/)
- Fulvio Paleari
- Alexandre Morlet

The code is at an ongoing stage of development, help us by sending bug reports, patches and suggestions!
