yambopy
=======

Create automatic workflows for yambo and quantum espresso using python.
Do pre/post-processing, data analysis and plotting for yambo and quantum espresso.

Yambo official website: http://www.yambo-code.eu/ 

Yambo Github page [download yambo & yambopy]: https://github.com/yambo-code/yambo


Documentation
-------------
You can find explained tutorials and a partial 
documentation on the Yambo wiki page: http://www.yambo-code.org/wiki

Features
--------
- Create Yambo and Quantum Espresso input files from python
- Collect and manipulate the human-readable output data for analysis
- Automatic submissions of calculations (e.g., custom workflows for convergence or multi-executable runs)
- Access Yambo netCDF databases and Quantum Espresso xml files
- Analyse, interpolate and plot the results in various ways using matplotlib
- Visualize advanced quantities such as:
  -  dielectric function, exciton weights in k and q-space, electron-phonon matrix elements...
- Tutorials

Requirements
------------
- yambo (>5.0.0): http://www.yambo-code.org/
- numpy: http://www.numpy.org/
- scipy: https://www.scipy.org/
- matplotlib: http://matplotlib.org/
- netCDF4: http://unidata.github.io/netcdf4-python/
- pyyaml: https://pyyaml.org/
- Quantum Espresso (optional): http://www.quantum-espresso.org/
- Abipy (optional): https://abinit.github.io/abipy/

TODO
----
- Review and update of all features
- Full documentation & tutorials
- Test suite
- Make it easier to add new features

Authors
------
Yambopy was started and initially developed by [Henrique Pereira Coutada Miranda](http://henriquemiranda.github.io/).

Current developers and maintainers:
- [Alejandro Molina Sanchez](http://alexmoratalla.github.io/)
- [Fulvio Paleari](http://palful.github.io)

Collaborators include(d)
- Matteo Zanfrognini
- Jorge Cervantes
- Riccardo Reho
- Alexandre Morlet
- Davide Romanin
- Michele Re Fiorentin
- You if you want to share your scripts!

The code is at an ongoing stage of development, help us by sending bug reports, patches and suggestions!
