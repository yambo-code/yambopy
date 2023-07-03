yambopy
=======

Create automatic workflows for yambo and quantum espresso using python. Work directly with netCDF databases.
Do pre/post-processing, data analysis and plotting for yambo and quantum espresso.

Yambo official website: http://www.yambo-code.eu/ 

Yambo Github page [download yambo & yambopy]: https://github.com/yambo-code/yambo


Documentation
-------------
You can find explained tutorials and a partial 
documentation on the Yambo wiki page: http://www.yambo-code.org/wiki
Additional information about capabilities and usage are available inside the tutorial folder and by reading the docstrings of the various classes.

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

Installation
------------
Follow the steps on the Yambo wiki: https://www.yambo-code.eu/wiki/index.php/Tutorials#Setting_up_Yambopy

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

Troubleshooting, bugs and questions
-----------------------------------
Please write a post in the yambopy subsection of the yambo forum: https://www.yambo-code.eu/forum/viewforum.php?f=35&sid=77b7f6076dea7cdf40432efbc035feb6

Current development goals
-------------------------
- General review and update of all features
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

How to contribute
-----------------
If you want to contribute, we suggest the following steps:
1. Fork this repository
2. Implement and test your new feature(s) in the forked repo
3. Create a pull request in order to include your development in the official code
