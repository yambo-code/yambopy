
![yambopy_text](docs/logos/yambopy_text.png)

yambopy
=======

Create automatic workflows for yambo and quantum espresso using python. Work directly with netCDF databases.
Do pre/post-processing, data analysis and plotting for yambo and quantum espresso.

- Yambo official website: http://www.yambo-code.eu/ 
- Yambo wiki: http://www.yambo-code.org/wiki
- Yambo Github page [download yambo & yambopy]: https://github.com/yambo-code/yambo


Documentation
-------------

The main usage of yambopy is by importing its modules in the user's own scripts, such as:
```
from yambopy import *
from qepy import *
```
There is also a command line interface feature. Just type
```
yambopy 
```
on terminal to see the options. Typing `yambopy [option]` will show the related help message.

You can find explained tutorials and a partial 
documentation on the Yambo wiki page: https://www.yambo-code.eu/wiki/index.php/First_steps_in_Yambopy

The tutorials contain examples scripts illustrating how to use some of the features: they are intended to be copied, modified and adapted to other use cases and to your ideas and needs. 

Additional information about capabilities and usage are available inside the tutorial folder and by reading the docstrings of the various classes. Keep in mind that a basic knowledge of python (`numpy` and `matplotlib` packages) will greatly help while using yambopy. 

Features
--------
- Create Yambo and Quantum Espresso input files from python
- Collect and manipulate the human-readable output data for analysis
- Automatic submissions of calculations (e.g., custom workflows for convergence or multi-executable runs)
- Access Yambo netCDF databases and Quantum Espresso xml files
- Analyse, interpolate and plot the results in various ways using matplotlib
- Visualize advanced quantities such as:
  -  dielectric function, exciton weights in k and q-space, electron-phonon matrix elements...
- [Aiida](https://github.com/aiidateam) plugin for Yambo-Aiida workflows 
- Tutorials

Installation
------------

Make sure that you have a suitable python environment (created for example with [conda](https://docs.conda.io/projects/miniconda/en/latest/) or [venv](https://docs.python.org/3/library/venv.html)).

#### Regular installation of released version
Type `pip install yambopy`

#### Local installation from this repository (for latest patches)
Clone this repository in your local machine or cluster, enter the directory and type `pip install .`
 
#### More information
Follow the installation steps on the [Yambo wiki](https://www.yambo-code.eu/wiki/index.php/First_steps_in_Yambopy).

Requirements
------------
- numpy: http://www.numpy.org/
- scipy: https://www.scipy.org/
- matplotlib: http://matplotlib.org/
- netCDF4: http://unidata.github.io/netcdf4-python/
- lxml: https://lxml.de/
- pyyaml: https://pyyaml.org/
- monty: https://pypi.org/project/monty/

Yambopy works for the following DFT/MBPT codes:
- yambo (>5.0.0): https://www.yambo-code.eu/
- Quantum Espresso (optional): http://www.quantum-espresso.org/

Troubleshooting, bugs and questions
-----------------------------------
Please write a post in the yambopy subsection of the [yambo forum](https://www.yambo-code.eu/forum/viewforum.php?f=35&sid=77b7f6076dea7cdf40432efbc035feb6).

Current development goals
-------------------------
- General review and update of all features
- Brillouin zone paths patch
- Full support for finite-momentum BSE postprocessing
- Full documentation & tutorials
- Increase efficiency of I/O for large database sizes and numbers
- Test suite
- Make it easier to add new features

Authors
------
Original author:
- [Henrique Pereira Coutada Miranda](http://henriquemiranda.github.io/).

Current developers and maintainers:
- [Fulvio Paleari](http://palful.github.io) (CNR - Nanoscience institute, Modena)
- [Alejandro Molina Sanchez](http://alexmoratalla.github.io/) (University of Valencia)
- José Castelo (University of Valencia) 

Active contributors:
- Claudio Attaccalite
- Miki Bonacci
- Jorge Cervantes-Villanueva
- Riccardo Reho
- Michele Re Fiorentin
- You if you want to share your scripts!

Past contributors:
- Matteo Zanfrognini
- Alexandre Morlet
- Davide Romanin
- Daniel Murphy

The code is at an ongoing stage of development, help us by sending bug reports, patches and suggestions!

How to contribute
-----------------
If you want to contribute, we suggest the following steps:
1. Fork this repository
2. Implement and test your new feature(s) in the forked repo
3. Create a pull request in order to include your development in the official code

Acknowledgements
----------------
- The [Abipy](https://abinit.github.io/abipy/) library developed for the Abinit code was the original inspiration for Yambopy. In particular, abipy's `SkwInterpolator` module for band structure interpolations has been directly imported into yambopy. 
- Yambopy logos by Claudia Cardoso
- University of Luxembourg
- University of Valencia
- Nanoscience Institute of the Italian National Research Council
- [MaX](https://www.max-centre.eu/): Materials at the eXascale EU center of excellence

![yambopy_logo](docs/logos/yambopy_square.png)
