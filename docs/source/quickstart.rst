Quickstart
==========

Were we will show some small snippets of code as examples of how to tuse the code.
The general idea is that yambopy provides a class for each step of the calculation.

qepy
------------------
``qepy`` is the module to handle the `Quantum espresso <http://www.quantum-espresso.org/>` part of the calculation
Here we will start with qepy to generate the input files for the scf and nscf calculations.

PwIn
~~~~~~~~~~~~~~~~~~

Here we show how to create the input file for a scf calculation

.. code-block:: python
    
    from qepy import *
    
    #create input file from scratch
    qe = PwIn()

    #input structure
    qe.atoms = [['Si',[0.125,0.125,0.125]],
                ['Si',[-.125,-.125,-.125]]]
    qe.atypes = {'Si': [28.086,"Si.pbe-mt_fhi.UPF"]}

    #control variables
    qe.control['prefix'] = "'si'" #strings need double "''"
    qe.control['wf_collect'] = '.true.' #logicals

    #system
    qe.system['celldm(1)'] = 10.3
    qe.system['ecutwfc'] = 30
    qe.system['occupations'] = "'fixed'"
    qe.system['nat'] = 2
    qe.system['ntyp'] = 1
    qe.system['ibrav'] = 2

    #electrons
    qe.electrons['conv_thr'] = 1e-8

    #write file
    qe.write('si.scf')


schedulerpy
------------------
``schedulerpy`` is a module to abstract the job execution of the environment.
Here is the example to run the scf and nscf calculation using the input files we just created.

.. code-block:: python
    
    from schedulerpy import *

    # scheduler 1 node and 4 cores
    shell = Scheduler.factory(nodes=1,cores=4)

    # scheduler of pbs type
    shell = Scheduler.factory(scheduler='pbs')

    #add commands to the shell
    shell.add_command("echo 'hello world'")

    #view commands on the screen
    print( shell ) 

    #write to a file
    shell.write("commands.sh") 

    #submit or run the job
    shell.run() 


yambopy
-----------
``yambopy`` is the module to read/write input files, and read output files from yambo.
More recently we started to read some netcdf databases created by yambo.
Here we will show how to create input files for a GW and BSE calculations.

YamboIn
~~~~~~~~~~

.. code-block:: python

    from yambopy import *

    #creathe input file in 'bse' folder with SAVE
    y = YamboIn('yambo -b -o b -k sex -y d -V all',folder='bse')

    # define variables
    y['FFTGvecs'] = [30,'Ry'] # scalar + units
    y['BndsRnXs'] = [1,30] # array with integers
    y['BSEBands'] = [3,6] # array with integers
    y['BEnRange'] = [[0,8],'eV'] # array + units
    y['BEnSteps'] = 500 # numbers
    y['KfnQPdb'] = 'E < yambo/ndb.QP' #strings

    #write the file 
     y.write('bse/yambo_run.in')

    #create ypp input file
    y = YamboIn('ypp -e -a -V all',filename='ypp.in')

    #read local file
    y = YamboIn(filename='bse/yambo_run.in')

    #analyse the data
    ya = YamboAnalyser(folder)
    print(ya)

    # plot eel and eps from BSE
    ya.plot_bse('eel')
    ya.plot_bse('eps')

YamboAnalyser
~~~~~~~~~~~~~~~~~~
Here we read the GW band-structure that we just calculated.
And the BSE absorption spectra.


yambopy (bash)
--------------------------
We also include a python script that can be added to the PATH and executed from
the unix command line and provides many features of yambopy in a direct way.

.. code-block:: bash

    $ yambopy #lists all possible commands
    $ yambopy plotem1s #help about this command

