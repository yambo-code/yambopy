Tutorial
==========

In the `tutorial` folder you can find some examples fo how to get started using `yambopy`.
The first step in any calculation with yambo is to calculate the ground state proprieties using either `abinit` or `pw.x`.
We don't have support to read and write `abinit` input files. To do that you should use the `abipy <https://github.com/gmatteo/abipy>`_ package:

Si
----

1. Ground State

Go to the `tutorial` folder and run the ground state calculation using the `gs_si.py` file:

.. code-block:: bash

    python gs_si.py

The script will run a relaxation of the structure, read the optimized cell parameter and create a new input file that is used
to run a self-consistent (scf) cycle and a non self-consistent (nscf) cycle using the charge density calculated on the previous run.

2. GW convergence

Afterwards you can run a GW calculation using the `gw_si.py` script and a Bethe-Salpether (BSE) calculation using the `bse_si.py`.
In the begining of each script (for GW or BSE) there is a check for the presence of the SAVE database. In case it is not present it will be generated.

In the `gw_conv_si.py` you will find an example of how ot use the `optimize()` function to converge the calculation parameters.

.. code-block:: python

    #create the yambo input file
    y = YamboIn('%s -d -g n -V all'%yambo,folder='gw_conv')
    y['QPkrange'][0][2:4] = [6,10]
    conv = { 'FFTGvecs': [[10,15,20],'Ry'],
             'NGsBlkXd': [[1,2,5], 'Ry'],
             'BndsRnXd': [[1,10],[1,20],[1,30]] }

    def run(filename):
        """ Function to be called by the optimize function """
        folder = filename.split('.')[0]
        print(filename,folder)
        os.system('cd gw_conv; yambo -F %s -J %s -C %s 2> %s.log'%(filename,folder,folder,folder))

    y.optimize(conv,run=run)

This code will run yambo as many times as variables specified in the `conv` dictionary.
The first calculation is always called `reference` and takes always the first element of each of the lists.
Then for each element of the list that is not the first one a calculation is made.

3. Collect the data

Once all the calculations are finished it's time to pack all the files in the `json` format for posterior analysis.
For this use the `YamboOut()` class:

.. code-block:: python

  #pack the files in .json files
  for dirpath,dirnames,filenames in os.walk('gw_conv'):
    #check if there are some output files in the folder
    if ([ f for f in filenames if 'o-' in f ]):
        y = YamboOut(dirpath)
        y.pack()


4. Plot the data

After this you should have a set of `.json` files in the folder, one for each calculation.
To make a plot of them all you just need to run:

.. code-block:: python

  #plot the results using yambmo analyser
  y = YamboAnalyser('gw_conv')
  y.plot_gw('qp')

You can add more plots by simply adding more files in the folder you give as input to the `YamboAnalyser()` class.

BN
----
In a similar fashion to the previous example you can run the ground state calculations for Boron Nitride using `gs_bn.py`.

.. code-block:: bash

    python gs_bn.py

MoS\ :sub:`2`
-----------------

The groundstate calculation for MoS\ :sub:`2` is made in a similar fashion as the previous examples.
If some of the steps are already calculated you can tell the script not to run them using for example:

.. code-block:: bash

    python gs_mos2.py -n2

The option `-n2` will tell the script not to run the double grid `nscf` calculation.
