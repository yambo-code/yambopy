yambopy
==========

Only the `yambo` run levels are hardcoded in `yambopy`. This means that any new
variable in `yambo` can me immediately used from `yambopy` without needing to
add new variables to the `yambopy` python code.

YamboIn
~~~~~~~~~~

Read, write, create and manipulate `yambo` input files with python.

The class can be initialized in two ways:
    - Specifying the run levels: `yambopy` will run `yambo`, read the generated input file and initialize the class with the variables
    - Without run levels: just initialize the class and then the user can specify the variables in the script

.. code-block:: python

    from yambopy import YamboIn

    yi = YamboIn('-b -o b -k sex -y d',folder='tutorial')

    #set some variables
    yi['FFTGvecs'] = [15,'Ry']
    yi['NGsBlkXs'] = [1,'Ry']
    yi['BndsRnXs'] = [1,30]

    #print the input file in the terminal
    print yi

    #this will write the file on the hard drive
    yi.write('input.in')

Keep in mind that the run levels are hardcoded in `yambopy` this means that each
new run level has to be added to the `_runlevels` list in the `YamboIn()` class.

An interesting feature of the `YamboIn()` class is the `optimize()` function that
helps to make convergence tests. After you defined the all the relevant variables
for the `yambo` input file you might want to run the code several times changing
some of them to see how the final result changes.
Using the initialization above you can use the following code to see how the
results change as you use more strict convergence parameters:

.. code-block:: python

  conv = { 'FFTGvecs': [[10,15,20],'Ry'],
           'NGsBlkXs': [[5,10,20], 'Ry'],
           'BndsRnXs': [[1,10],[1,20],[1,30]] }

  def run(filename):
      folder = filename.split('.')[0]
      os.system('yambo -F %s -J %s'%(filename,folder))

  yi.optimize(conv,run=run)

YamboOut
~~~~~~~~~

This class is used to read the output files of a typical `yambo` calculation and
pack the results in a `.json` file for posterior analysis using `YamboAnalyser`.
Currently we save the `o-` data and the input file that is written at the end.
`YamboOut()` also tries to get information about the positions of the atoms and
lattice from the `SAVE` directory.
This is only possible if you have netCDF4 support in your python installation.

.. code-block:: python

    from yambopy import YamboOut

    yo = YamboOut('tutorial')
    yo.pack()

    #print the output file in the terminal
    print yo

YamboAnalyser
~~~~~~~~~~~~~~~~~~

This class is used to read the `.json` files generated with `YamboOut()` and plot them.

.. code-block:: python

    from yambopy import YamboAnalyser

    ya = YamboAnalyser('tutorial')
    ya.plot_bse() #for the case of a bse calculation
    ya.plot_gw() #for the case of a gw calculation

    #print loaded files in the terminal
    print ya

