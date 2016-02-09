Quickstart
==========

Only the `yambo` runlevels are hardcoded in `yambopy`. This means that any new variable in `yambo` can me immediately used from `yambopy` without needing to add new variables to the `yambopy` python code.

YamboIn
-------

Read, write, create and manipulate yambo input files with python.

The class can be initialized in two ways:
    - Specifying the runlevels: `yambopy` will run yambo, read the generated input file and initialize the class with the variables
    - Without runlevels: just initialize the class and then the user can specify the variables in the script

.. code-block:: python

    from yambopy import YamboIn

    yi = YamboIn('-b -o b -k sex -y d',folder='tutorial')

    #set some variables
    yi['FFTGvecs'] = [15,'Ry']
    yi['NGsBlkXs'] = [1,'Ry']
    yi['BndsRnXs'] = [[1,30],'']

    #print the input file in the terminal
    print yi

    #this will write the file on the hard drive
    yi.write('input.in')

Keep in mind that the runlevels are hardcoded in `yambopy` this means that each new runlevel has to be added to the `_runlevels` list in the YamboIn class.

An interesting feature of the `YamboIn()` class is the `optimize()` funciton that helps to make convergence tests.
After you defined the all the relevant variables for the yambo input file you might want to run the code several times changing some of them to see how the final result changes.
Using the initialization above you can use the following code to see how the results change as you use more sctrict convergence paramaters:

.. code-block:: python

  conv = { 'FFTGvecs': [[10,15,20],'Ry'],
           'NGsBlkXs': [[5,10,20], 'Ry'],
           'BndsRnXs': [[1,10],[1,20],[1,30]] }

  def run(filename):
      folder = filename.split('.')[0]
      os.system('yambo -F %s -J %s'%(filename,folder))

  yi.optimize(conv,run=run)

YamboOut
--------

This class is used to read the output files of a typical yambo calculation and pack the results in a `.json` file for posterior analysis using `YamboAnalyser`.
Currently we save the `o-` data and the input file that is written at the end. YamboOut also tries to get information about the positions of the atoms and lattice from the `SAVE` directory.
This is only possible if you have netCDF4 support in your python installation.

.. code-block:: python

    from yambopy import YamboOut

    yo = YamboOut('tutorial')
    yo.pack()

    #print the output file in the terminal
    print yo

YamboAnalyser
-------------

This class is used to read the `.json` files generated with YamboOut and plot them.

.. code-block:: python

    from yambopy import YamboAnalyser

    ya = YamboAnalyser('tutorial')
    ya.plot_bse() #for the case of a bse calculation
    ya.plot_gw() #for the case of a gw calculation

    #print the output file in the terminal
    print yo


pw.x
-----

`yambopy` provides a class `PwIn()` to create and edit input files for `pw.x` from the `Quantum Espresso <http://www.quantum-espresso.org/>`_ suite.
This class works in a similar way as `YamboIn()` so you can start it either by reading a file from the hard drive
or specifying the variables in a python script.

The `input <http://www.quantum-espresso.org/wp-content/uploads/Doc/INPUT_PW.html>`_ file for `pw.x` is split into different sections.
youc an acess the variables for each section using :code:`.<section>['variable_name']`.

Here is an example of how to create an input file for Silicon.

.. code-block:: python

    from qepy import PwIn

    qe = PwIn()
    qe.atoms = [['Si',[0.125,0.125,0.125]],
                ['Si',[-.125,-.125,-.125]]]
    qe.atypes = {'Si': [28.086,"Si.pbe-mt_fhi.UPF"]}

    qe.control['prefix'] = "'si'"
    qe.control['wf_collect'] = '.true.'
    qe.system['celldm(1)'] = 10.3
    qe.system['ecutwfc'] = 60
    qe.system['occupations'] = "'fixed'"
    qe.system['nat'] = 2
    qe.system['ntyp'] = 1
    qe.system['ibrav'] = 2
    qe.kpoints = [4, 4, 4]
    qe.electrons['conv_thr'] = 1e-8

    #print the output file in the terminal
    print qe

    #write the input file
    qe.write('qe.in')


ph.x
-----

`yambopy` provides a class `PhIn()` to write input files for `ph.x` from the  `Quantum Espresso <http://www.quantum-espresso.org/>`_ suite.

.. code-block:: python

    from qepy import PhIn

    ph = PhIn()
    ph['nq1'],ph['nq2'],ph['nq3'] = [1,1,1]
    ph['tr2_ph'] = 1e-12
    ph['prefix'] = "'si'"
    ph['epsil'] = ".false."
    ph['trans'] = ".true."
    ph['fildyn'] = "'si.dyn'"
    ph['fildrho'] = "'si.drho'"
    ph['ldisp'] = ".true."

    print ph
    ph.write('si.ph')

dynmat.x
--------

`yambopy` provides a class `DynmatIn()` to write input files for `dynmat.x` from the  `Quantum Espresso <http://www.quantum-espresso.org/>`_ suite.

.. code-block:: python

    from qepy import DynmatIn

    md = DynmatIn()
    md['asr'] = "'simple'"
    md['fildyn'] = "'si.dyn1'"
    md['filout'] = "'si.modes'"

    print md
    md.write('si.dynmat'%folder)
