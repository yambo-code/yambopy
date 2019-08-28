qepy
==========

PwIn
~~~~~~~~~~~~~~~~~~

`yambopy` provides a class `PwIn()` to create and edit input files for `pw.x`
from the `Quantum Espresso <http://www.quantum-espresso.org/>`_ suite.
This class works in a similar way as `YamboIn()` so you can start it either by
reading a file from the hard drive
or specifying the variables in a python script.

The `input <http://www.quantum-espresso.org/wp-content/uploads/Doc/INPUT_PW.html>`_
file for `pw.x` is split into different sections.
you can access the variables for each section using :code:`.<section>['variable_name']`.

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


PhIn
~~~~~~~~~

`yambopy` provides a class `PhIn()` to write input files for `ph.x` from the
`Quantum Espresso <http://www.quantum-espresso.org/>`_ suite.

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

DynmatIn
~~~~~~~~~~~~~

`yambopy` provides a class `DynmatIn()` to write input files for `dynmat.x`
from the  `Quantum Espresso <http://www.quantum-espresso.org/>`_ suite.

.. code-block:: python

    from qepy import DynmatIn

    md = DynmatIn()
    md['asr'] = "'simple'"
    md['fildyn'] = "'si.dyn1'"
    md['filout'] = "'si.modes'"

    #write the input file in the terminal
    print md
    md.write('si.dynmat'%folder)

Unfolding
~~~~~~~~~~~~~
The class `Unfolding()` is useful to unfold the electronic structure calculated in a supercell into the original primitive cell of
the material. Currently it generates and reads Quantum Espresso XML files. The class is based in the work of Popescu and Zunger published in `Phys. Rev. B 85, 085201 (2012) <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.85.085201>`_.

There is an example adapted to hBN tutorial/bn-folding. Currently there are several additional options:

- write to file: If True it prints in the file projection.dat the results of the unfolding.
- spin: "none" or "spinor".
- band_min: To avoid the processing in core levels we can set the starting band for the unfolding.




