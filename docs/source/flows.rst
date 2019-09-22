Flows
=====

The flows structure (or function) take care of handling the tasks.


Yambo Tasks
~~~~~~~~~~~

Tasks contain interdependent works.

YamboQPBSETasks
---------------

This task run a GW and Bethe-Salpeter calculation.

YamboQPTask
-----------



Quantum Espresso Tasks
~~~~~~~~~~~~~~~~~~~~~~

In addition to the Yambo-related tasks, yambopy has also pw-related tasks to perform self-consistent
non-selfconsistent calculations, band structure calculations and cell optimization.


PwRelaxTasks
------------

The relaxation task performs three concatenated calculations. First, the atomic relaxation is performed. The second calculation reads the
new atomic positions and it performs a cell relaxation. The third and last calculation is the self-consistent calculations of the density
with the optimized cell parameters and atomic positions.

atomic relaxation >> cell relaxation >> self-consistent calculation

This flow includes some specific variables as inputs:

.. code-block:: bash

    cell_dofree
    pseudo_dir
    spinor 
    pseudo_dir


You can find examples for silicon and hexagonal BN in the folder ``tutorials/si`` and ``tutorials/bn``, respectively. The example runs with the following command:

.. code-block:: bash
    
    python flow-pw.py -r


PwNscfTasks
-----------

The Nscf task performs a self-consistent and a non-self consistent calculation plus the ``p2y`` runs to prepare the QE output file in the Yambo format . This is the preliminar calculation before using Yambo.

This flow includes some specific variables as inputs:

.. code-block:: bash

    nscf_bands
    nscf_kpoints
    spinor
    pseudo_dir

   
You can find examples for silicon and hexagonal BN in the folder ``tutorials/si`` and ``tutorials/bn``, respectively. The example runs with the following command:

.. code-block:: bash
    
    python flow-pw.py -n

PwBandsTasks
------------

This taks performs a self-consisten and a band calcualtion using QE. The options are similar to the options of PwNscfTasks with the exception of the variable ``path_kpoints``. This variable is defined using the class ``Path``. In the tutorial for silicon we have defined the path as follows:

.. code-block:: bash

    p = Path([ [[1.0,1.0,1.0],'$\Gamma$'],
               [[0.0,0.5,0.5],'$X$'],
               [[0.0,0.0,0.0],'$\Gamma$'],
               [[0.5,0.0,0.0],'$L$']], [20,20,20])

The example runs with the command:

.. code-block:: bash
    
    python flow-pw.py -b
    
Optionally is possible to plot the band structure using the class ``PwXML``:
 
.. code-block:: bash
    
    python flow-pw.py -p
