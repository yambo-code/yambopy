Flows
=====

The flows structure (or function) take care of handling the tasks. The tasks are
interdependent works. For example, the calculation of the BSE spectra, the calculations of the ground-state and band structure, etc. 
We have created already some flows with common yambo calculations, named as ``Yambopy Factories`` and hosted in ``yambopy/io/factories.py``.

YamboTask
~~~~~~~~~~~~~~~~~

This is one of the basic and main tasks
The basic tasks are one-shot calculations using Yambo. They are the building-blocks of combined tasks. 

BSE Task
--------

We have created the example ``flow-bse.py`` in the silicon folders to demonstrate how to create a task. The same
can be used for any other run level of Yambo such as GW run levels, non-linear run levels or real-time run levels.

The variables are set in dictionaries like ``yambo_dict`` and we create a list of task. The first task is to set
the location of the SAVE folder.

.. code-block:: bash

    p2y_task = P2yTask.from_folder('nscf_flow/t2')

We define the usual group of variables using dictionaries:

.. code-block:: bash

    # Coulomb-cutoff and RIM dictionary
    cutoffdict = dict(RandQpts=1000000,RandGvec=[1,'RL'])

    # Parallel Environment dictionary
    paradict = dict(X_all_q_ROLEs="q",X_all_q_CPU="2")

    # BSE variables dictionary
    bse_dict = dict(BEnSteps=1000,  FFTGvecs=[10,'Ry'], BEnRange=[[0,5],'eV'], BndsRnXp=[1,10],
                     NGsBlkXp=[1,'Ry'], BSENGexx=[10,'Ry'], BSENGBlk=[1,'Ry'], BSEBands=[2,7])
                                                                                                                
    # Merge all dict variables
    yamboin_dict = {**yamboin_dict,**cutoffdict,**paradict}

Once we have all variables we can define the BSE task (option ``from_runlevel``)

.. code-block:: bash

    bse_task = YamboTask.from_runlevel([p2y_task],'-r -o b -b -k sex -y h -V all',yamboin_dict)

Once we have all the tasks defined we create a list of task:

.. code-block:: bash

    tasks.append(bse_task)

Now the list of tasks defines the Yambopy Flow:

.. code-block:: bash

    yambo_flow = YambopyFlow.from_tasks('bse_flow',tasks)

And we can create and run the flow.

.. code-block:: bash

    yambo_flow.create(agressive=True)
    yambo_flow.run()

If all was done correctly, running the example:

.. code-block:: bash
    python flow_bse.py

We will obtain the following message:

.. code-block:: bash

   ======================YambopyFlow.run=======================
   t0  YamboTask  ready
   ========================YambopyFlow=========================
   t0  YamboTask  done

Note that by default we obtain the results in the folder ``bse_flow/t0`` with the jobname ``run``. We have only set one
task and the corresponding folder is ``t0``. In the situation of multiple tasks the results will be separated
according to the task order.


Yambopy Factories
~~~~~~~~~~~~~~~

The ``factories`` are usually frequent interdependent Yambo tasks. For example, we have created some interdependent 
Yambo tasks like convergence tests, QP+BSE calculations.

PwNscfYamboIPChiTasks
---------------------

YamboIPChiTask
---------------

This factory run the calculation of the dielectric function at the independent-particle 
approximation.

YamboQPTask
-----------

This factory run a GW calculation.

YamboQPBSETasks
---------------

This factory run a GW and Bethe-Salpeter calculation.

Quantum Espresso Factories
~~~~~~~~~~~~~~~~~~~~~~~~~~

In addition to the Yambo-related tasks, yambopy has also pw-related tasks to perform self-consistent non-selfconsistent calculations, band structure calculations and cell optimization.


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

PhPhononTasks
------------

ABINIT Factories
~~~~~~~~~~~~~~~~~~~~~~~~~~

AbinitNscfTasks
---------------

AbinitNscfTasksFromAbinitInput
---------------

