The yambopy databases (DBs) are a set of classes to interact with the modules of 
Yambo. The classes read the netCDF files that contain dipoles, exciton states, quasi-particles states, 
dielectric function, Green's functions, etc. They are useful to access and plot
results, check intermediate results. They are located in the folder ``yambopy/dbs``.

YamboExcitonDB
~~~~~~~~~~~~~~

Located in ``yambopy/dbs/excitondb.py``, 
read the excitonic states database from yambo. It is useful to plot excitonic weigths on
top of the band structure or in a map of the Brillouin zone.

There is a short example: ``tutorials/bn/plot-excitondb.py``. Previously one needs to
run a Bethe-Salpeter calculation with the option diagonalization and with the flag
``WRbsWF``.
We have defined the common path along the Brillouin zone for hexagonal lattices:

.. code-block:: bash
   path = Path([ [[  0.0,  0.0,  0.0],'$\Gamma$'],
                 [[  0.5,  0.0,  0.0],'M'],
                 [[1./3.,1./3.,  0.0],'K'],
                 [[  0.0,  0.0,  0.0],'$\Gamma$']], [int(npoints*2),int(npoints),int(sqrt(5)*npoints)] )

We have selected the ground state excitonic state. In order to read and plot the excitonic state we also need to charge
the information of the structure:

.. code-block:: bash
   save = YamboSaveDB.from_db_file(folder='bse_flow/t0/SAVE')
   lat  = YamboLatticeDB.from_db_file(filename='bse_flow/t0/SAVE/ns.db1')
   yexc = YamboExcitonDB.from_db_file(lat,filename='ndb.BS_diago_Q01',folder='bse_flow/t0/run')

In order to plot the bands without interpolation we select the function ``get_exciton_bs``

.. code-block:: bash
   exc_bands = yexc.get_exciton_bs(save,path,states,size=1.0)

Usually k-grids of Bethe-Salpeter calculation are not enough dense to obtain a smooth band structure. There is the option
of performing an interpolation using the function ``interpolate``.

.. code-block:: bash
   exc_bands_inter = yexc.interpolate(save,path,states,lpratio=5,f=None,size=0.5,verbose=True)

.. figure:: figures/exciton-band-not-interpolated.png
    :width: 200px
    :align: center
    :height: 100px
    :alt: alternate text
    :figclass: align-center

.. image:: figures/exciton-band-interpolated.png
   :width: 50

YamboQPDB
~~~~~~~~~

Located in ``yambopy/dbs/qpdb.py``, this class reads quasi-parcticle data files
generated by Yambo ``ndb.QP``. These files describe the quasiparticle states,
such as the quasi-particle energies, the lifetimes and the Z factors. There is an
example available in ``tutorials/bn/plot-qp``. A run of GW is needed before running
the script.

The class ``YamboQPDB``

.. image:: figures/gw-scissor.png
   :width: 3%

.. image:: figures/gw-bands-not-interpolated.png
   :width: 3%

.. image:: figures/gw-bands-interpolated.png
   :width: 3%


YamboSaveDB
~~~~~~~~~~~

Reads the information from the SAVE database in Yambo. The arguments are:

.. code-block:: bash
   ``save``: Path with the save folder (default:SAVE)
   ``filename``: name of the filename of the ns.db1 database created with yambo (default:ns.db1)

YamboLatticeDB
~~~~~~~~~~~~~~

Class to read the lattice information from the netcdf file ``ns.db1``.
