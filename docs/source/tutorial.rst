Tutorial
==========

In the `tutorial` folder you can find some examples fo how to get started using `yambopy`.
The first step in any calculation is to calgulate the ground state properties using either `abinit` or `pw.x`.
All the examples presented here use `pw.x` to calculate the groudstate.

Si
----
To calculate the ground state run the `gs_si.py` file:

.. code-block:: bash

    python gs_si.py

The script will run a relaxation of the structure, read the optimized cell parameter and create a new inpuf file that is used
to run a Self-consistent (scf) cycle and a Non Self-consistent (nscf) cycle using the charge density calculated on the previous run.

Aftewards you can run a GW calculation using the `gw_si.py` script, a Bethe-Salpether (BSE) calculation using the `bse_si.py`.
In the begining of each script (for GW or BSE) there is a check for the presence of the SAVE database. In case it is not present it will be generated.

If the `gw_conv_si.py` you will find an example of how ot use the `optimize()` function to converge the calculation parameters.

BN
----
In a similar fashion to the previous example you can run the ground state calculations for Boron Nitride using `gs_bn.py`.

.. code-block:: bash

    python gs_bn.py

MoS\ :sub:`2`
-----------------

The groundstate calculation is made in a similr fashion as the previous examples.
If you some of the steps is already calculated you can tell the script not to run them using for example:

.. code-block:: bash

    python gs_mos2.py -n2

The option `-n2` will tell the script not to run the double grid `nscf` calculation.
