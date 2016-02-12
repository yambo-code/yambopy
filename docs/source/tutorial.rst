Tutorial
==========

In the `tutorial` folder you can find some examples of how to get started using `yambopy`.
The first step in any calculation with `yambo` is to calculate the ground state proprieties using either `abinit` or `pw.x`.
We don't have support to read and write `abinit` input files. To do that you should use the `abipy <https://github.com/gmatteo/abipy>`_ package.
Included in the `yambopy` package we include some basic scripts to generate Quantum Espresso input files.

GW convergence (Si)
--------------------

**1. Ground State**

Go to the `tutorial` folder and run the ground state calculation using the `gs_si.py` file:

.. code-block:: bash

    python gs_si.py

The script will run a relaxation of the structure, read the optimized cell parameter and create a new input file that is used
to run a self-consistent (scf) cycle and a non self-consistent (nscf) cycle using the charge density calculated on the previous run.

**2. GW convergence**

Afterwards you can run a GW calculation using the `gw_si.py` script and a Bethe-Salpether (BSE) calculation using the `bse_si.py`.
In the beginning of each script (for GW or BSE) there is a check for the presence of the SAVE database. In case it is not present it will be generated.

In the `gw_conv_si.py` you will find an example of how to use the `optimize()` function to converge the calculation parameters.

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

This code will run `yambo` as many times as variables specified in the `conv` dictionary.
The first calculation is called `reference` and uses as parameters the first element of each of the lists.
For each of the other elements of the list a calculation is made.

**3. Collect the data**

Once all the calculations are finished it's time to pack all the files in the `json` format for posterior analysis.
For this use the `YamboOut()` class:

.. code-block:: python

  #pack the files in .json files
  for dirpath,dirnames,filenames in os.walk('gw_conv'):
    #check if there are some output files in the folder
    if ([ f for f in filenames if 'o-' in f ]):
        y = YamboOut(dirpath,save_folder=dirpath)
        y.pack()

This snippet of code can be called using the function:

.. code-block:: python

    pack_files_in_folder('gw_conv',save_folder='gw_conv')

**4. Plot the data**

After this you should have a set of `json` files in the folder, one for each calculation.
To make a plot of them all you just need to run:

.. code-block:: python

  #plot the results using yambo analyser
  y = YamboAnalyser('gw_conv')
  y.plot_gw('qp')
  y.plot_gw_path('qp')

You can add more plots by simply adding more files in the folder you give as input to the `YamboAnalyser()` class.
At the end you should obtain a plot like this:

.. image:: figures/gw_si.png

Coulomb-cutoff (BN)
-------------------------------

In this example we will test the convergence of the coulomb truncation for a BSE calculation in single layer Boron Nitride.
For that we define a loop where we perform a self-consistent ground state calculation, non self-consistent calculation, create the databases
and run `yambo` with increasing vacuum and plot the absorption spectra.

**2. Coulomb truncation convergence**

In the folder `tutorials/bn/` you find the python script `bse_cutoff.py`.
You can run this script with:

.. code-block:: bash

    python bse_cutoff.py -r

The main loop changes the `layer_separation` variable using values from a list.
In the script you can find how the functions `scf`, `ncf` and `database` are defined.

.. code-block:: python

    #for each separation run the ground state calculation and
    for layer_separation in layer_separations:

      root_folder = "%s/%d"%(work_folder,layer_separation)
      if not os.path.isdir(root_folder):
          os.makedirs(root_folder)

      # run the ground state calculation
      print("scf cycle")
      scf(layer_separation,folder="%s/scf"%root_folder)
      os.system("cd %s/scf; pw.x < %s.scf > scf.log"%(root_folder,prefix))

      # run the non self consistent calculation
      print("nscf cycle")
      src ='%s/scf/%s.save'%(root_folder,prefix)
      dst ='%s/nscf/%s.save'%(root_folder,prefix)
      nscf(layer_separation,folder="%s/nscf"%root_folder)
      os.system( 'cp -r %s %s'%(src,dst) )
      os.system("cd %s/nscf; pw.x < %s.nscf > nscf.log"%(root_folder,prefix))

      # generate the database
      database('%s'%root_folder,nscf_folder="%s/nscf"%root_folder)

      # calculate the absorption spectra using yambo
      y = YamboIn('yambo -r -b -o b -k sex -y d -V all',folder=root_folder)

      y['FFTGvecs'] = [30,'Ry']
      y['NGsBlkXs'] = [1,'Ry']
      y['BndsRnXs'] = [1,30]

      y['CUTGeo'] = 'box z'
      y['CUTBox'] = [0,0,layer_separation-1]

      y['KfnQP_E']  = [1.0,1.0,1.0] #scissor operator
      y['BSEBands'] = [3,6]
      y['BEnSteps'] = 500
      y['BEnRange'] = [[1.0,6.0],'eV']
      y.write('%s/yambo_run.in'%root_folder)
      os.system('cd %s; %s -F yambo_run.in -J %d'%(root_folder,yambo,layer_separation))

**3. Plot the convergence**

You can plot the results using:

.. code-block:: bash

    python bse_cutoff.py -p

You should obtain a plot like this:

.. image:: figures/bse_cutoff.png


Parallel Bethe-Salpeter (MoS\ :sub:`2`)
-----------------------------------------------------------------

In this tutorial we will show how you can paralellize the dielectric function calculation in
separate jobs for a BSE optical absorption spectra calculation.

The idea is that in certain clusters its advantageous to split the dielectric function calculation
in smaller jobs (one for each q-point) that can run at the same time.
Using the `yambo` parallelization you can separate the dielectric function calculation among many cpus
using the variable `q` in `X_all_q_CPU` and `X_all_q_ROLEs`. The issue is that you still need to make a big reservation
and in some cases there is load imbalance (some nodes end up waiting for others). Splitting in smaller jobs
can help your jobs to get ahead in the queue and avoid the load imbalance.
If there are many free nodes you might end up running all the q-points at the same time.

**The idea is quite simple:** you create an individual input file for each q-point, submit each job separately, collect
the results and do the final BSE step (this method should also apply for a GW calculation).

**1. Ground State**

The ground state calculation for MoS\ :sub:`2` is made in a similar fashion as the previous examples.
If some of the steps are already calculated you can tell the script not to run them using for example:

.. code-block:: bash

    python gs_mos2.py -n2

The option `-n2` will tell the script not to run the double grid `nscf` calculation.

**2. Parallel Dielectric function**

Here we tell `yambo` to calculate the dielectric function. We read the number of q-points the system has
and generate one input file per q-point. Next we tell `yambo` to calculate the first q-point. `Yambo` will calculate the dipoles
and the dielectric function at the first q-point.
Once the calculation is done we copy the dipoles to the SAVE directory. After that we can run each q-point calculation
as a separate job.
Here the user can decide to submit one job per q-point on a cluster or use a program like
`gnuparallel <http://www.gnu.org/software/parallel/>`_ to schedule the jobs according to the
available resources. In this example we use the second option.

.. code-block:: python

    #create the yambo input file
    y = YamboIn('yambo -r -b -o b -V all',folder='bse_par')
    y['FFTGvecs'] = [15,'Ry']
    y['NGsBlkXs'] = [1,'Ry']
    y['BndsRnXs'] = [[1,40],'']
    y.write('bse_par/yambo_run.in')

    #get the number of q-points
    _,nkpoints = y['QpntsRXs'][0]

    #prepare the q-points input files
    f = open('jobs.sh','w')
    for nk in xrange(1,int(nkpoints)+1):
        y['QpntsRXs'] = [[nk,nk],'']
        y.write('bse_par/yambo_q%d.in'%(nk))
        if nk != 1:
            f.write('cd bse_par; %s -F yambo_q%d.in -J %d\n'%(yambo,nk,nk))
    f.close()

    #calculate first q-point and dipoles
    os.system('cd bse_par; %s -F yambo_q1.in -J 1'%yambo)
    #copy dipoles to save
    os.system('cp bse_par/1/ndb.dip* bse_par/SAVE')
    #run jobs using gnuparallel
    os.system('parallel :::: jobs.sh')

**3. BSE**

Once the dielectric function is calculated its time to collect the data in one folder and
do the last step of the calculation: solve the Bethe-Salpeter equation.
The `merge_eps.py` gathers the dielectric function databases into one folder with the correct numeration.
In the next version of `yambo` this step won\'t be necessary (yambo will write the databases with the correct numeration).
The the last step is to generate the BSE Hamiltonian and diagonalize it.

.. code-block:: python

    #gather all the files
    os.system('cp merge_eps.py bse_par')
    os.system('cd bse_par; python merge_eps.py')

    y = YamboIn('yambo -r -b -o b -k sex -y d -V all',folder='bse_par')
    y['FFTGvecs'] = [15,'Ry']
    y['NGsBlkXs'] = [1,'Ry']
    y['BndsRnXs'] = [[1,40],'']
    y['BSEBands'] = [8,11]
    y['BEnSteps'] = 500
    y['BEnRange'] = [[1.0,6.0],'eV']
    y.arguments.append('WRbsWF')

    y.write('bse_par/yambo_run.in')
    os.system('cd bse_par; %s -F yambo_run.in -J yambo'%yambo)


**3. Collect and plot the results**

You can plot the data much in the same way as you did for the GW calculation.

.. code-block:: python

    #collect the data
    pack_files_in_folder('bse_par')

    #plot the results using yambo analyser
    y = YamboAnalyser('bse_par')
    y.plot_bse('eps')

You should obtain a plot like this:

.. image:: figures/bse_mos2.png
