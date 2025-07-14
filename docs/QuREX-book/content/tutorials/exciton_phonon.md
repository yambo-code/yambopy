# Exciton-Phonon
Exciton-phonon coupling refers to the interaction between excitons (bound states of an electron and a hole) and phonons (quanta of lattice vibrations) in a material. This coupling can significantly affect the optical and electronic properties of the material. When an exciton interacts with phonons, it can scatter, leading to changes in its energy and momentum. This process is crucial for understanding phenomena such as thermalization, relaxation, and recombination of excitons in semiconductors and other materials.
This workflow is split up into several sections:


0. [DFT](#step-0-dft-ground-state-calculation-using-pwx)
Calculate ground state wavefunctions of the system using Quantum Espresso `pw.x`.
1. [DFPT](#step-1-dfpt-calculation-using-phx)
Calculate the phonons and $\Delta V_{scf}$ for the electron-phonon matrix elements using Quantum Espresso `ph.x`.
2. [ELPH](#step-2-dvscf-calculation-using-letzelph-yambopy)
Compute the electron-phonon matrix elements $g_{kkp}$ using the [LetzElPhC code](https://github.com/muralidhar-nalabothula/LetzElPhC/). 
3. [GW-BSE](#step-3-gw-bse-calculation-using-yambo)
Compute BSE kernel and diagonalize on top of GW calculation, using the Yambo code to get the exciton wavefunctions.
4. [EXC-PH](#step-4-exciton-phonon-calculation-using-yambopy)
Combine all of the above to compute the exciton-phonon matrix elements that are used to compute photoluminesensce spectrum with python [scripts](https://github.com/muralidhar-nalabothula/PhdScripts/).

We will go through these steps with a simple system of hBN monolayer that can be readily run with relatively small computational cost. Step 1 of calculating the phonons takes some time, but one can parallelize and calculate at the same time the GW-BSE kernel at step 3.
For more detailed information on exciton-phonon coupling and its effects on luminescence in hexagonal boron nitride, you can refer to the paper: {cite}`PhysRevMaterials.7.024006`.


## Step 0: DFT ground state calculation using pwx
Details about the DFT ground state calculation using pw.x.
`force_symmorphic=.true.`
The DFT ground state calculations serve 2 purposes, to acquire the starting wave functions for the calculations of the phonons, but also to acquire the Kohn-Sham energies that are used for the MBT calculations to get the exciton wavefunctions. 



## Step 1: DFPT calculation using ph.x
Details about the DFPT calculation using ph.x.

## Step 2: DVSCF calculation using Letzelph-yambopy
Details about the DVSCF calculation using LetzElPhC-yambopy.
For documentation on how to compile and run the LetzELPhC code using yambopy have a look at the [LetzElPhC documentation](../external_files/LetzElPhC_documentation.pdf).
We run `yambopy l2y`:
```
yambopy l2y -ph /path/of/ph_input.in -b n_i n_f -par n_qpools n_kpools -lelphc /path/to/lelphc_exe -D
```
Keep the `-D` flag to make sure we keep the required `ndb.elph` and `ndb.Dmat` which are the electron-phonon matrix elements, and the Dynamical matrices respectively. 

We should now have the electron-phonon database `ndb.elph`. 

## Step 3: GW-BSE calculation using Yambo
Details about the GW-BSE calculation using Yambo.
Construct the BSE kernel and diagonalize to get the exciton wavefunctions.
Make sure to turn on the flag: `WRbsWF` to write the exciton wavefunctions, which is commented out by default.
For the theoretical background on this section please look at [GW](../theory/GW) and [BSE](../theory/bse_equation).

## Step 4: Exciton-Phonon calculation using yambopy
Details about the Exciton-Phonon calculation using yambopy.
For this step we use these [scripts](https://github.com/muralidhar-nalabothula/PhdScripts/) to calculate the exciton-phonon matrix elements but also to calculate the photoluminescense spectrum of the material. The script we will be using is called `ex_ph_program.py`.

It takes a few basic inputs:
```
- calc_folder: /Path/ #where the calculations took place.
- SAVE_dir: calc_folder + '/path/to/SAVE'
- BSE_dir: calc_folder + '/path/of/job/bse/'
- elph_file: calc_folder + '/path/to/ndb.elph'
- Dmat_file: calc_folder + '/path/to/ndb.Dmats'
- nstates: 19*4             # Number of states to include in the PL, `nq * ntransitions` (e.g., 19*4)
- lumin: True               # Compute luminescence if set to `True`
- Exph: True 
- Temp: 20                  # Temperature used in luminescence (in Kelvin, e.g., 20)
- ome_range:[1, 8, 1000]    # Range for omega in the format `(min, max, numpoints)` (in eV, e.g., [1, 8, 1000])
- broading: 0.005           # Broadening (in eV, e.g., 0.005)
- npol: 2                   # Polarization, set to 2 for 2D materials
```
Run the python script.

Now you obtain `Ex-ph.npy` and `luminescence_intensities.dat`. Which we can plot:
![Photoluminescence Spectrum of 2D hBN](/plots/2D_hBN_PL.png)

Follow the [Jupyter Notebook](../../notebooks/exciton_phonon.ipynb) for a pratical tutorial on how to compute photoluminescence and Raman including exciton-phonon coupling.


# References

```{bibliography}
    :cited:
```