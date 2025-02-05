# LiF tutorial - From ground state DFT to excited state MBPT via Wannerization methods
This tutorial guides you through calculating the electronic ground state properties of Lithium Fluoride (LiF) using Density Functional Theory (DFT)[^1]. For our simulations, we'll be utilizing the Quantum Espresso code package[^6]. We'll further enhance the electronic band structure analysis by employing a Wannierization procedure[^2] via Fourier interpolation. This step aims to increase the sampling grid while minimizing computational expenses.

## Advanced analysis with Yambo

Upon establishing the ground state properties, our attention shifts to excited state phenomena using Yambo[^4][^5]. Our focus here extends to assessing the absorption spectra and delineating the exciton band structure. The Wannierization procedure applied to the excitonic Hamiltonian mirrors that of the electronic scenario, providing a basis for comparison with existing literature findings[^3].

## Tools and libraries
For the effective pre- and post-processing of our data, we rely on [yambopy](https://github.com/rreho/yambopy/tree/devel-tbwannier), a robust tool available at yambopy GitHub repository. The specific branch `devel-tbwannier` is dedicated to the input/output (IO) and interpolation tasks related to both the electronic and excitonic Hamiltonians.

## Computational details
All simulations within this tutorial adhere to the use of LDA pseudopotentials, intentionally excluding Spin-Orbit coupling (SOC) to streamline our computational approach.

## Crystal structure LiF
Lithium Fluoride (LiF) is an insualtor with a wide band-gap and the following characteristics:
- FCC lattice with 48 symmetry operations
- Two atoms per cell, Li and F (8 electrons)
- Lattice constant 7.61 [a.u.]
- Plane waves cutoff 80 Ry

For the Yambo analysis we will see that any independent particle approximation drastically fails to describe the complex absorption structures that are observed experimentally.
Indeed, the spectrum of bulk LiF is dominated by a strongly bound exciton (3 eV binding energy)

<div style="text-align: center;">
    <img src="https://hackmd.io/_uploads/BJrSA7ZxA.png" alt="LiF" width="70%">
</div>

## DFT
This step is foundational, setting the stage for all subsequent analyses by establishing a precise understanding of the system's electronic distribution.
### self-consistent field calculation
The initial step involves calculating the electronic charge density, occupation, Fermi level, and so on. This is achieved through a self-consistent field (SCF) calculation on the auxiliary Kohn-Sham system. 
1. Create an input file ```scf.in```:
```bash
  &control
    calculation = 'scf',
    verbosity='high'
    pseudo_dir   = "../psps/",
    outdir       = ".",
    wf_collect=.TRUE.,
    prefix='LiF',
 /
 &system
    ibrav = 2, nat = 2, ntyp = 2,
    force_symmorphic=.TRUE. ,
    ecutwfc = 80.0,
    nbnd = 10 ,
    celldm(1)=7.703476,
 /
 &electrons
    diago_full_acc = .TRUE.,
    conv_thr =  1.0d-10
 /
 ATOMIC_SPECIES
    Li    6.941       Li.LDA.cpi.UPF
     F    18.998403   F.LDA.cpi.UPF
 ATOMIC_POSITIONS crystal
Li  0. 0. 0.
 F  0.5 -0.5 -0.5
K_POINTS automatic
 4  4  4  0  0  0
```
2. run with:
```
mpirun -np 4 pw.x < scf.in |tee log_scf.out
```
### non-self-consistent calculation
A **non-self-consistent calculation (nscf)** generate a set of Kohn-Sham eigenvalues and eigenvectors for both occupied and unoccupied bands. The output of this calculation are needed for generating the Yambo databases accurately.
1. Create and go in a new directory `$ mkdir nscf`
2. Copy the `.save` folder from the previous `scf` run `$ cp -r ../scf/LiF.save .` 
3. Create an `nscf.in` input file:
```bash
 &control
    calculation = 'nscf',
    verbosity='high'
    pseudo_dir   = "../psps"
    outdir       = "."
    wf_collect=.TRUE.,
    prefix='LiF',
 /
 &system
    ibrav = 2, nat = 2, ntyp = 2,
    force_symmorphic=.TRUE. ,
    ecutwfc = 80.0,
    nbnd = 50 ,
    celldm(1)=7.703476,
 /
 &electrons
    diago_full_acc = .TRUE.,
    conv_thr =  1.0d-10
 /
 ATOMIC_SPECIES
    Li    6.941       Li.LDA.cpi.UPF
     F    18.998403   F.LDA.cpi.UPF
 ATOMIC_POSITIONS crystal
Li  0. 0. 0.
 F  0.5 -0.5 -0.5
K_POINTS automatic
 4  4  4  0  0  0
```
4. run:
```
$ mpirun -np 4 pw.x < nscf.in |tee log_nscf.out
```
### Band structure
In order to gain a deeper understanding of the electronic states of the system, we proceed to compute the band structure along a path of high-symmetry points in the Brillouin Zone (BZ), specifically tracing the $W-L-\Gamma-X-W$ line.
We detailed our analysis, employing the `projwfc.x` utility to compute wavefunctions projected onto atomic orbitals. This allows us to examine the band structure with an emphasis on orbital contributions, providing a more detailed view of the electronic states and their interactions with the atomic structure of the material.
<div style="text-align: center;">
    <img src="https://hackmd.io/_uploads/r1EIqE-xC.png" alt="bz-path" width=40% pos=center>
</div> 
1. Create and move to a new `bands` folder
2. Copy the `scf` `.save` folder in the current directory `$ cp -r ../scf/LiF.save .`
3. create a new input file `nscf.in` (in QE a `bands` and `nscf` calculation are equivalent)

```
  &control
    calculation = 'nscf',
    verbosity='high'
    pseudo_dir   = "../psps"
    outdir       = "."
    wf_collect=.TRUE.,
    prefix='LiF',
 /
 &system
    ibrav = 2, nat = 2, ntyp = 2,
    force_symmorphic=.TRUE. ,
    ecutwfc = 80.0,
    nbnd = 10 ,
    celldm(1)=7.703476,
 /
 &electrons
    diago_full_acc = .TRUE.,
    conv_thr =  1.0d-10
 /
 ATOMIC_SPECIES
    Li    6.941       Li.LDA.cpi.UPF
     F    18.998403   F.LDA.cpi.UPF
 ATOMIC_POSITIONS crystal
Li  0. 0. 0.
 F  0.5 -0.5 -0.5
K_POINTS crystal_b
121
0.5    0.25    0.75    1
0.5    0.25833333333333336    0.7416666666666667    1
...
0.5    0.25    0.75    1
```
Note that we created a list of k-points along the high-symmetry path using **yambopy**. After importing the correct libraries you can create an instance of the `Path` class along the given path with `npoints` between each high-symmetry point. 
```python=
# Define path in reduced coordinates using Class Path
npoints = 30
path_kpoints = Path([[[  0.5,  0.250,  0.750],'W'],
                     [[0.5,0.5,  0.5],'L'],  
              [[  0.0,  0.0,  0.0],'$\Gamma$'],
              [[  0.5,  0.0,  0.5],'X'],
              [[  0.5,  0.250,  0.750],'W']],[npoints,npoints,npoints,npoints] )
```
2. run: 
``` $ mpirun -np 4 pw.x < nscf.in |tee log_nscf.out```
3. Compute the orbital 
 wavefunction. Create a new `bands.in` file
```bash
&projwfc
    prefix = 'LiF',
    outdir = './',
    ngauss=0, degauss=0.036748
    kresolveddos=.true. ! optional
    DeltaE=0.01
    filpdos = 'LiF-pdos.dat'
 /
```
an run 
```
$ mpirun -np 4 projwfc.x < bands.in |tee projwfc.log
```
**Note** the output file has to be named `projwfc.log` for yambopy

4. Inspect the output file `projwfc.log` to get the indices of the state belonging to a given atom's orbital
```
     Atomic states used for projection
     (read from pseudopotential files):

     state #   1: atom   1 (Li ), wfc  1 (l=0 m= 1)
     state #   2: atom   1 (Li ), wfc  2 (l=1 m= 1)
...
```
5. Plot the orbital projected band structure via `yambopy`:
```python=
# import libraries
from yambopy import *
import matplotlib.pyplot as plt
from yambopy.lattice import car_red, red_car
import matplotlib.pylab as pylab
from matplotlib.ticker import MultipleLocator
import matplotlib.lines as mlines

# create a figure and ax
fig,ax = plt.subplots(dpi=300)

dotsize = 10 #size of the dots, change at will
# Read the indices from `projwfc.log`. Remember python counting starts from 0

atom_Li_s = [0]
atom_Li_p = [1,2,3]
atom_Li_d = [4,5,6,7,8]
atom_Li_f = [9,10,11,12,13,14,15,15]
atom_F_s = [16]
atom_F_p = [17,18,19]
atom_F_d = [20,21,22,23,24,25,26]
atom_F_f = [27,28,29,30,31]
# get ticks and labels, need previous istance of Path
ticks, labels =list(zip(*path_kpoints.get_indexes()))

# create an instance of ProjwfcXML class
band = ProjwfcXML(prefix='LiF',path=f'./bands/',qe_version='6.7')
nelectrons = 8
Li_s = band.plot_eigen(ax,path_kpoints=path_kpoints,selected_orbitals=atom_Li_s,color='pink',size=dotsize)
Li_p = band.plot_eigen(ax,path_kpoints=path_kpoints,selected_orbitals=atom_Li_p,color='yellow',size=dotsize)
Li_d = band.plot_eigen(ax,path_kpoints=path_kpoints,selected_orbitals=atom_Li_d,color='orange',size=dotsize)
Li_f = band.plot_eigen(ax,path_kpoints=path_kpoints,selected_orbitals=atom_Li_f,color='red',size=dotsize)
F_s = band.plot_eigen(ax,path_kpoints=path_kpoints,selected_orbitals=atom_F_s,color='blue',size=dotsize)
F_p = band.plot_eigen(ax,path_kpoints=path_kpoints,selected_orbitals=atom_F_p,color='indigo',size=dotsize)
F_d = band.plot_eigen(ax,path_kpoints=path_kpoints,selected_orbitals=atom_F_d,color='cyan',size=dotsize)
F_f = band.plot_eigen(ax,path_kpoints=path_kpoints,selected_orbitals=atom_F_f,color='green',size=dotsize)

# setting labels and legend of the plot
Li_s.set_label(r'$Li_s$')
Li_p.set_label(r"$Li_p$")
Li_d.set_label(r'$Li_d$')
Li_f.set_label(r"$Li_f$")
F_s.set_label(r'$F_s$')
F_p.set_label(r'$F_p$')
F_d.set_label(r'$F_d$')
F_f.set_label(r'$F_f$')
ax.set_ylabel('E [eV]')
#ax.set_ylim([-5,20])
lgnd = plt.legend(loc=(1.04, 0), scatterpoints=1, markerscale=20, fontsize=15)
lgnd.legendHandles[0]._sizes = [dotsize]
lgnd.legendHandles[1]._sizes = [dotsize]
lgnd.legendHandles[2]._sizes = [dotsize]
lgnd.legendHandles[3]._sizes = [dotsize]
lgnd.legendHandles[4]._sizes = [dotsize]
lgnd.legendHandles[5]._sizes = [dotsize]
lgnd.legendHandles[6]._sizes = [dotsize]
lgnd.legendHandles[7]._sizes = [dotsize]
#save figure
plt.savefig(f'./bands/full_orbbands.pdf',bbox_inches='tight')
```
You can play with the aesthetic and plotting options. The final result should look like this

<img src="https://hackmd.io/_uploads/HJA3REbeR.png" alt="full_orbbands" style="width: 48%; margin-right: 2%;">
<img src="https://hackmd.io/_uploads/ry-8Jr-lA.png" alt="orbbands" style="width: 48%;">


Examining the two figures reveals that the top three valence bands predominantly consist of Fluorine p orbitals, whereas the initial conduction band exhibits a composition of mixed orbitals, including $Li_p$, $Li_s$, $F_s$, and $F_p$. This detail is crucial for the Wannierization process and presents a slight deviation from the analysis conducted in reference [^3], which focused solely on $Li_s$ and $F_p$ orbitals.
### Wannierization
**Requirement**: Navigate to the `.save` directory within the [nscf](#non-self-consistent calculation) section and initialize the Yambo database using the `p2y` and `yambo` commands.

Before the Wannier90 analysis, it's necessary to perform an [nscf](#non-self-consistent calculation) employing a uniform k-grid. We will explore calculations using two different grids: a k-grid of 8x8x8 and a q-grid of 4x4x4, focusing our discussion on the q-grid.

1. To create a uniform k-grid, you could use the kmesh.pl utility available in the `wannier90` package. Yet, to maintain consistency in the k-grid utilized by both Yambo and QE, we opt for yambopy to generate a comprehensive list of k-points across the full Brillouin Zone (BZ).
```python=
# 
savedb_q = YamboSaveDB.from_db_file(f'{YAMBO_TUT_PATH}/unshifted-grid/SAVE')
lat_q = YamboLatticeDB.from_db_file(f'{YAMBO_TUT_PATH}/unshifted-grid/SAVE/ns.db1')
full_kpoints, kpoints_indexes, symmetry_indexes=savedb_q.expand_kpts()
full_kpoints_red = car_red(full_kpoints, lat_q.rlat)

# Generate list of q-points with weights

for i in range (len(full_kpoints_red)):
    print(f'{full_kpoints_red[i][0]:.10f}   {full_kpoints_red[i][1]:.10f}   {full_kpoints_red[i][2]:.10f} {1/len(full_kpoints_red):.10f}')
```
2. Create an `nscf.in` input file
```
  &control
    calculation = 'nscf',
    verbosity='high'
    pseudo_dir   = "../psps"
    outdir       = "."
    wf_collect=.TRUE.,
    prefix='LiF',
 /
 &system
    ibrav = 2, nat = 2, ntyp = 2,
    force_symmorphic=.TRUE. ,
    ecutwfc = 80.0,
    nbnd = 10 ,
    celldm(1)=7.703476,
 /
 &electrons
    diago_full_acc = .TRUE.,
    conv_thr =  1.0d-10
 /
 ATOMIC_SPECIES
    Li    6.941       Li.LDA.cpi.UPF
     F    18.998403   F.LDA.cpi.UPF
 ATOMIC_POSITIONS crystal
Li  0. 0. 0.
 F  0.5 -0.5 -0.5
K_POINTS crystal
64
  0.00000000  0.00000000  0.00000000  1.562500e-02
  0.00000000  0.00000000  0.25000000  1.562500e-02
  0.00000000  0.00000000  0.50000000  1.562500e-02
  0.00000000  0.00000000  0.75000000  1.562500e-02
  0.00000000  0.25000000  0.00000000  1.562500e-02
  0.00000000  0.25000000  0.25000000  1.562500e-02
  0.00000000  0.25000000  0.50000000  1.562500e-02
  0.00000000  0.25000000  0.75000000  1.562500e-02
  0.00000000  0.50000000  0.00000000  1.562500e-02
  0.00000000  0.50000000  0.25000000  1.562500e-02
  0.00000000  0.50000000  0.50000000  1.562500e-02
  0.00000000  0.50000000  0.75000000  1.562500e-02
  0.00000000  0.75000000  0.00000000  1.562500e-02
  0.00000000  0.75000000  0.25000000  1.562500e-02
  0.00000000  0.75000000  0.50000000  1.562500e-02
  0.00000000  0.75000000  0.75000000  1.562500e-02
  0.25000000  0.00000000  0.00000000  1.562500e-02
  0.25000000  0.00000000  0.25000000  1.562500e-02
  0.25000000  0.00000000  0.50000000  1.562500e-02
  0.25000000  0.00000000  0.75000000  1.562500e-02
  0.25000000  0.25000000  0.00000000  1.562500e-02
  0.25000000  0.25000000  0.25000000  1.562500e-02
  0.25000000  0.25000000  0.50000000  1.562500e-02
  0.25000000  0.25000000  0.75000000  1.562500e-02
  0.25000000  0.50000000  0.00000000  1.562500e-02
  0.25000000  0.50000000  0.25000000  1.562500e-02
  0.25000000  0.50000000  0.50000000  1.562500e-02
  0.25000000  0.50000000  0.75000000  1.562500e-02
  0.25000000  0.75000000  0.00000000  1.562500e-02
  0.25000000  0.75000000  0.25000000  1.562500e-02
  0.25000000  0.75000000  0.50000000  1.562500e-02
  0.25000000  0.75000000  0.75000000  1.562500e-02
  0.50000000  0.00000000  0.00000000  1.562500e-02
  0.50000000  0.00000000  0.25000000  1.562500e-02
  0.50000000  0.00000000  0.50000000  1.562500e-02
  0.50000000  0.00000000  0.75000000  1.562500e-02
  0.50000000  0.25000000  0.00000000  1.562500e-02
  0.50000000  0.25000000  0.25000000  1.562500e-02
  0.50000000  0.25000000  0.50000000  1.562500e-02
  0.50000000  0.25000000  0.75000000  1.562500e-02
  0.50000000  0.50000000  0.00000000  1.562500e-02
  0.50000000  0.50000000  0.25000000  1.562500e-02
  0.50000000  0.50000000  0.50000000  1.562500e-02
  0.50000000  0.50000000  0.75000000  1.562500e-02
  0.50000000  0.75000000  0.00000000  1.562500e-02
  0.50000000  0.75000000  0.25000000  1.562500e-02
  0.50000000  0.75000000  0.50000000  1.562500e-02
  0.50000000  0.75000000  0.75000000  1.562500e-02
  0.75000000  0.00000000  0.00000000  1.562500e-02
  0.75000000  0.00000000  0.25000000  1.562500e-02
  0.75000000  0.00000000  0.50000000  1.562500e-02
  0.75000000  0.00000000  0.75000000  1.562500e-02
  0.75000000  0.25000000  0.00000000  1.562500e-02
  0.75000000  0.25000000  0.25000000  1.562500e-02
  0.75000000  0.25000000  0.50000000  1.562500e-02
  0.75000000  0.25000000  0.75000000  1.562500e-02
  0.75000000  0.50000000  0.00000000  1.562500e-02
  0.75000000  0.50000000  0.25000000  1.562500e-02
  0.75000000  0.50000000  0.50000000  1.562500e-02
  0.75000000  0.50000000  0.75000000  1.562500e-02
  0.75000000  0.75000000  0.00000000  1.562500e-02
  0.75000000  0.75000000  0.25000000  1.562500e-02
  0.75000000  0.75000000  0.50000000  1.562500e-02
  0.75000000  0.75000000  0.75000000  1.562500e-02
```
3. Create an input file for wannier90 `LiF.win`:
```bash=
num_bands         =   10
num_wann          =   8


DIS_WIN_MIN = -25
DIS_WIN_MAX = 24
dis_num_iter      = 1000
dis_mix_ratio     = 0.4

num_iter          = 3000
iprint    = 3
!num_dump_cycles = 10
!num_print_cycles = 1000


bands_num_points = 201
bands_plot_format = gnuplot

write_rmn = true
write_tb  = true
write_xyz = .true.
wannier_plot_supercell = 3
use_ws_distance = .true.
translate_home_cell = .true.

bands_plot = true
write_hr = true
Fermi_energy = 3.63060
Begin Atoms_Frac
F            0.5000000000       -0.5000000000      -0.5000000000
Li           0.0000000000       0.0000000000       0.0000000000
End Atoms_Frac

Begin Projections
Li :l=0;l=1
F : l=0;l=1
End Projections

Begin kpoint_path
W  0.500    0.250    0.750 L  0.500    0.500    0.500
L  0.500    0.500    0.500 Γ  0.000    0.000    0.00
Γ  0.000    0.000    0.00 X  0.500    0.000    0.500
X  0.500    0.000    0.500 W  0.500    0.250    0.750
End kpoint_path

Begin Unit_Cell_Cart
Bohr
     -3.851738  0.0000000000      3.851738
     0.000000000000000   3.851738 3.851738
     -3.851738  3.851738 0.000000000000000
End Unit_Cell_Cart

mp_grid      = 4 4 4

!exclude_bands 1,6,7,8,9,10
Begin kpoints
  0.00000000  0.00000000  0.00000000  1.562500e-02
  0.00000000  0.00000000  0.25000000  1.562500e-02
  0.00000000  0.00000000  0.50000000  1.562500e-02
  0.00000000  0.00000000  0.75000000  1.562500e-02
  0.00000000  0.25000000  0.00000000  1.562500e-02
  0.00000000  0.25000000  0.25000000  1.562500e-02
  0.00000000  0.25000000  0.50000000  1.562500e-02
  0.00000000  0.25000000  0.75000000  1.562500e-02
  0.00000000  0.50000000  0.00000000  1.562500e-02
  0.00000000  0.50000000  0.25000000  1.562500e-02
  0.00000000  0.50000000  0.50000000  1.562500e-02
  0.00000000  0.50000000  0.75000000  1.562500e-02
  0.00000000  0.75000000  0.00000000  1.562500e-02
  0.00000000  0.75000000  0.25000000  1.562500e-02
  0.00000000  0.75000000  0.50000000  1.562500e-02
  0.00000000  0.75000000  0.75000000  1.562500e-02
  0.25000000  0.00000000  0.00000000  1.562500e-02
  0.25000000  0.00000000  0.25000000  1.562500e-02
  0.25000000  0.00000000  0.50000000  1.562500e-02
  0.25000000  0.00000000  0.75000000  1.562500e-02
  0.25000000  0.25000000  0.00000000  1.562500e-02
  0.25000000  0.25000000  0.25000000  1.562500e-02
  0.25000000  0.25000000  0.50000000  1.562500e-02
  0.25000000  0.25000000  0.75000000  1.562500e-02
  0.25000000  0.50000000  0.00000000  1.562500e-02
  0.25000000  0.50000000  0.25000000  1.562500e-02
  0.25000000  0.50000000  0.50000000  1.562500e-02
  0.25000000  0.50000000  0.75000000  1.562500e-02
  0.25000000  0.75000000  0.00000000  1.562500e-02
  0.25000000  0.75000000  0.25000000  1.562500e-02
  0.25000000  0.75000000  0.50000000  1.562500e-02
  0.25000000  0.75000000  0.75000000  1.562500e-02
  0.50000000  0.00000000  0.00000000  1.562500e-02
  0.50000000  0.00000000  0.25000000  1.562500e-02
  0.50000000  0.00000000  0.50000000  1.562500e-02
  0.50000000  0.00000000  0.75000000  1.562500e-02
  0.50000000  0.25000000  0.00000000  1.562500e-02
  0.50000000  0.25000000  0.25000000  1.562500e-02
  0.50000000  0.25000000  0.50000000  1.562500e-02
  0.50000000  0.25000000  0.75000000  1.562500e-02
  0.50000000  0.50000000  0.00000000  1.562500e-02
  0.50000000  0.50000000  0.25000000  1.562500e-02
  0.50000000  0.50000000  0.50000000  1.562500e-02
  0.50000000  0.50000000  0.75000000  1.562500e-02
  0.50000000  0.75000000  0.00000000  1.562500e-02
  0.50000000  0.75000000  0.25000000  1.562500e-02
  0.50000000  0.75000000  0.50000000  1.562500e-02
  0.50000000  0.75000000  0.75000000  1.562500e-02
  0.75000000  0.00000000  0.00000000  1.562500e-02
  0.75000000  0.00000000  0.25000000  1.562500e-02
  0.75000000  0.00000000  0.50000000  1.562500e-02
  0.75000000  0.00000000  0.75000000  1.562500e-02
  0.75000000  0.25000000  0.00000000  1.562500e-02
  0.75000000  0.25000000  0.25000000  1.562500e-02
  0.75000000  0.25000000  0.50000000  1.562500e-02
  0.75000000  0.25000000  0.75000000  1.562500e-02
  0.75000000  0.50000000  0.00000000  1.562500e-02
  0.75000000  0.50000000  0.25000000  1.562500e-02
  0.75000000  0.50000000  0.50000000  1.562500e-02
  0.75000000  0.50000000  0.75000000  1.562500e-02
  0.75000000  0.75000000  0.00000000  1.562500e-02
  0.75000000  0.75000000  0.25000000  1.562500e-02
  0.75000000  0.75000000  0.50000000  1.562500e-02
  0.75000000  0.75000000  0.75000000  1.562500e-02
end kpoints
```
For the initial projections see the discussion in [Band structure](#Band-structure)

3. Copy the `scf` `.save` folder `$ cp -r ../scf/LiF.save .` (for the k-grid you might want to run another scf run with an 8x8x8 grid to ensure consistency)
4. Run an `nscf` calculation
```
$ mpirun -np 4 pw.x < nscf.in |tee log_nscf.out
```
5. Run the pre-processing steps of Wannier90
```
$ mpirun -np 1 wannier90.x -pp LiF
```
6. Create an input file `pw2wan.in` for `pw2wannier90.x` and run `pw2wannier90.x`:
```bash 
&inputpp
outdir = './'
prefix = 'LiF'
seedname = 'LiF'
write_amn = .true.
write_mmn = .true.
write_unk = .true. !optional
/
```
```
$ mpirun -np 1 pw2wannier90.x < pw2wan.in |tee log_pw2wan.out
```
7. Finally run `wannier90.x`
```
$ mpirun -np 1 wannier90.x LiF
```
You should now have a `seedname_hr.dat` file containing the real-space Hamiltonian in the MLWF basis.
We can compare the electronic band structure obtained with `Wannier90` and `QE` via the following Python scripting.
```python=
import tbmodels
from yambopy import *
import matplotlib.pyplot as plt
from yambopy.lattice import car_red, red_car
import matplotlib.pylab as pylab
from pylab import rcParams
from matplotlib.ticker import MultipleLocator
import matplotlib.lines as mlines

WORK_PATH='./'
unshifted_gridpath='yambo_tutorials/unshifted-grid/'
# Read xml file from a bands calculation
xml = ProjwfcXML(prefix='LiF',path=f'{unshiftedgrid_path}/bands/',qe_version='6.7')

# Compare DFT vs Wannier band structure
from qepy.lattice import Path, calculate_distances 
# Class PwXML. QE database reading
xml = ProjwfcXML(prefix='LiF',path=f'{unshiftedgrid_path}/bands/',qe_version='6.7')
wann_bands = np.loadtxt(f'{unshiftedgrid_path}/nscf-wannier-qgrid/LiF_band.dat',usecols=(0,1))
# Class PwXML. QE database reading
fig, ax = plt.subplots()
kpoints_dists = calculate_distances(xml.kpoints[:xml.nkpoints])
#xml.plot_eigen_ax(ax, path_kpoints, y_offset=-1., lw =1, ylim=(-4,17))
for ib in range(xml.nbands):
  ax.plot(kpoints_dists, xml.eigen[:,ib], c='red',lw=1.0)
ax.scatter(wann_bands[:,0]/np.max(wann_bands[:,0])*np.max(kpoints_dists), wann_bands[:,1], c='black',s=0.5)
# tb_eb = tb_ebands.add_kpath_labels(ax)
legend_entries = [
    mlines.Line2D([], [], color='black', ls ='-', label='Wannier'),
    mlines.Line2D([], [], color='red', label='DFT'),
]

# # Add custom legend outside the loop
ax.legend(handles=legend_entries, loc='upper center', bbox_to_anchor=(0.5, 1.20), ncol=3)
# ax.plot(qe_en, qe_bands)
ax.set_ylabel('E [eV]')
ax.set_ylim(-5,26)
#plt.show()
plt.savefig(f'{unshiftedgrid_path}/bands/DFTvsWann_q.pdf',bbox_inches='tight')

```
We can now compare the results obtained Wannierizing the k- and q-grid with the DFT computed band structure.

<img src="https://hackmd.io/_uploads/BJlrTvZxC.png" alt="dftvswann_k" style="width: 48%; margin-right: 2%;">
<img src="https://hackmd.io/_uploads/ryiSaw-x0.png" alt="dftvswann_q" style="width: 48%;">

The two plots are reasonable in agreement but not ideal. We recover the degeneracy at W in the valence bands manifold. 
It is possible to increase the accurarcy of Wannier interpolation increasing the size of the grid. 
Notably, if I were to include only $Li_s$ and $F_p$ orbitals as in ref[^3] is difficult to achieve a reasonable band structure and it might be necesssary to increase the `num_iter=1000000`. From the orbital projected band structure it's evident that these orbitals are not enough to describe the first conduction band and upmost three valence band. 

# Excited state - YAMBO
**Requirements**: Run the DFT [scf and nscf](#DFT) calculations. You can download all input and references file [here](https://media.yambo-code.eu/educational/tutorials/files/LiF.tar.gz)

Now, we transition to examining the excited state properties calculated via the ab initio Many-Body Perturbation Theory (Ai-MBPT) approach, facilitated by the Yambo code. Contrary to the approach detailed [previously](#DFT), here we employ shifted k-grids, as demonstrated below: 
```
...
K_POINTS automatic
 4  4  4  1  1  1
```
This phase involves a comparative analysis between results derived from *shifted* and *unshifted* grids. Assuming the [preceding steps](#DFT) were executed without issues, you should now possess a `LiF.save` directory:
```
$ ls LiF.save
charge-density.dat  data-file-schema.xml  F.LDA.cpi.UPF  Li.LDA.cpi.UPF  wfc10.dat  wfc1.dat  wfc2.dat	wfc3.dat  wfc4.dat  wfc5.dat  wfc6.dat	wfc7.dat 
```
## Conversion to Yambo format
The PWscf `LiF.save` output is converted to the Yambo format using the `p2y` executable, found in the yambo `bin` directory. Enter `LiF.save` and launch `p2y`:
```
$ cd LiF.save
$ p2y
...
 <--->  XC potential           : Slater exchange(X)+Perdew & Zunger(C)
 <--->  Atomic species        :  2
 <--->  Max atoms/species     :  1
 <---> == DB1 (Gvecs and more) ...
 <---> ... Database done
 <---> == DB2 (wavefunctions)  ...
 <---> [p2y] WF I/O |                                        | [000%] --(E) --(X)
 <---> [p2y] WF I/O |########################################| [100%] --(E) --(X)
 <---> == DB3 (PseudoPotential) ...
 <--->  == P2Y completed ==
```
This output repeats some information about the system and generates a `SAVE` directory that contains among others the `ns.db1` database with information about the crystal structure of the system.
## Initialization
The initialization procedure is common to all yambo runs.
1. Enter the `YAMBO` directory (where you have previously created the `SAVE` folder) and type
```bash
$ yambo -F Inputs/01_init -J 01_init
```
and you will see:
```
 <---> [03.01] X indexes
 <---> X [eval] |                                        | [000%] --(E) --(X)
 <---> X [eval] |########################################| [100%] --(E) --(X)
 <---> X[REDUX] |                                        | [000%] --(E) --(X)
 <---> X[REDUX] |########################################| [100%] --(E) --(X)
 <---> [03.01.01] Sigma indexes
 <---> Sigma [eval] |                                        | [000%] --(E) --(X)
 <---> Sigma [eval] |########################################| [100%] --(E) --(X)
 <---> Sigma[REDUX] |                                        | [000%] --(E) --(X)
 <---> Sigma[REDUX] |########################################| [100%] --(E) --(X)
 <---> [04] Timing Overview
 <---> [05] Memory Overview
 <---> [06] Game Over & Game summary
```
this run produced an `r_setup` file with a lot of useful information about your system, like the Fermi level, the crystal structure, the number of k-points, the gap etc...

## Random-Phase approximation

In this section, we calculate the absorption spectrum of LiF under various levels of approximation for the kernel. The most fundamental approach to estimating LiF's absorption spectrum is the independent particle approximation (`02_RPA_no_LF`). Following this, we introduce a kernel based on the Hartree approximation (`03_RPA_LF`), aligning with the Random-Phase-Approximation (RPA) framework. Moreover, to closely align with experimental absorption peak data, we employ a rigid scissor shift of the electronic states(`03_RPA_LF_QP`)
1. Create an input file called `02_RPA_no_LF`
```
optics                       # [R OPT] Optics
chi                          # [R LR] Linear Response.
% QpntsRXd
  1 | 1 |                   # [Xd] Transferred momenta
%
% BndsRnXd
   1 |  10 |                 # [Xd] Polarization function bands
%
NGsBlkXd= 1              RL  # [Xd] Response block size
% EnRngeXd
  7.50000 | 25.00000 |   eV  # [Xd] Energy range
%
% DmRngeXd
  0.10000 |  0.30000 |   eV  # [Xd] Damping range
%
ETStpsXd= 300                # [Xd] Total Energy steps
% LongDrXd
 1.000000 | 0.000000 | 0.000000 |      # [Xd] [cc] Electric Field
%
```
2. Run yambo:
```
$ yambo -F 02_RPA_no_LF -J 02_RPA_no_LF
```
The optional *-J* flag is used to label the output/report/log files.

3. Create the input files `03_RPA_LF` and `03_RPA_LF_QP`:
```
optics                       # [R OPT] Optics
chi                          # [R LR] Linear Response.
% QpntsRXd
  1 | 1 |                   # [Xd] Transferred momenta
%
% BndsRnXd
   1 |  10 |                 # [Xd] Polarization function bands
%
NGsBlkXd=51              RL  # [Xd] Response block size
% EnRngeXd
  7.50000 | 25.00000 |   eV  # [Xd] Energy range
%
% DmRngeXd
  0.10000 |  0.30000 |   eV  # [Xd] Damping range
%
ETStpsXd= 300                # [Xd] Total Energy steps
% LongDrXd
 1.000000 | 0.000000 | 0.000000 |      # [Xd] [cc] Electric Field
%
```
<span style="color: red;">
% XfnQP_E 
    
 5.190000 | 1.000000 | 1.000000 |      # [EXTQP Xd] E parameters (c/v), scissor
    
%
    
</span> 

and run yambo
```
$ yambo -F 03_RPA_LF -J -03_RPA_LF
$ yambo -F 03_RPA_LF -J -03_RPA_LF_QP
```

4. Plot the output data via
```python 
# read files
data_02_RPA_no_LF = np.loadtxt(f'{YAMBO_TUT_PATH}/shifted-grid/o-02_RPA_no_LF.eps_q1_inv_rpa_dyson', usecols=(0,1))
data_03_RPA_LF = np.loadtxt(f'{YAMBO_TUT_PATH}/shifted-grid/o-03_RPA_LF.eps_q1_inv_rpa_dyson', usecols=(0,1))
data_03_RPA_LF_QP = np.loadtxt(f'{YAMBO_TUT_PATH}/shifted-grid/o-03_RPA_LF_QP.eps_q1_inv_rpa_dyson', usecols=(0,1))
data_exp = np.loadtxt(f'{YAMBO_TUT_PATH}/e2_experimental.dat', usecols=(0,1))
#plot
fig,ax =plt.subplots()
ax.plot(data_02_RPA_no_LF[:,0], data_02_RPA_no_LF[:,1], label='02_RPA_no_LF')
ax.plot(data_03_RPA_LF[:,0], data_03_RPA_LF[:,1], label='03_RPA_LF')
ax.plot(data_03_RPA_LF_QP[:,0], data_03_RPA_LF_QP[:,1], label='03_RPA_LF_QP')
ax.plot(data_exp[:,0], data_exp[:,1], label='Experiment')
ax.set_ylabel('Absorption [a.u.]')
ax.set_xlabel('Energy [eV]')
ax.legend()
plt.savefig(f'{YAMBO_TUT_PATH}/shifted-grid/absorption_tut_v1.png', bbox_inches='tight')
```

![absorption_tut_v1](https://hackmd.io/_uploads/SJnNyFZl0.png)
 
As we can see, the energy shifting helps with matching the experimental energy onset but there is an experimental peak below the fundamental peak which is not captured by our current level of approximation. This state is an **EXCITON**.

## Adiabatic LDA (ALDA) approximation
Our first attempt to go beyond RPA is using TDDFT in the Adiabatic LDA approximation. Depending on the exchange-correlation (xc) functional used in the ground state calculation yambo will produce a corresponding xc-kernel. In this case, we are using a Perdew Wang parametrization.
Create a `04_alda_g_space` file
```bash
optics                       # [R OPT] Optics
chi                          # [R LR] Linear Response.
alda_fxc                     # [R TDDFT] The ALDA TDDFT kernel
% QpntsRXd
 1 | 1 |                     # [Xd] Transferred momenta
%
% BndsRnXd
  1 | 10 |                   # [Xd] Polarization function bands
%
NGsBlkXd=  51            RL  # [Xd] Response block size
% EnRngeXd
  7.50000 | 25.00000 |   eV  # [Xd] Energy range
%
% DmRngeXd
 0.100000 | 0.300000 |   eV  # [Xd] Damping range
%
ETStpsXd= 300                # [Xd] Total Energy steps
% LongDrXd
 1.000000 | 0.000000 | 0.000000 |      # [Xd] [cc] Electric Field
%
FxcGRLc= 51              RL  # [TDDFT] XC-kernel RL size
```
and run `yambo` again
```
$ yambo -F 04_alda_g_space -J 04_alda_g_space
```
If we plot again, the result does not improve considerably.

![absorption_tut_v2](https://hackmd.io/_uploads/S1ZlGYWgA.png)

## The Statically screened Electron-electron interaction
As simple LDA and ALDA approximations have not described correctly the experimental spectrum we have tom ove towads more elaborate techniques like the Bethe-Salpeter equation (BSE).
A key ingredient in the BSE kernel is the *electron-electron interaction* commonly evaluated in the static approximation. The input file `05_W` describes how to calculate it
```
em1s                         # [R Xs] Static Inverse Dielectric Matrix
% QpntsRXs
  1 | 19 |                   # [Xs] Transferred momenta
%
% BndsRnXs
  1 | 50 |                   # [Xs] Polarization function bands
%
NGsBlkXs=51              RL  # [Xs] Response block size
% LongDrXs
 1.000000 | 0.000000 | 0.000000 |      # [Xs] [cc] Electric Field
%
```
The variables used in htis input file have the same physical meaning of those used in the optical absorption calculation. The only difference is that, in general, the response function dimension obtained in the examples (03-04), gives an upper bound to the number of RL vectors needed here. This is because the size of response matrix in an RPA calculation defines also the size of the Hartree potential, whose short-range components are not screened. In the present case, instead, the electron-electron interaction is screened and, for this reason, the RL vectors are considerably smaller than in the RPA case.
Run this calculation via
```
$ yambo -F 05_W
```
to create the `ndb.em1s` database which will be read in the BSE calculation.
## The BSE equation, Excitons
In this subsection, we calculate the excitonic absorption spectrum, solving the BSE equation.
1. Create an input file `06_BSE`
```
optics                       # [R OPT] Optics
bse                          # [R BSK] Bethe Salpeter Equation.
bss                          # [R BSS] Bethe Salpeter Equation solver
% KfnQP_E
 5.80000 | 1.000000 | 1.000000 |      # [EXTQP BSK BSS] E parameters (c/v), scissor
%
BSKmod= "SEX"              # [BSK] Resonant Kernel mode. (`x`;`c`;`d`)
% BSEBands
  2 |  7 |                   # [BSK] Bands range
%
BSENGBlk=  51            RL  # [BSK] Screened interaction block size
BSENGexx=  773           RL  # [BSK] Exchange components
BSSmod= "d"                  # [BSS] Solvers `h/d/i/t`
% BEnRange
  7.50000 | 25.00000 |   eV  # [BSS] Energy range
%
% BDmRange
  0.15000 |  0.30000 |   eV  # [BSS] Damping range
%
BEnSteps= 300                # [BSS] Energy steps
% BLongDir
 1.000000 | 0.000000 | 0.000000 |      # [BSS] [cc] Electric Field
%
% BSEQptR
 1 | 19 |                             # [BSK] Transferred momenta range
%
```
Note that we are applying a 5.8 eV scissor and including both exchange and correlations terms in the kernel via the screened exchange (SEX) approximation for the kernel (`BSKmod="SEX"`)
Remember that the BSE kernel is written in Bloch space and its size is given by
```
BSE kernel size = Valence Bands x Conduction bands x K-points in the whole BZ
```

2. Run yambo:
```
$ yambo -F 06_BSE -J 06_BSE
```
3. Plot the output files and compare the results.
```python
...
data_06_BSE = np.loadtxt(f'{YAMBO_TUT_PATH}/shifted-grid/o-06_BSE.eps_q1_diago_bse', usecols=(0,1))
...
ax.plot(data_06_BSE[:,0], data_06_BSE[:,1], label='06_BSE')
```

![absorption_tut_v3](https://hackmd.io/_uploads/r1f3rFbgC.png)

We see that we have much better agreement with experiment and the BSE equation is able to describe the bound electron-hole state responsible for the peak observed experimentally below the QP gap.

## shifted vs unshifted grid
We can now repeat the same steps starting from an unshifted grid (`4 4 4 0 0 0`).
Note that in the unshifted case you have 64 k-points in the BZ, while in the shifted case you have 4 shifts for a total of 64*4=256 points in the BZ. 
![absorptionshitvsunshift](https://hackmd.io/_uploads/BydV9YWeR.png)

The two results differ only slightly.

**Notes for developers:**
1. Question: In the `r_setup` files I can see that the K-grid has 10 IBZ points and 256 full BZ while the Q-grid has 19 IBZ and 256 points. Why this difference?

# Exciton band structure
In this section, we will try to reproduce the exciton band structure from *Neaton et al.* [^3]


<a id="excbandsref"></a>
<figure>
    <img src="https://hackmd.io/_uploads/r1zOnKWxA.png" alt="Alt text for the image" style="width:60%">
    <figcaption>Singlet (left panel) and triplet (right panel) exciton dispersion for LiF. Linear interpolation of exciton eigenergeies explicitly explicitly obtained by diagonalizing the BSE at 72 Q points along the high-symmetry path (blue curve). Wannier-Fourier interpolated exciton dispersion starting from a centered 5 × 5 × 5 Q mesh (dotted red curve). The inset of the left panel shows how Wannier-Fourier interpolation performs without isolating the long-range contribution.
    </figcaption>
</figure>

We will work with the 2 *unshifted-grids* a (`4 4 4 0 0 0`) and an (`8 8 8 0 0 0`).
After solving the BSE equation for all 19 q-points in the BZ you will have a list of files `ndb.BS_diag_Q*` containing the eigenvalues and eigenvectors for each q-transferred momenta.
In order to plot the exciton band structure:
1. Create a `ypp` input file `ypp-excbands.in` for interpolation:
```
excitons                         # [R] Excitonic properties
interpolate                      # [R] Interpolate
INTERP_mode= "BOLTZ"                # Interpolation mode (NN=nearest point, BOLTZ=boltztrap aproach)
INTERP_Shell_Fac= 20.00000       # The bigger it is a higher number of shells is used
INTERP_NofNN= 1                  # Number of Nearest sites in the NN method
BANDS_steps= 30                  # Number of divisions
cooIn= "rlu"                     # Points coordinates (in) cc/rlu/iku/alat
cooOut= "rlu"                    # Points coordinates (out) cc/rlu/iku/alat
States= "1 - 4"                  # Index of the BS state(s)
% INTERP_Grid
-1 |-1 |-1 |                             # Interpolation BZ Grid
%
#PrtDOS                        # Print Exciton Density of States
% DOSERange
 1.000000 |-1.000000 |         eV    # Energy range
%
DOSESteps=  500                  # Energy steps
DOS_broad= 0.100000        eV    # Broadening of the DOS
%BANDS_kpts                      # K points of the bands circuit
0.5|0.5|0.5|
0.0|0.0|0.0|
0.5|0.0|0.5|
0.5|0.250|0.750|
0.5|0.50|0.50|
%
```
Note that, in my experience the `NN` interpolation mode does not work.
2. Run `ypp`
```
mpirun -np 4 ypp -F ypp-excbands.in -J 'bse-bands'
```
3. Plot the exciton band structure from the `o-bse-bands.excitons_interpolated` output file
```python=
...
#data_excbands_shifted = np.loadtxt(f'{YAMBO_TUT_PATH}/shifted-grid/bse-bands/o-bse-bands.excitons_interpolated', usecols=(0,1,2,3,4,5,6))
data_excbands_unshifted = np.loadtxt(f'{YAMBO_TUT_PATH}/unshifted-grid/bse-bands/o-bse-bands.excitons_interpolated', usecols=(0,1,2,3,4,5,6))
fig,ax =plt.subplots()
#ax.plot(data_excbands_shifted[:,0], data_excbands_shifted[:,1], label='06_BSE_shifted', c='blue')
ax.plot(data_excbands_unshifted[:,0], data_excbands_unshifted[:,1], label='06_BSE_unshifted',c='orange')
#ax.plot(data_excbands_shifted[:,0], data_excbands_shifted[:,2], c ='blue')
ax.plot(data_excbands_unshifted[:,0], data_excbands_unshifted[:,2], c='orange')
#ax.plot(data_excbands_shifted[:,0], data_excbands_shifted[:,3], c ='blue')
ax.plot(data_excbands_unshifted[:,0], data_excbands_unshifted[:,3], c='orange')
#ax.plot(data_excbands_shifted[:,0], data_excbands_shifted[:,4], c ='blue')
ax.plot(data_excbands_unshifted[:,0], data_excbands_unshifted[:,4], c='orange')
#ax.plot(data_excbands_shifted[:,0], data_excbands_shifted[:,5], c ='blue')
ax.plot(data_excbands_unshifted[:,0], data_excbands_unshifted[:,5], c='orange')
#ax.plot(data_excbands_shifted[:,0], data_excbands_shifted[:,6], c ='blue')
ax.plot(data_excbands_unshifted[:,0], data_excbands_unshifted[:,6], c='orange')
ticks, labels =list(zip(*path_kpoints.get_indexes()))
kpoints_dists= np.array(kpoints_dists)
ax.set_xticks(kpoints_dists[np.array(ticks)]/np.max(kpoints_dists[np.array(ticks)])*np.max(data_excbands_unshifted[:,0]), labels)
ax.set_ylabel('Energy [eV]')
ax.set_xlabel('$q_x$')
ax.legend(loc='upper left', fontsize=8)
plt.savefig(f'{YAMBO_TUT_PATH}/excbands_unshifted.pdf', bbox_inches='tight')
```
<figure>
    <img src="https://hackmd.io/_uploads/rJ_gY_feC.png" alt="Alt text for the image" style="width:100%">
    <figcaption>Exciton band structure computed starting with a k-grid = 4 4 4 0 0 0.
    </figcaption>
</figure>

If we compare this plot with the one above from Haber et al[^3] we see that the result do not match. Notably, the degeneracy at W is not recovered.
We can try to increase the `k-grid = 8 8 8 0 0 0` and plot the result as above
<figure>
    <img src="https://hackmd.io/_uploads/rybqcdMxC.png" alt="Alt text for the image" style="width:100%">
    <figcaption>Exciton band structure computed starting with a k-grid = 8 8 8 0 0 0.
    </figcaption>
</figure>

The result change noticeably and the degeneracy at W is almost recovered.
These discrepancies might be due to underconverged parameters in the BSE calculation.
We can compare the results obtained with shifted and unshifted grid

<figure>
    <img src="https://hackmd.io/_uploads/SyzkTOGeC.png" alt="Alt text for the image" style="width:100%">
    <figcaption>Exciton band structure computed starting with a k-grid = 8 8 8 0 0 0 (orange) and k-grid = 4 4 4 1 1 1 (blue).
    </figcaption>
</figure>

**Note for developers:** 
1. **Error:** the distances computed with `yambopy` and `ypp` differs. As you note in line 19 of the script above I had to rescale the `$q_x$` distances.
2. **Error:** Even though I asked in the `ypp` input for `BANDS_steps=30` I end up with 93 q-points instead of 121.

## Exciton band structure: singlets vs triplets
Yambo, by defaults, cmomputes singlets. In order to plot triplets we can set to "0 RL" the cutoff of the ex exchange in the BSE kernel.
1. Change the <span style="color: red;"> BSENGexx=  0 RL</span> input flag in the `06_BSE` input file and run again the BSE calculation
2. Compute the exciton band structure with `ypp` as above
3. Plot the exciton band structure, comparing singlets vs triplets results.

<figure>
    <img src="https://hackmd.io/_uploads/r1TBaFGeR.png" alt="Alt text for the image" style="width:100%">
    <figcaption>Exciton band structure for singlets (blue) and triplets(orange) components. The plots were generated computed starting with a k-grid = 4 4 4 1 1 1 (orange) and k-grid = 4 4 4 1 1 1 (blue).
    </figcaption>
</figure>

Note that now we recover the degenercy at $W$ and the result is not so different from the [reference](#excbandsref). There is still some discrepancy on the first two bands between $L$ and $\Gamma$ that should be degenerate.

# Useful Link
1. Download input files QE + Yambo [here](https://media.yambo-code.eu/educational/tutorials/files/LiF.tar.gz):
2. Access additional input files for Wannierization and Python notebook from branch `devel tb-wannier` of [this repo](https://github.com/rreho/yambopy/tree/master)

# TB-wannier

<figure>
    <img src="https://hackmd.io/_uploads/rkbDh9zlA.png" alt="Alt text for the image" style="width:100%">
    <figcaption> Utilities implemented in yambopy for exciton Wannierization.
    </figcaption>
</figure>

# Bibliography
[^1]: Kohn, W. (2019). Density functional theory. Introductory quantum mechanics with MATLAB: for atoms, molecules, clusters, and nanocrystals.
[^2]: Marzari, N., Mostofi, A. A., Yates, J. R., Souza, I., & Vanderbilt, D. (2012). Maximally localized Wannier functions: Theory and applications. Reviews of Modern Physics, 84(4), 1419.
[^3]: Haber, J. B., Qiu, D. Y., da Jornada, F. H., & Neaton, J. B. (2023). Maximally localized exciton Wannier functions for solids. Physical Review B, 108(12), 125118
[^4]: Sangalli, D., Ferretti, A., Miranda, H., Attaccalite, C., Marri, I., Cannuccia, E., ... & Marini, A. (2019). Many-body perturbation theory calculations using the yambo code. Journal of Physics: Condensed Matter, 31(32), 325902.
[^5]: Marini, A., Hogan, C., Grüning, M., & Varsano, D. (2009). Yambo: an ab initio tool for excited state calculations. Computer Physics Communications, 180(8), 1392-1403.
[^6]: Giannozzi, P., Baseggio, O., Bonfà, P., Brunato, D., Car, R., Carnimeo, I., ... & Baroni, S. (2020). Quantum ESPRESSO toward the exascale. The Journal of chemical physics, 152(15).