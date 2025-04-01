Database Tutorial - Qepy
========

Here you find a basic postprocessing tutorial for quantum espresso databases.

Inside the directory there are several scripts showing some of the functionalities of these classes. The examples provided are not exhaustive and you can explore the code in qepy/\* (or write to the yambo forum) to investigate about all the features.

In particular we treat the following classes:
Scripts and classes DESCRIPTION.

# Download databases
You can download the relevant databases for this tutorial [here](www.yambo-code.org/educational/tutorials/files/databases\_qepy.tar) from the yambo website, or by simply typing:
> wget www.yambo-code.org/educational/tutorials/files/databases\_qepy.tar

# Yambo wiki
You can find an in-depth discussion of this tutorial on the yambo wiki [here](http://www.yambo-code.org/wiki/index.php?title=Yambopy_tutorial:_band_structures).

# Tutorial 1. BN (semiconductor). Band structure
==============================================

Folder 'bn-semiconductor'

1. Plot band structure

python plot-qe-bands.py

2. Plot the atomic orbital projected band structure

python plot-qe-orbitals.py

# Tutorial 2. Iron (metal). Band structure
==============================================

Folder 'iron-metal'

0. Calculate scf density and bands

python flow-iron.py

1. Plot band structure

python plot-qe-bands.py

2. Plot the atomic orbital projected band structure. The dot size is correlated
with the weight of the atomic orbitals.

python plot-qe-orbitals-size.py

3. Plot the atomic orbital projected band structure. We can set a colormap
   related with the partial weights of two sets of atomic orbitals.

python plot-qe-orbitals-colormap.py

# Tutorial 3. Unfolding supercell band structure
==============================================

Folder 'bn-semiconductor'

1. Calculate scf density and bands for pristine (pc) case.

2. Calculate scf density for supercell (sc) case - In this case, 'ecutwfc' and 'nbnd' have to be the same between pc and sc cases. However, lattice parameters for the sc have to be exactly a factor X with respect to the pc lattice parameter. If 'celldm(1) = 5.952637' for pc case, 'celldm(1) = 11.905274' for 2x2 sc case. In the case of a monolayer, if "celldm(1) = 5.919710945" for pc case, then "celldm(3) = 2.959855472" for 2x2 sc case. Besides, the atomic positions of repeated atoms between pc and sc have to be exactly the same.

2. Calculate bands for sc case - The BZ path has to be exactly a factor X with respect to the pc BZ path. If G-M is (0.0, 0.0, 0.0) - (0.0, 0.5, 0.0) in pc case, (0.0, 0.0, 0.0) - (0.0, 1.0, 0.0) for 2x2 sc case.

(Previous database can be downloaded using the link above)

2. Plot the pc band structure and unfolded sc band structure.

python plot-unfolding.py

# Tutorial 4. BN (semiconductor). GW Band structure
==============================================

1. Find stretching coefficients

2. Plot DFT and GW band structures non-interpolated

3. Plot DFT and GW band structures interpolated

4. Compare non-interpolated and interpolated band structures

