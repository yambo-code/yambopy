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

# Tutorial 3. BN (semiconductor). GW Band structure
==============================================

1. Find stretching coefficients

2. Plot DFT and GW band structures non-interpolated

3. Plot DFT and GW band structures interpolated

4. Compare non-interpolated and interpolated band structures

