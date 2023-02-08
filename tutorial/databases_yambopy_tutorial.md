Database Tutorial - Yambopy
========

Here you find a basic postprocessing tutorial for yambo databases.

Inside the directory there are several scripts showing some of the functionalities of these classes. The examples provided are not exhaustive and you can explore the code in yambopy/dbs (or write to the yambo forum) to investigate about all the features.

In particular we treat the following classes:
1. bz\_plot.py: YamboLatticeDB (managing lattice information inside ns.db1)
2. elph\_plot.py: YamboElectronPhononDB (managing el-phon matrix elements from ndb.elph\_gkkp\*)
3. dipoles\_plot.py YamboDipolesDB (managing dipole el-light matrix elements from ndb.dipoles)
4. exc\_read.py, exc\_kspace\_plot.py, exc\_abs\_plot.py: YamboExcitonDB (managing exciton data from ndb.BS\_diagoQ\*)

# Download databases
The scripts provided are self-explanatory. In order to run them on a simple system, you can download the relative databases for monolayer hBN (12x12x1 k-grid and q-grid for electron-phonon and exciton data) [here](www.yambo-code.org/educational/tutorials/files/databases\_yambopy.tar) from the yambo website, or by simply typing:
> wget www.yambo-code.org/educational/tutorials/files/databases\_yambopy.tar

# Yambo wiki
You can find an in-depth discussion of this tutorial on the yambo wiki [here](http://www.yambo-code.org/wiki/index.php?title=Yambopy_tutorial:_Yambo_databases).
