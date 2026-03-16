Electron-Phonon Tutorial - Yambopy
========

Here you find a tutorial for dealing with electron-phonon coupling in yambopy.

Inside the directory there is a script showing some of the functionalities of the electron-phonon classes. The examples provided are not exhaustive and you can explore the code in yambopy/dbs and yambopy/letzelph\_interface (or write to the yambo forum) to investigate about all the features.

In particular we treat the following classes:
1. elph\_plot.py: YamboElectronPhononDB (managing el-phon matrix elements from ndb.elph\_gkkp\* -- the Yambo default format)
2. elph\_plot.py: LetzElphElectronPhononDB (managing el-phon matrix elements from ndb.elph -- the LetzElPhC code)

In addition, there is a tutorial on how to run the LetzElPhC code to produce the matrix elements: this requires compilation of both Quantum Espresso and LetzElPhC on the machine.

# Download databases
The scripts provided are self-explanatory. In order to run them on a simple system, you can download the relative databases for monolayer MoS2 (6x6x1 k and q-grid, including spin-orbit interaction) [here](www.yambo-code.org/educational/tutorials/files/yambopy\_electron\_phonon.tar.gz) from the yambo website, or by simply typing:
> wget www.yambo-code.org/educational/tutorials/files/databases\_yambopy.tar

# Yambo wiki
You can find an in-depth discussion of this tutorial on the yambo wiki [here](https://wiki.yambo-code.eu/wiki/index.php?title=Yambopy_tutorial:_electron-phonon_coupling).
