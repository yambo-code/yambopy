# Copyright (C) 2015 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
#
"""
Create, read and write yambo input files
Read, modify and write yambo databases
Analyse results from yambo calculations

Modules:
    io
        - inputfile: read, write and manipulate yambo input files
        - outputfile: read yambo output files and save in .json
    dbs
        - savedb: read information in the nx.db1
        - latticedb: read lattice parameters, symmetries and k-points from ns.db1
        - dipolesdb: dipole matrix elements from ndb.dip*
        - em1sd: static dielectric screening from ndb.em1s*
        - electronsdb: read the electronic states from ns.db1
        - qpdb: read the quasiparticle energies db ndb.QP
        - electronsdb: read the electronic states from ns.db1

    bse
        - excitonwf: read the excitonic
        - excitonweight: read the excitonic weights from the ypp output file
        - bse_absorption: generate a .json file with the bse absorption calculation (including information about the excitons)

    analyse:
        - analyse: read .json files generated with yamboout and plot them together
        - recipes: user contibuted scripts
"""
import numpy as np
from yambopy.jsonencoder import *
from yambopy.netcdf import *
from yambopy.plot import *
from yambopy.units import *

#lattice stuff
from yambopy.lattice import *

#yambo databases
from yambopy.dbs.savedb import *
from yambopy.dbs.dipolesdb import *
from yambopy.dbs.qpdb import *
from yambopy.dbs.em1sdb import *
from yambopy.dbs.greendb import *
from yambopy.dbs.latticedb import *
from yambopy.dbs.electronsdb import *

#input/output files
from yambopy.io.inputfile import *
from yambopy.io.outputfile import *

#bse/excitons files
from yambopy.bse.excitonwf import *
from yambopy.bse.excitonweight import *
from yambopy.bse.bse_absorption import *

#analyse stuff
from yambopy.analyse import *
from yambopy.recipes import *
