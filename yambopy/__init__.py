#
# License-Identifier: GPL
#
# Copyright (C) 2024 The Yambo Team
#
# Authors: HPC, FP
#
# This file is part of the yambopy project
#
"""
Create, read and write yambo input files
Read, modify and write yambo databases
Analyse results from yambo calculations

Modules:
    io
        - YamboIn: read, write and manipulate yambo input files
        - YamboOut: read yambo output files and save in .json
    dbs
        - YamboSaveDB: read information in the ns.db1
        - YamboLatticeDB: read lattice parameters, symmetries and k-points from ns.db1
        - YamboDipolesDB: dipole matrix elements from ndb.dip*
        - YamboStaticScreeningDB: static dielectric screening from ndb.em1s*
        - YamboElectronsDB: read the electronic states from ns.db1
        - YamboQPDB: read the quasiparticle energies db ndb.QP
        - YamboGreenDB: read the green's functions calculated using yambo
        - YamboExcitonsDB: read excitonic properties from a BSE calculation
        - YamboKernelDB: read excitonic kernel in transition space
        - YamboNLDB: read nonlinear response calculated with yambo_nl

    bse
        - YamboExcitonWaveFunctionXSF: read the excitonic
        - YamboExcitonWeight: read the excitonic weights from the ypp output file
        - YamboBSEAbsorptionSpectra: generate a .json file with the bse absorption calculation (including information about the excitons)

    em1s
        - YamboEm1sRotate: rotate em1s from IBZ to BZ
        - YamboEm1sExpand: expand em1s from unit cell to supercell [IN DEVELOPMENT]
    analyse:
        - YamboAnalyser: read .json files generated with yamboout and plot them together
"""
import numpy as np

class yambopyenv():
    YAMBO = "yambo"
    P2Y = "p2y"
    E2Y = "e2y"
    YPP = "ypp"
    SCHEDULER = "bash"
    YAMBO_RT = "yambo_rt"
    YPP_RT = "ypp_rt"
    YAMBO_NL = "yambo_nl"
    YPP_NL = "ypp_nl"

#tools and units
from yambopy.tools.jsonencoder import *
from yambopy.units import *

#lattice-related operations
from yambopy.lattice import *

#kpoint mesh operations
from yambopy.kpoints import *

#skw interpolator (adapted from abipy version)
from yambopy.tools.skw import *

#yambo databases
#from yambopy.dbs.savedb import *
from yambopy.dbs.dipolesdb import *
from yambopy.dbs.qpdb import *
from yambopy.dbs.hfdb import *
from yambopy.dbs.em1sdb import *
from yambopy.dbs.greendb import *
from yambopy.dbs.latticedb import *
from yambopy.dbs.electronsdb import *
from yambopy.dbs.rtdb import *
from yambopy.dbs.nldb import *
from yambopy.dbs.excitondb import *
from yambopy.dbs.wfdb import *
from yambopy.dbs.elphondb import *
from yambopy.dbs.bsekerneldb import *
from yambopy.dbs.excphondb import *
from yambopy.dbs.kqgridsdb import *

#input/output files
from yambopy.io.inputfile import *
from yambopy.io.outputfile import *
from yambopy.io.jsonfile import *
from yambopy.io.iofile import *
from yambopy.io.xsffile import *

#bse/excitons files
from yambopy.bse.excitonwf import *
from yambopy.bse.excitonweight import *
from yambopy.bse.bse_absorption import *
from yambopy.bse.bse_dispersion import *
from yambopy.bse.excitonradiativelifetimes import *

#em1s/static screening operations files
from yambopy.em1s.em1s_rotate import *

#ndb.QP operations
from yambopy.quasiparticles.QP_rotate import *

#LetzElPhC interface
from yambopy.letzelphc_interface.lelphcdb import *
from yambopy.letzelphc_interface.lelph2y import *

#analyse stuff
from yambopy.analyse import *

#workflow files
from yambopy.common.save_generation import *
from yambopy.common.workflow import *
from yambopy.common.calculation_manager import *
from yambopy.common.transform_matrix_element import *

#realtime files
from yambopy.rt.rt_movie import *
from yambopy.rt.rt_timestep_optimize import *

#non-linear files
from yambopy.nl.linear_optics import *
from yambopy.nl.fft_interp import *
from yambopy.nl.external_efield import *
from yambopy.nl.damp_it import *
from yambopy.nl.harmonic_analysis import *
from yambopy.nl.hhg_tools import *

#doublegrid files
from yambopy.double_grid.dg_convergence import *

#gkkp files
from yambopy.gkkp.compute_gkkp import *
from yambopy.gkkp.refine_gkkp import *

#data
from yambopy.data import *
