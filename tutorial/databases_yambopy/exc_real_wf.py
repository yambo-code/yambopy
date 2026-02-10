import numpy as np
from yambopy.dbs.excitondb import YamboExcitonDB
from yambopy.dbs import excitondb
from yambopy.dbs.latticedb import YamboLatticeDB
from yambopy.dbs.wfdb import YamboWFDB
import os

iqpt = 1 # qpt index of exciton
calc_path = '.'
BSE_dir = 'GW_BSE'

# load lattice db
lattice = YamboLatticeDB.from_db_file(os.path.join(calc_path, 'SAVE','ns.db1'))

# load exciton db
# note in case too many excitons, load only first few with `neigs' flag
# DO NOT forget to include all degenerate states when giving neigs flag !
#
filename = 'ndb.BS_diago_Q%d' % (iqpt)
excdb = YamboExcitonDB.from_db_file(lattice, filename=filename,
                                    folder=os.path.join(calc_path, BSE_dir),
                                    neigs = -1)

#Load the wavefunction database
wfdb = YamboWFDB(path='.', latdb=lattice,
                      bands_range=[np.min(excdb.table[:, 1]) - 1,
                      np.max(excdb.table[:, 2])])

## plot the exciton wavefunction with hole fixed at [0,0,0]
# in a [1,1,1] supercell with 80 Ry wf cutoff. (give -1 to use full cutoff)
# I want to set the degeneracy threshold to 0.01 eV
# For example I want to plot the 3rd exciton, so iexe = 2 (python indexing )
#
excdb.real_wf_to_cube(iexe=2, wfdb=wfdb, fixed_postion=[0.0      ,  0.0 ,       0.0],
                    supercell=[1,1,1], degen_tol=0.01, wfcCutoffRy=-1, fix_particle='h')
# fixed_postion is in reduced units
# in case, you want to plot hole density by fixing electron, set fix_particle = 'e'

## .cube will be dumped and use vesta to visualize it !

