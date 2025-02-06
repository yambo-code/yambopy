import os
import numpy as np
from yambopy.dbs.excitondb import YamboExcitonDB
from yambopy.dbs.latticedb import YamboLatticeDB
from yambopy.dbs.wfdb import YamboWFDB
from .exciton_matrix_elements import exciton_X_matelem

def compute_exciton_spin(path='.',bse_dir='SAVE',iqpt=1, nstates=-1, 
                        sz = 0.5*np.array([[1, 0], [0, -1]]) ):
    ## Computing <S'|S_z|S>
    # Note that it also computes off-diagonal.
    ## One needs to diagonalize in degnerate subspace to get the spin of excitons 
    filename = 'ndb.BS_diago_Q%d'%(iqpt)
    #
    lattice  = YamboLatticeDB.from_db_file(os.path.join(path,'SAVE','ns.db1'))
    #
    excdb = YamboExcitonDB.from_db_file(lattice,filename='ndb.BS_diago_Q1',
                                        folder=os.path.join(path,bse_dir),
                                        Load_WF=True, neigs=nstates)
    #
    wfdb = YamboWFDB(path=path,bands_range=[np.min(excdb.table[:,1])-1,
                                            np.max(excdb.table[:,2])]  )
    #
    print(np.min(excdb.table[:,1])-1,np.max(excdb.table[:,2]))
    assert wfdb.nspinor == 2, "Makes sense only for nspinor = 2"
    #
    excdb.convert_to_kcv()
    #
    elec_sz = wfdb.get_spin_m_e_BZ(s_z=sz)
    #
    excQpt = excdb.car_qpoint
    # convert to crystal units 
    excQpt = lattice.lat@excQpt
    #
    exe_Sz = exciton_X_matelem(excQpt, np.array([0,0,0]), excdb.eigenvectors,
                               excdb.eigenvectors, elec_sz[None,...], wfdb.kBZ,
                               diagonal_only=False)
    print("Note : This is a spin matrix. Diagonalize the matrix in degenerate subspace to get spin values.")
    return exe_Sz
