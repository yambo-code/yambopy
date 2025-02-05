import numpy as np
from yambopy.wannier.wann_asegrid import ase_Monkhorst_Pack
from yambopy.wannier.wann_utils import *
from yambopy.wannier.wann_dipoles import TB_dipoles
from yambopy.wannier.wann_occupations import TB_occupations
from yambopy.dbs.bsekerneldb import *
from time import time
from yambopy.wannier.wann_H2p import H2P
class ExcBerry(H2P):
    def __init__(self, nk, nb, nc, nv,eigv, eigvec,bse_nv, bse_nc, T_table, latdb, kmpgrid, qmpgrid,excitons=None, \
                  kernel_path=None, excitons_path=None,cpot=None,ctype='v2dt2',ktype='direct',bsetype='resonant', method='model',f_kn=None, \
                  TD=False,  TBos=300): 
        super().__init__(nk, nb, nc, nv,eigv, eigvec,bse_nv, bse_nc, T_table, latdb, kmpgrid, qmpgrid,excitons=None, \
                  kernel_path=None, excitons_path=None,cpot=None,ctype='v2dt2',ktype='direct',bsetype='resonant', method='model',f_kn=None, \
                  TD=False,  TBos=300)
    
    
    
    
    
    def  find_surface_points(self, tolerance = 1e-6):
        tmp_qcar_points = self.qmpgrid.car_kpoints()
        idx_x0 = np.where(np.abs(tmp_qcar_points[:,0]) < tolerance)[0]
        idx_y0 = np.where(np.abs(tmp_qcar_points[:,1]) < tolerance)[0]
        idx_z0 = np.where(np.abs(tmp_qcar_points[:,2]) < tolerance)[0]
        surface_indices = {'x': idx_x0, 'y': idx_y0, 'z': idx_z0}
        self.surface_indices = surface_indices
    
    def _get_elec_overlap(self):
        Mmn = np.zeros((self.nb, self.nb,self.nk, self.nk), dtype=np.complex128)
        # here l stands for lambda, just to remember me that there is a small difference between lambda and transition index
        for n in range(self.nb):
            for m in range(self.nb):   
                for ik in range(self.nk):
                    for ikp in range(self.nk):
                        Mmn[n,m,ik, ikp] = np.vdot(self.eigvec[ik,:,n],self.eigvec[ikp,:,m])
        self.Mmn = Mmn