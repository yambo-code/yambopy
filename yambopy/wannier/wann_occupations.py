import numpy as np
from yambopy.wannier.wann_utils import fermi_dirac_T, fermi_dirac

class TB_occupations():
    '''Class that handles occupations. The eigenvalues are assumed to be sorted'''
    def __init__(self, eigv, Tel, Tbos, elphg=0.0, Eb=0, fermie=0.0):
        self.eigv = eigv
        self.fermie = fermie
        self.Tel = Tel
        self.Tbos = Tbos
        self.Eb = Eb
        #self.elphg = self.elphg
    
    def _get_fkn(self, method):
        nk = self.eigv.shape[0]
        nb = self.eigv.shape[1]
        f_kn = np.zeros((nk,nb), dtype=np.float64)
        if (method =='FD'):
            if (self.Tel==0):
                f_kn = fermi_dirac(self.eigv, self.fermi)
            else:
                f_kn = fermi_dirac_T(self.eigv, self.T, self.fermie)
        if (method == 'Boltz'):
            kb = 8.617333262*10**-5
            f_kn = np.exp(-self.Eb/(kb*self.Tbos))   
            return  f_kn     
