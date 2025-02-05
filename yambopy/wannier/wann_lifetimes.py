import numpy as np
from yambopy.wannier.wann_dipoles import TB_dipoles
from yambopy.wannier.wann_utils import HA2EV

class TB_lifetimes(TB_dipoles):
    ''' compute the lifetimes of excitons for 3D, 2D and 1D systems
    $$
    \tau_{3 D, \alpha, \beta}^n=\frac{3 c^2 \hbar^3 N_{\mathbf{k}}}{4 \chi\left(E_0^n\right)^3 F_{\alpha, \beta}^{n, B S E}}
    $$
    $$
    \tau_{2 D, \alpha, \beta}^n=\frac{\hbar A_{u c} N_{\mathbf{k}}}{8 \pi \chi E_0^n F_{\alpha, \beta}^{n, B S E}}
    $$
    $$
    \tau_{1 D, \alpha, \beta}^n=\frac{c \hbar^2 l_1 N_{\mathbf{k}}}{2 \pi \chi\left(E_0^n\right)^2 F_{\alpha, \beta}^{n, B S E}}
    $$
    '''
    
    def __init__(self , dim, h2peigv, latdb, ntransitions, nc, nv, nkpoints, eigv, eigvec, \
                 eta, hlm, T_table, h2peigvec = None, method = 'real'):
        super().__init__(ntransitions, nc, nv, nkpoints, eigv, eigvec, \
                 eta, hlm, T_table, h2peigvec = None, method = 'real')
        self.dim = dim
        self.latdb = latdb
        try:
            self.h2peigvec = h2peigvec
        except ValueError:
            print('\nError! Before computing lifetimes you need the excitonic energies, therefore to solve the BSE\n')
        if (self.dim == '3D'):
            self.tau3D = self._get_tau3D()
        elif(self.dim == '2D'):
            self.tau2D = self._get_tau2D()
        elif(self.dim == '1D'):
            self.tau1D = self._get_tau1D()
    
    @classmethod
    def _get_tau3D(cls):
        tau_3D = np.zeros(cls.ntransitions,3,3)
        F_kcv = cls.F_kcv
        #h2peigvec are in eV-> convert in atomic units for formula, return seconds
        h2peigvec = cls.h2peigvec/HA2EV
        for t in range(0,cls.ntransitions):
            tau_3D[t,:,:] = 3*cls.nkpoints/(4*cls.h2peigvec[t]**3*F_kcv[t,:,:])
        
        return tau_3D

    @classmethod
    def _get_tau2D(cls):
        tau_2D = np.zeros(cls.ntransitions,3,3)
        F_kcv = cls.F_kcv
        # compute area of unit cell (assuming in plane are first and second lattice vector)
        vc = np.linalg.norm(np.cross(cls.latdb.alat[0],cls.latdb.alat[1]))
        #h2peigvec are in eV-> convert in atomic units for formula, return seconds
        h2peigvec = cls.h2peigvec/HA2EV
        for t in range(0,cls.ntransitions):
            tau_2D[t,:,:] = vc*cls.nkpoints/(8*np.pi*cls.h2peigvec[t]*F_kcv[t,:,:])
        
        return tau_2D
    
    @classmethod
    def _get_tau1D(cls):
        tau_1D = np.zeros(cls.ntransitions,3,3)
        F_kcv = cls.F_kcv
        # compute length of 1D system (assuming lattice vector to be the relevant one)
        vc = np.linalg.norm(cls.latdb.alat[0])
        #h2peigvec are in eV-> convert in atomic units for formula, return seconds
        h2peigvec = cls.h2peigvec/HA2EV
        for t in range(0,cls.ntransitions):
            tau_1D[t,:,:] = vc*cls.nkpoints/(2*np.pi*cls.h2peigvec[t]**2*F_kcv[t,:,:])
        
        return tau_1D
