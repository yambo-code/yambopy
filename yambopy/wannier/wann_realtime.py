import numpy as np
from yambopy.wannier.wann_utils import *
from yambopy.wannier.wann_dipoles import TB_dipoles
from yambopy.wannier.wann_ode import SecondOrderRungeKutta
from time import time
from scipy.integrate import solve_ivp

class RealTime():
    '''
    RealTime Simulations class
    needs h2p class to be given as an argumet
    Note that his class works with time in femtoseconds
    '''
    def __init__(self, h2p, t0, h,  model, E, eta) :
        self.h2p = h2p # BSE hamiltonian object
        self.E = E # electric field
        self.bse_nv = self.h2p.bse_nv
        self.bse_nc = self.h2p.bse_nc
        self.nb = self.h2p.nb
        self.model = model
        self.nc = self.h2p.nc
        self.nv = self.h2p.nv
        self.eta = eta # smearing
        self.dip_eta = eta
        self.f_kn = self.h2p.f_kn
        self.latdb = self.h2p.latdb
        self.t0 = t0/AU2FS
        self.h = h/AU2FS
        self.rho_nmk = self._get_rho_nmk()
        self.e_nmk = self._get_enmk()

        # I think for now I should give the normal dipoles not weighted by the BSE
        # if(self.h2p.nq != 1): 
        #     self.h2p.h2peigvec_vck=self.h2p.h2peigvec_vck[self.h2p.q0index]
        #     self.h2p.h2peigv_vck = self.h2p.h2peigv_vck[self.h2p.q0index]
        #     self.h2p.h2peigvec = self.h2p.h2peigvec[self.h2p.q0index]
        #     self.h2p.h2peigv = self.h2p.h2peigv[self.h2p.q0index]
        self.tb_dipoles = TB_dipoles(self.nc, self.nv,  self.bse_nc, self.bse_nv, self.h2p.nk, self.h2p.eigv, \
                                  self.h2p.eigvec, self.dip_eta, self.model.hlm, self.h2p.T_table, \
                                  self.h2p.BSE_table )
        self.d_knm = self.tb_dipoles.d_knm

    def _F_nmk(self, t, rho_nmk_t):
        tmp_sum1 = np.einsum('knpi, pmk -> nmki', self.d_knm, rho_nmk_t)
        tmp_sum2 = np.einsum('npk, kpmi -> nmki', rho_nmk_t, self.d_knm)
        F_nmk = np.einsum('i, nmkj -> nmk', self.E.E_t(t),(tmp_sum1-tmp_sum2))
        return F_nmk
    
    # Differential equation for rho_nm
    def drho_dt(self, t, rho_nmk_flat):
        rho_nmk_t = rho_nmk_flat.reshape(self.nb, self.nb, self.h2p.nk)
        F_nmk = self._F_nmk(t, rho_nmk_t)
        # Compute the full derivative array using broadcasting
        drho_dt_val = -1j * (self.e_nmk - 1j * self.eta) * rho_nmk_t - 1j * F_nmk
        return drho_dt_val.flatten()  # Flatten the derivatives before returning


    def _get_enmk(self):
        e_nmk = np.zeros((self.nb, self.nb, self.h2p.nk), dtype=np.complex128)
        e_nmk = self.h2p.eigv[:, :, np.newaxis] - self.h2p.eigv[:,np.newaxis,:]
        e_nmk = np.transpose(e_nmk, (1,2,0))
        return e_nmk/HA2EV # energies are given in eV, convert to Hartree
    
    def _get_rho_nmk(self):
        rho_nmk = np.zeros((self.nb, self.nb, self.h2p.nk), dtype = np.complex128)
        for k in range(self.h2p.nk):
            np.fill_diagonal(rho_nmk[...,k], self.f_kn[k,:] )
        return rho_nmk/ (self.latdb.lat_vol) # here density is number of occupations/direct volume in meter
    
    def _step(self, n, m, ik):
        k1 = self.h * self.drho_dt(self.t, n, m, ik)
        k2 = self.h * self.drho_dt(self.t + 0.5*self.h, n, m, ik)
        self.rho_nmk[n,m,ik] += k2
        self.t += self.h

    def solve(self, final_t):
        """
        Solves the ODE from the initial to the final time across all band indices.

        :param final_t: The final value of t for the simulation.
        :return: Updated rho_nmk with time-evolved values.
        """
        # Re-initialize for rerun
        self.rho_nmk = self._get_rho_nmk()
        final_t = final_t/AU2FS # working in a.u. of time
        
        ts = np.linspace(self.t0, final_t, int(np.floor((final_t-self.t0)/self.h)))
        # Flatten rho_nmk for y0
        y0_flat = self.rho_nmk.flatten()

        # Solve the ODE
        sol = solve_ivp(
            self.drho_dt,        # Function that returns the derivative
            t_span=(self.t0, final_t),  # Time span in atomic units
            y0=y0_flat,          # Initial condition
            t_eval=ts,           # Times at which to store the solution
            method='RK45'        # The integration method to use
        )

        # Reshape the result to the original shape
        sol_reshaped = sol.y.reshape(self.nb, self.nb, self.h2p.nk, -1)

        return sol.t * AU2FS, sol_reshaped*(self.latdb.lat_vol) # get in occupations and not as density over volume