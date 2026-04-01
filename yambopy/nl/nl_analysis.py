# Copyright (c) 2023-2026, Claudio Attaccalite,
#                          Myrta Gruening
# All rights reserved.
#
# This file is part of the yambopy project
# Calculate linear response from real-time calculations (yambo_nl)
# Modified to allow generalisation
#
import numpy as np
from yambopy.units import ha2ev,fs2aut, Junit, EFunit,SVCMm12VMm1,AU2VMm1
from yambopy.nl.external_efield import Divide_by_the_Field
from tqdm import tqdm
from scipy.optimize import least_squares
import sys
import os
from abc import ABC,abstractmethod
#
#
# Template class
#
class Xn_from_signal(ABC):
    def __init__(self,nldb,X_order=4,T_range=[-1, -1],l_out_current=False,nsamp=-1,samp_mod='',solver='',tol=1e-10,debug_mode=False): 
        self.time = nldb.IO_TIME_points # Time series
        self.T_step =self.time[1] - self.time[0]     # Time step of the simulation
        self.T_deph =12.0/nldb.NL_damping 
        self.efield = nldb.Efield[0] # External field of the first run
        self.pumps = nldb.Efield_general[1:] 
        self.efields = nldb.Efield
        self.n_runs = len(nldb.Polarization)     # Number of external laser frequencies
        self.polarization = nldb.Polarization     # Array of polarizations for each laser frequency
        self.current     =nldb.Current # Array of currents for each laser frequency
        self.l_eval_current=nldb.l_eval_CURRENT
        self.l_out_current = l_out_current and self.l_eval_current
        self.X_order = X_order
        SOLVE_MODES = ['', 'full', 'lstsq', 'lstsq_opt', 'svd']
        if solver not in SOLVE_MODES: 
            raise ValueError("Invalid solver mode. Expected one of: %s" % SOLVE_MODES)
        self.solver = solver
        SAMP_MODES = {'','linear', 'log', 'random'}    
        if samp_mod not in SAMP_MODES:
            raise ValueError(f"Invalid sampling mode. Expected one of: {SAMP_MODES}")
        self.samp_mod = samp_mod
        self.nsamp = nsamp
        self.T_urange = T_range 
        self.freqs = np.array([efield["freq_range"][0] for efield in nldb.Efield], dtype=np.double)     # Harmonic frequencies
        self.prefix = f'-{nldb.calc}' if nldb.calc != 'SAVE' else ''
        self.out_dim = 0
        self.tol = tol
        self.debug_mode = debug_mode
        super().__init__()
        
    def __str__(self): 
        """
        Print info of the class
        """
        s="\n * * *  Xn from signal class  * * * \n\n"
        s+="Max time: "+str(self.time[-1])+"\n"
        s+="Time step : "+str(self.T_step)+"\n"
        s+="Type Efield    : "+str(self.efield["name"])+"\n"
        s+="Number of runs   : "+str(self.n_runs)+"\n"
        s+="Current evaluated    : "+str(self.l_eval_current)+"\n"
        s+="Max harmonic order   : "+str(self.X_order)+"\n"
        if self.samp_mod != '':
            s+="Sampling mode    : "+str(self.samp_mode)+"\n"
        if self.solver != '':
            s+="Solver           : "+str(self.solver)+"\n"
        if self.nsamp > 0:
            s+="Sampling points  : "+str(self.nsamp)+"\n"
        if self.T_urange!=[-1, -1]:
            s+="User time range      : "+str(self.T_urange)+" [au] \n"
        s+="Frequency range: ["+str(self.freqs[0])+","+str(self.freqs[-1])+"] [au] \n"
        if  self.l_out_current:
            s+="Output is set to current."
        else:
            s+="Output is set to polarization."
        return s
    
    @abstractmethod
    def set_defaults(self):
        pass

    @abstractmethod
    def set_sampling(self,ifrq): 
       pass
   
    @abstractmethod
    def define_matrix(self,samp,ifrq):
        pass

    @abstractmethod
    def update_time_range(self):
        pass

    @abstractmethod
    def get_Unit_of_Measure(self,i_order):
        pass

    def get_sampling(self,T_range,idir,ifrq): 
        i_t_start = int(np.round( T_range[0] / self.T_step)) 
        i_deltaT  = int(np.round((T_range[1]-T_range[0])/self.T_step)/self.nsamp)
        if self.samp_mod=='linear':
            i_t = i_t_start + i_deltaT * np.arange(self.nsamp)
        elif self.samp_mod=='log':
            T_i = np.geomspace(i_t_start * self.T_step, T_range[1], self.nsamp, endpoint=False)
            i_t = np.array(np.round(T_i/self.T_step),dtype=int)
        elif self.samp_mod=='random':
            i_t = np.random.uniform(i_t_start,int(np.round(t / self.T_step)), self.nsamp)
        T_i = i_t*self.T_step - self.efield["initial_time"] 
        if self.l_out_current:
            S_i = np.array([self.current[ifrq][idir,i] for i in i_t])
        else:
            S_i = np.array([self.polarization[ifrq][idir,i] for i in i_t]) 
        return T_i,S_i
    
    def solve_lin_system(self,mat,samp,init=None):
        mat_dim = int(mat.shape[1])
        out=np.zeros(mat_dim,dtype=np.cdouble)
        if self.solver=="full" and ((not self.IsSquare(mat)) or (not self.IsWellConditioned(mat))):
            print(f'WARNING: solver changed to least square since square:{self.IsSquare(mat)} and well-conditioned:{self.IsWellConditioned(mat)}')
            self.solver = "lstsq"
        if self.solver=="full":
            out = np.linalg.solve(mat,samp)
        if self.solver=="lstsq":
            out = np.linalg.lstsq(mat,samp,rcond=self.tol)[0]
        if self.solver=="svd":
            inv = np.linalg.pinv(mat,rcond=self.tol)
            for i_n in range(mat_dim):
                out[i_n]=out[i_n]+np.sum(inv[i_n,:]*samp[:])
        if self.solver=="lstsq_opt":
            if(init is None):
                x0_cmplx = np.linalg.lstsq(mat, samp, rcond=tol)[0]
            else:
                x0_cmplx = init
            x0 = np.concatenate((x0_cmplx.real, x0_cmplx.imag))
            res = least_squares(residuals_func, x0, ftol=1e-11,gtol=1e-11,xtol=1e-11,verbose=0,x_scale='jac',args=(mat,samp))
            out = res.x[0:int(res.x.size/2)] + 1j * res.x[int(res.x.size/2):res.x.size]
        return out
    
    def perform_analysis(self):
        _ = self.set_defaults()
        out = np.zeros((self.out_dim, self.n_runs, 3), dtype=np.cdouble) 
        for i_f in tqdm(range(self.n_runs)):
            T_r = self.set_sampling(i_f)
            for i_d in range(3):
                samp_time, samp_sig= self.get_sampling(T_r,i_d,i_f)
                matrix = self.define_matrix(samp_time,i_f)
                raw = self.solve_lin_system(matrix,samp_sig)
                out[:, i_f, i_d] = raw[:self.out_dim]
                if self.debug_mode:
                    print(f"Freq #{i_f}, direction: {i_d}:")
                    print("***Sampling:")
                    print(samp_time, samp_sig)
                    print("***Matrix:")
                    print(matrix)
                    print("***Solution:")
                    print(raw)
        return out

    def get_Unit_of_Measure(self,i_order): # not sure if this is a general or specific method - let it here for the moment
        linear = 1.0
        ratio = SVCMm12VMm1 / AU2VMm1 # From AU to statVolt/cm ...is there a better way than this?
        if self.l_out_current:
            linear = Junit/EFunit
            ratio = 1.0/EFunit
        if i_order == 0:
            return np.power(ratio, 1, dtype=np.float64)*linear
        return np.power(ratio, i_order - 1, dtype=np.float64)*linear

    @abstractmethod
    def output_analysis(self,out,to_file):
        pass

    @abstractmethod
    def reconstruct_signal(self,out,to_file):
        pass

    def append_runinfo(self,T):
        s = "# Field details:"
        s+= "# Type field            "+str(self.efield["name"])+"\n"
        s+= "# Field intensity       "+str(self.efield["intensity"])+"\n"
        s+= "# Field versor          "+str(self.efield["versor"])+"\n"
        if self.efield["name"]=='QSSIN':
            s+= "Centre Gaussian     "+str(self.T0/fs2aut)+"[fs] \n"
            s+= "Gaussian sigma      "+str(self.sigma/fs2aut)+"[fs] \n"
        s+= "# Analysis details:"
        s+="# Max harmonic order   : "+str(self.X_order)+"\n"
        s+="# Sampling mode        : "+str(self.samp_mod) +"\n"
        s+="# Solver               : "+str(self.solver)+"\n"
        s+="# Sampling points      : "+str(self.nsamp)+"\n"
        s+="# Start sampling time  : "+str(T/fs2aut)+" [fs] \n"
        return s
        

### some maths auxiliary functions
    def IsSquare(self,m):
        return m.shape[0] == m.shape[1]

    def IsWellConditioned(self,m): # with this I am trying to avoid inverting a matrix
        return np.linalg.cond(m) < 1/sys.float_info.epsilon

    def residuals_func(x,M,S_i):
        x_cmplx=x[0:int(x.size/2)] + 1j * x[int(x.size/2):x.size]
        return np.linalg.norm(np.dot(M, x_cmplx) - S_i)
