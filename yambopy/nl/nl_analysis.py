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
import datetime
from abc import ABC,abstractmethod
#
#
# Template class
#
class Xn_from_signal(ABC):
    def __init__(self,nldb,X_order=4,T_range=[-1, -1],l_out_current=False,nsamp=-1,samp_mod='',solver='',tol=1e-10,debug_mode=False,
                 current_type='berry'):
        """
        current_type (used only with l_out_current=True), length gauge:
           'berry'   : NL_J(1:3), Berry-phase current of Souza et al. PRB 69, 085106 (finite differences in k
                       of the time-dependent valence states, NL_Berry_current.F)
           'dipoles' : NL_J(4:6), J = Tr[rho v] with the velocity matrix elements DIP_v (no k-derivative).
                       Recommended for the DC (order 0, shift current): on coarse grids the k-derivative of
                       the resonant first-order states makes the Berry DC current converge very slowly.
        """
        CURRENT_TYPES = {'berry':0, 'dipoles':3}
        if current_type not in CURRENT_TYPES:
            raise ValueError(f"Invalid current_type. Expected one of: {list(CURRENT_TYPES)}")
        self.current_type = current_type
        self.i_curr0 = CURRENT_TYPES[current_type]   # first row of NL_J_t (6 rows: berry x,y,z, dipoles x,y,z)
        self.nldb = nldb                # kept for the description of the simulation in the output files
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
            s+="Output is set to current ("+self.current_type+")."
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
            S_i = np.array([self.current[ifrq][idir+self.i_curr0,i] for i in i_t])
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

    # ------------------------------------------------------------------
    # Description of the output files (header) and of the simulation (footer)
    # ------------------------------------------------------------------
    @staticmethod
    def _order_name(n):
        names={0:"DC / rectification (order 0)",1:"linear response (order 1)",2:"second harmonic (order 2)",
               3:"third harmonic (order 3)"}
        return names.get(n,f"harmonic of order {n}")

    def response_units(self,n):
        """Units of the order-n output, consistent with get_Unit_of_Measure."""
        if self.l_out_current:
            if n==0: return "A/V^2"
            if n==1: return "S/m (= A/m^2 / (V/m))"
            return f"A/m^2 / (V/m)^{n}"
        if n==0: return "esu (cm/statV)"
        if n==1: return "dimensionless (a.u.)"
        if n==2: return "esu (cm/statV)"
        return f"esu (cm/statV)^{n-1}"

    def response_description(self,n):
        """Name, definition and notes of the response function extracted at order n (one monochromatic field)."""
        v=np.asarray(self.efield["versor"],dtype=float)
        vlab=self._versor_label(v)
        if self.l_out_current:
            src=f"current J(t) ({self.current_type} current: " + \
                ("NL_J(1:3), Berry-phase/length-gauge current" if self.current_type=='berry' else "NL_J(4:6), Tr[rho v] with the dipoles DIP_v") + ")"
            if n==0:
                name="sigma^(2)(0;w,-w) - DC photocurrent (shift + injection-type, bulk photovoltaic)"
                defi="J_DC = sigma(0;w,-w)E(w)E(-w) + sigma(0;-w,w)E(-w)E(w) = 2 sigma |E(w)|^2"
                note=["Output S_0 = -J_DC/|E(w)|^2 = -2 sigma (factor -4/E0^2, see Divide_by_the_Field):",
                      "multiply by -1/2 to get sigma in the convention of Lumen nl_response (Shift_Current.F)."]
            elif n==1:
                name="sigma^(1)(w) - linear optical conductivity"; defi="J(w) = sigma(w) E(w)"; note=[]
            elif n==2:
                name="sigma^(2)(-2w;w,w) - second-harmonic current"; defi="J(2w) = sigma(-2w;w,w) E(w)E(w)"; note=[]
            else:
                name=f"sigma^({n})(-{n}w;w,...,w) - harmonic current of order {n}"; defi=f"J({n}w) = sigma E(w)^{n}"; note=[]
        else:
            src="polarization P(t) (Berry-phase polarization, NL_P)" if not self.nldb.l_use_DIPOLES else "polarization P(t) (from the dipoles)"
            if n==0:
                name="chi^(2)(0;w,-w) - optical rectification"
                defi="P_DC = chi(0;w,-w)E(w)E(-w) + chi(0;-w,w)E(-w)E(w) = 2 chi |E(w)|^2"
                note=["Output X_0 = -P_DC/|E(w)|^2 = -2 chi (factor -4/E0^2, see Divide_by_the_Field)."]
            elif n==1:
                name="chi^(1)(w) - linear susceptibility"; defi="P(w) = chi(w) E(w)"
                note=["Dielectric function: eps(w) = 1 + 4 pi chi(w)."]
            elif n==2:
                name="chi^(2)(-2w;w,w) - second-harmonic generation (SHG)"; defi="P(2w) = chi(-2w;w,w) E(w)E(w)"
                note=["1 esu (cm/statV) = 4.189e-4 m/V = 4.189e8 pm/V  (SI chi = 4 pi/(3e4) x esu)."]
            elif n==3:
                name="chi^(3)(-3w;w,w,w) - third-harmonic generation (THG)"; defi="P(3w) = chi(-3w;w,w,w) E(w)^3"; note=[]
            else:
                name=f"chi^({n})(-{n}w;w,...,w) - harmonic of order {n}"; defi=f"P({n}w) = chi E(w)^{n}"; note=[]
        lines=[ "YamboPy: nonlinear response from a real-time simulation (Lumen / yambo_nl)",
                f"Quantity    : {name}",
                f"Definition  : {defi}",
                f"Field       : E(t) = E0 sin(w t) e_v, E(w) = i E0/2 (coefficient of e^(-iwt)), e_v = {vlab}",
                f"Signal      : {src}",
                f"Tensor      : column 'a' is the Cartesian component a of the response to the field along e_v,",
                f"              i.e. the tensor element (a; v,...,v), e.g. X^(2)_a = chi^(2)_(a v v)",
                f"Units       : {self.response_units(n)}" ]
        lines.append("Normalization: per unit volume of the simulation cell (2D systems: rescale by L_z/t_eff)")
        for x in note: lines.append("Note        : "+x)
        return lines

    @staticmethod
    def _versor_label(v):
        lab=[c for c,x in zip('xyz',v) if abs(x)>1e-6]
        if len(lab)==1 and abs(abs(v['xyz'.index(lab[0])])-1)<1e-6: return f"{lab[0]} {np.round(v,4).tolist()}"
        return str(np.round(v,4).tolist())

    def make_header(self,description,qlabel):
        """Header for np.savetxt: description lines, column list and the column-name line."""
        cols=["Columns:","   1  E [eV]   photon energy of the field (hbar w)"]
        ic=2
        for d in ('x','y','z'):
            cols.append(f"  {ic:2d}  Im {qlabel}_{d}"); ic+=1
            cols.append(f"  {ic:2d}  Re {qlabel}_{d}"); ic+=1
        names="E[eV] "+" ".join([f"Im({qlabel}_{d}) Re({qlabel}_{d})" for d in ('x','y','z')])
        return "\n".join(description+[""]+cols+["",names])

    def append_runinfo(self,T):
        """Footer: simulation and analysis parameters."""
        nl=self.nldb
        def fld(e,tag):
            s =f"# {tag}\n"
            s+=f"#   Type                  : {e['name']}\n"
            s+=f"#   Intensity [a.u.]      : {e['intensity']:.6e}\n"
            s+=f"#   Amplitude E0          : {e['amplitude']:.6e} [a.u.]  = {e['amplitude']*AU2VMm1:.6e} [V/m]\n"
            s+=f"#   Versor                : {np.round(np.asarray(e['versor'],dtype=float),6).tolist()}\n"
            s+=f"#   Frequency (1st run)   : {e['freq_range'][0]*ha2ev:.6f} [eV]\n"
            s+=f"#   Initial time          : {e['initial_time']/fs2aut:.4f} [fs]\n"
            if e['name'] in ('QSSIN','GAUSS','PULSE','SPULSE'):
                s+=f"#   Width/damping         : {e['damping']/fs2aut:.4f} [fs]\n"
            return s
        s ="#\n# ===================== Simulation parameters (ndb.Nonlinear) =====================\n"
        s+=f"# Database                : {nl.nl_path}\n"
        s+=f"# Gauge                   : {nl.Gauge}\n"
        s+=f"# Correlation             : {nl.Correlation}\n"
        s+=f"# Integrator              : {nl.Integrator}\n"
        s+=f"# Bands                   : {list(np.asarray(nl.RT_bands).tolist())}\n"
        s+=f"# Damping (dephasing)     : {nl.NL_damping*ha2ev:.6f} [eV]\n"
        rad = "= dephasing" if nl.RAD_LifeTime<=0 else f"{nl.RAD_LifeTime/fs2aut:.4f} [fs]"
        s+=f"# Radiative life-time     : {rad}\n"
        s+=f"# Time step               : {nl.RT_step/fs2aut:.6f} [fs]\n"
        s+=f"# Time steps / total time : {nl.NE_steps} / {nl.NE_steps*nl.RT_step/fs2aut:.4f} [fs]\n"
        if nl.n_frequencies>0:
            s+=f"# Frequencies             : {nl.n_frequencies} in [{nl.NL_er[0]*ha2ev:.4f}, {nl.NL_er[1]*ha2ev:.4f}] [eV]\n"
        if nl.n_angles>0:
            s+=f"# Angles                  : {nl.n_angles}\n"
        s+=f"# Use dipoles (no Berry)  : {nl.l_use_DIPOLES}\n"
        s+=f"# Force 2nd order k-deriv.: {nl.l_force_SndOrd}\n"
        s+=f"# Current evaluated       : {nl.l_eval_CURRENT}\n"
        if nl.Correlation not in ('IPA',):
            s+=f"# G-vectors Hartree / Sx  : {nl.QP_ng_SH} / {nl.QP_ng_Sx}\n"
        s+=fld(self.efield,"Probe field")
        for i_p,e in enumerate(self.pumps):
            if e['name']!='none': s+=fld(e,f"Pump field {i_p+1}")
        s+="# ========================== Analysis parameters ==========================\n"
        s+=f"# Analysis                : {type(self).__name__}\n"
        s+=f"# Signal                  : {'current ('+self.current_type+')' if self.l_out_current else 'polarization'}\n"
        s+=f"# Max harmonic order      : {self.X_order}\n"
        s+=f"# Sampling mode           : {self.samp_mod}\n"
        s+=f"# Solver                  : {self.solver}\n"
        s+=f"# Sampling points         : {self.nsamp}\n"
        s+=f"# Start sampling time     : {T/fs2aut:.4f} [fs]\n"
        if self.T_urange!=[-1, -1]:
            s+=f"# User time range         : {self.T_urange} [a.u.]\n"
        s+=f"# Written                 : {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
        return s
        

### some maths auxiliary functions
    def IsSquare(self,m):
        return m.shape[0] == m.shape[1]

    def IsWellConditioned(self,m): # with this I am trying to avoid inverting a matrix
        return np.linalg.cond(m) < 1/sys.float_info.epsilon

    def residuals_func(x,M,S_i):
        x_cmplx=x[0:int(x.size/2)] + 1j * x[int(x.size/2):x.size]
        return np.linalg.norm(np.dot(M, x_cmplx) - S_i)
