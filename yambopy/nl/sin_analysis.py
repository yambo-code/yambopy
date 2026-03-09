# Copyright (c) 2023-2025, Claudio Attaccalite,
#                          Myrta Gruening, Mao Yunchen
# All rights reserved.
#
# This file is part of the yambopy project
# Calculate linear response from real-time calculations (yambo_nl)
# Modified to allow generalisation
#
import numpy as np
from yambopy.units import ha2ev,fs2aut
from yambopy.nl.external_efield import Divide_by_the_Field
from tqdm import tqdm
import sys
import os
from abc import ABC,abstractmethod
from yambopy.nl.nl_analysis import Xn_from_signal
#
#
# Derived class for monochromatic signal
#    
class Xn_from_sine(Xn_from_signal):
        def set_defaults(self):
            EFIELDS = ["SIN","SOFTSIN"]
            if self.efield["name"] not in EFIELDS:
                raise ValueError(f"Invalid electric field for frequency mixing analysis. Expected one of: {EFIELDS}")
            for i_n in range(len(self.pumps)):
                if self.pumps[i_n]["name"] != 'none':
                    raise ValueError("This analysis is for one monochromatic field only.")
            if self.solver == '':
                self.solver = 'full'
            self.out_dim = self.X_order + 1
            if self.samp_mod== "":
                self.samp_mod="linear"
            return

        def set_sampling(self,ifrq):
            samp_order = 2*self.X_order + 1
            if self.nsamp == -1:
                self.nsamp = samp_order
            T_period = 2.0 * np.pi / self.freqs[ifrq]
            T_range, out_of_bounds = self.update_time_range(T_period)
            if (out_of_bounds):
                print(f'User range redifined for frequency {self.freqs[ifrq]* ha2ev:.3e} [eV]')
            print(f"Time range: {T_range[0] / fs2aut:.3f} - {T_range[1] / fs2aut:.3f} [fs]")
            return T_range
        
        def update_time_range(self,T_period): # not sure if this is a general or specific method - let it here for the moment
            T_range = self.T_urange
            out_of_bounds = False
            if T_range[0] <= 0.0:
                T_range[0] = self.time[-1] - T_period
            if T_range[1] > 0.0:
                T_range[1] = T_range[0] + T_period
            else:
                T_range[1] = self.time[-1]                
            if T_range[1] > self.time[-1]:
                T_range[1] = slef.time[-1]
                T_range[0] = T_range[1] - T_period
                out_of_bound = True
            return T_range, out_of_bounds

        def define_matrix(self,T_i,ifrq):
            M_size = len(T_i)
            M = np.zeros((M_size, M_size), dtype=np.cdouble)
            M[:, 0] = 1.0
            W = self.freqs[ifrq]
            for i_n in range(1, self.X_order+1):
                exp_neg = np.exp(-1j * i_n* W * T_i, dtype=np.cdouble)
                exp_pos = np.exp(1j * i_n * W * T_i, dtype=np.cdouble)
                M[:, i_n] = exp_neg
                M[:, i_n +self.X_order] = exp_pos
            return M

        def output_analysis(self,out,to_file=True):
            for i_order in range(self.X_order + 1):
                T = 10000.0
                for i_f in range(self.n_runs):
                    out[i_order, i_f, :] *= Divide_by_the_Field(self.efields[i_f], i_order)
                    T_period = 2.0 * np.pi / self.freqs[i_f]
                    Trange, _ = self.update_time_range(T_period)
                    T = min(T,Trange[0])
                out[i_order,:,:]*=self.get_Unit_of_Measure(i_order)
                run_info = self.append_runinfo(T)
                if (to_file):
                    output_file = f'o{self.prefix}.YamboPy-X_probe_order_{i_order}'
                    header = "E[eV] " + " ".join([f"X{i_order}/Im({d}) X{i_order}/Re({d})" for d in ('x','y','z')])
                    if self.l_out_current:
                        output_file = f'o{self.prefix}.YamboPy-Sigma_probe_order_{i_order}'
                        header = "E[eV] " + " ".join([f"S{i_order}/Im({d}) S{i_order}/Re({d})" for d in ('x','y','z')])
                    values = np.column_stack((self.freqs * ha2ev, out[i_order, :, 0].imag, out[i_order, :, 0].real,
                                      out[i_order, :, 1].imag, out[i_order, :, 1].real,
                                      out[i_order, :, 2].imag, out[i_order, :, 2].real))
                    np.savetxt(output_file, values, header=header, delimiter=' ', footer="Harmonic analysis results")
                    outf = open(output_file, "a")
                    outf.writelines(run_info)
                    outf.close()
                else:
                    print(run_info)
                    return (self.freqs, out)

        def reconstruct_signal(self,out,to_file=True):
            Seff = np.zeros((self.n_runs, 3, len(self.time)), dtype=np.cdouble)
            for i_f in tqdm(range(self.n_runs)):
                for i_d in range(3):
                    for i_order in range(self.X_order + 1):
                        freq_term = np.exp(-1j * i_order * self.freqs[i_f] * self.time)
                        Seff[i_f, i_d, :] += out[i_order, i_f, i_d] * freq_term
                        Seff[i_f, i_d, :] += np.conj(out[i_order, i_f, i_d]) * np.conj(freq_term)    
            for i_f in tqdm(range(self.n_runs)):
                values = np.column_stack((self.time / fs2aut, Seff[i_f, 0, :].real, Seff[i_f, 1, :].real, Seff[i_f, 2, :].real))
                output_file = f'o{self.prefix}.YamboPy-pol_reconstructed_F{i_f + 1}'
                header="[fs] Px Py Pz"
                if self.l_out_current:
                    output_file = f'o{self.prefix}.YamboPy-curr_reconstructed_F{i_f + 1}'
                    header="[fs] Jx Jy Jz"
                if (to_file):
                    np.savetxt(output_file, values, header=header, delimiter=' ', footer="Reconstructed signal")
                else:
                   return values

