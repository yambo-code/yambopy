# Copyright (c) 2025, Ignacio M. Alliati and Myrta Grüning
# All rights reserved.
#
# This file is part of the yambopy project
# Floquet analysis of the time-dependent bands from real-time Berry phase 
#
from yambopy import *
from yambopy.plot import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as pltcol
import sys
import os
from datetime import datetime
from yambopy.units import hbar_eVfs,ha2ev,fs2aut
#
# This class postprocess data from the ndb.V_bands and the ndb.V_bands_K_section databases 
# 
# 
class VbPP():

    def __init__(self,vb_db,kpt,band,vb_pars=None): 
        """vb_db: VBdb object (from YamboVbandsDB class)
           freq : float, field frequency
           vb_pars: list of floats
        """
        self.n_timesteps = vb_db.n_timesteps -1
        self.tvecs = self.read_tvecs(vb_db,kpt,band)
        times_= []
        for i in range(self.n_timesteps):
            times_.append(vb_db.times[i]/fs2aut) 
        self.times = times_ 
        self.kpt   = kpt
        self.band  = band
        self.basis_index = vb_db.basis_index
        self.basis_size = vb_db.basis_size
        self.ks_ev = get_bands(kpt=kpt,band=band)
        freq =  get_frequency(vb_db.vb_path)
        self.freq  =freq*ha2ev 
        self.period = 2.0*np.pi/freq/fs2aut #in fs
        self.qe_thrs,self.err_thrs,self.step,self.max_iter,self.max_fl_mode = self.set_conv_parameters(vb_db,vb_pars=None)
        self.tot_fl_modes = self.max_fl_mode * 2 + 1
        self.exp_mat_long,self.exp_mat_m1 = self.get_exp_matrix()

    def __str__(self):
        """
        Print all info of the class
        """
        s="\n * * * Vbands PP class  * * * \n\n"
        s+="Selected Kpt  : "+str(self.kpt)+"\n"
        s+="Selected Band : "+str(self.band)+"\n"
        s+="KS energy     : "+str(self.ks_ev)+" [eV]\n"
        s+="Basis index   : "+str(self.basis_index)+"\n"
        s+="\n =================================\n\n"
        s+="Field freq    : "+str(self.freq)+" [eV] \n"
        s+="Field period  : "+str(self.period)+" [fs] \n"
        s+="\n =================================\n\n"
        s+="QEnergy thrsh : "+str(self.qe_thrs)+"\n"
        s+="Error thrsh   : "+str(self.err_thrs)+"\n"
        s+="Conv step     : "+str(self.step)+"\n"
        s+="Max iteration : "+str(self.max_iter)+"\n"
        s+="\n =================================\n\n"
        s+="Max Floq mode : "+str(self.max_fl_mode)+"\n"
        s+="Tot Floq modes: "+str(self.tot_fl_modes)+"\n"
        s+="\n =================================\n\n"
        for i,v in enumerate(self.tvecs):
            s+="Time: "+str(self.times[i])+" fs:\n"
            s+='bnd   c.real             c.imag\n'
            for j in range(len(v)):
                idx = self.basis_index[0] + j
                s+=str(idx)+"   "+str(v[j].real)+"   "+str(v[j].imag)+"\n"
        return s

    def read_tvecs(self,vb_db,kpt,band):
        list_of_evecs = []
        for it in range(vb_db.n_timesteps-1):
            evec=np.zeros(vb_db.basis_size,dtype=complex)
            
            for i in range(vb_db.basis_size):
                evec[i] = vb_db.tvecs[it][kpt-1,band-1,i]
            list_of_evecs.append(evec)
        return list_of_evecs

    def set_conv_parameters(self,vb_db,vb_pars=None): # would be nice to be a dictionary
        if vb_pars is None:
            qe_thrs=1e-7
            err_thrs=1e-8
            step=0.01
            max_iter=300
            max_fl_mode = vb_db.fl_order
        else:
            qe_thrs,err_thrs,step,max_iter,max_fl_mode = list_pars    
            max_fl_mode = min(max_fl_mode,self.fl_order)
        list_pars = [qe_thrs,err_thrs,step,max_iter,max_fl_mode]
        
        return list_pars
        
    def get_exp_matrix(self):
        """uses variables of class to set up a call to the
           external build_exp_matrix
        """
        mat, mat_m1 = build_exp_matrix(
                 listof_times=self.times,
                 max_fl_mode=self.max_fl_mode,
                 freq=self.freq)
        return mat, mat_m1

# Processing part
    
    def calc_pvecs(self,qe_ev):
        """Returns the periodic part of the
           Floquet basis functions
        """
        mat_of_pvecs = np.zeros((len(self.tvecs),self.basis_size),dtype=complex)
        for i,v in enumerate(self.tvecs):
            mat_of_pvecs[i,:] = np.exp(+1j * qe_ev * self.times[i] / hbar_eVfs) * v

        return mat_of_pvecs


    def calc_fvecs(self,qe_ev,mat_of_pvecs=None):
        """this will calculate the Floquet vectors for a given qe_evx
           and return """
        if self.exp_mat_m1 is None:
            raise AttributeError("You need to initialize the FL space first")
        if mat_of_pvecs is None:
            mat_of_pvecs = self.calc_pvecs(qe_ev)

        mat_of_fvecs = np.matmul(self.exp_mat_m1,mat_of_pvecs[:self.tot_fl_modes,:])

        return mat_of_fvecs

    def recalc_pvecs_via_fl(self,qe_ev=None,mat_of_pvecs=None,mat_of_fvecs=None):
        """Function to calculate pVecs via the obtained fVecs at all times steps,
           including those not used to generate those fVecs, i.e., outside the
           first period considered. This allow to determine whether the qe used
           truly makes the tVecs periodic
        """
        if qe_ev is None and mat_of_pvecs is None:
            raise ValueError("Provide either qe or pVecs")
        if qe_ev is None and mat_of_fvecs is None:
            raise ValueError("Provide either qe or fVecs")
        if mat_of_pvecs is None:
            mat_of_pvecs = self.calc_pvecs(qe_ev)
        if mat_of_fvecs is None:
            mat_of_fvecs = self.calc_fvecs(qe_ev,mat_of_pvecs=mat_of_pvecs)

        mat_of_recalc_pvecs = np.matmul(self.exp_mat_long,mat_of_fvecs)

        return mat_of_recalc_pvecs


    def run_nl2fl(self,qe_ev=None,tag=None,iter_num=None): #tag to be reviewed
        """run #TODO
        """
        if qe_ev is None:
            qe_ev = self.ks_ev

        if tag is None:
            tag = datetime.today().strftime('%Y%m%d-%H.%M.%S')
            if iter_num is not None:
                tag += '_iter'+str(iter_num)

        fl_eigenvectors = FLeigenvectors(self.period,self.freq,self.max_fl_mode,self.times,tag)

        nl_in = self.calc_pvecs(qe_ev=qe_ev)
        fl_out = self.calc_fvecs(qe_ev,mat_of_pvecs=nl_in)
        nl_out = self.recalc_pvecs_via_fl(mat_of_pvecs=nl_in,mat_of_fvecs=fl_out)
        err = np.sum(np.abs(nl_out - nl_in))
        
        fl_eigenvectors.store_results(nl_in,nl_out,fl_out,qe_ev,err)

        return fl_eigenvectors

    def find_qe(self,qe_ev=None,tag=None):
        """ Secant solver to iterate over executions of run_nl2fl
            and minimize the error between nl_in and nl_out
        """
        if qe_ev   is None:
            qe_ev = self.ks_ev

        lof_qe  = []
        lof_err = []

        # Iteration with guess
        evecs = self.run_nl2fl(qe_ev,tag=tag,iter_num=0)
        lof_qe.append(evecs.FL_qe)
        lof_err.append(evecs.err)
        if evecs.err < self.qe_thrs:
            evecs.nr_it = 0
            evecs.nr_acc = self.qe_thrs
            return evecs
        # Iteration with step
        evecs = self.run_nl2fl(qe_ev-self.step,tag=tag,iter_num=1)
        lof_qe.append(evecs.FL_qe)
        lof_err.append(evecs.err)
        _delta_qe = abs(lof_qe[-1]-lof_qe[-2])

        # Loop
        iter_num = 1
        while ((_delta_qe > self.qe_thrs) or (lof_err[-1] > self.err_thrs)) and (iter_num <= self.max_iter):
            iter_num += 1
            _derivative = (lof_err[-1]-lof_err[-2])/(lof_qe[-1]-lof_qe[-2])
            _qe = lof_qe[-1] - lof_err[-1]/_derivative
            evecs = self.run_nl2fl(_qe,tag=tag,iter_num=iter_num)
            lof_qe.append(evecs.FL_qe)
            lof_err.append(evecs.err)
            _delta_qe = abs(lof_qe[-1]-lof_qe[-2])

        evecs.nr_it = iter_num
        evecs.nr_acc = abs(lof_qe[-1]-lof_qe[-2])

        return evecs
#
# This class handles Floquet eigenvectors
# 
# 
class FLeigenvectors:

    def __init__(self,T,f,max_eta,listof_times,tag):
        self.period = T
        self.freq = f
        self.max_fl_mode = max_eta
        self.tot_fl_modes = 2 * max_eta + 1
        self.times = np.array(listof_times)
        self.tag = tag
        self.dir = 'figs-'+str(self.tag)

    @classmethod
    def FromArrayOnly(cls,array,tag):
        cls.FL_vecs = array
        cls.dir = 'figs-'+str(tag)
        cls.tot_fl_modes = array.shape[0]
        cls.max_fl_mode = int((cls.tot_fl_modes-1)/2)
        return cls

    def store_results(self,NL_in,NL_out,FL_vecs,FL_qe,err):
        self.NL_in = NL_in
        self.NL_out = NL_out
        self.FL_vecs = FL_vecs
        self.FL_qe = FL_qe
        self.err = err

    def calc_all_times_pVecs(self,t_step=0.0025,n_steps=669):
        M,alltimes = build_exp_matrix(
                     t0_fs=self.times[0],
                     t_step=t_step,
                     n_steps=n_steps,
                     max_fl_mode=self.max_fl_mode,
                     freq=self.freq,
                     l_inv=False)

        NL_out_alltimes = np.matmul(M,self.FL_vecs)
        return alltimes,NL_out_alltimes

    def plot_realtime(self,band_to_plot=1,t_step=0.0025,n_steps=669):
        X,Y=self.calc_all_times_pVecs(t_step=t_step,n_steps=n_steps)

        fig,axes=plt.subplots(2)
        fig.set_size_inches(8.3,5.8)
        fig.suptitle(f'Time-dependent projection over Kohn-Sham state: {band_to_plot+1}')

        axes[0].plot(X,Y[:,band_to_plot].real,label='FL_calculated',color='tab:blue')
        axes[0].plot(self.times,self.NL_in[:,band_to_plot].real ,label='NL_in' ,marker='o',linestyle='',color='tab:olive',ms=12)
        axes[0].plot(self.times,self.NL_out[:,band_to_plot].real,label='NL_out',marker='.',linestyle='',color='tab:blue',ms=11)
        axes[1].plot(X,Y[:,band_to_plot].imag,label='FL_calculated',color='tab:blue')
        axes[1].plot(self.times,self.NL_in[:,band_to_plot].imag ,label='NL_in' ,marker='o',linestyle='',color='tab:olive',ms=12)
        axes[1].plot(self.times,self.NL_out[:,band_to_plot].imag,label='NL_out',marker='.',linestyle='',color='tab:blue',ms=11)
        axes[1].set_xlabel('Time (fs)')
        axes[0].set_ylabel(f'Re[ d_{band_to_plot+1} ]')
        axes[1].set_ylabel(f'Im[ d_{band_to_plot+1} ]')
        os.system(f'if [ ! -d {self.dir} ]; then mkdir {self.dir};fi')
        plt.legend()
        plt.savefig(f'{self.dir}/fig-real_time_projection_over_KS_state_{band_to_plot+1}.pdf')
        plt.close()

    def plot_floquet(self,labels='+1'):
        _fl_vec = self.FL_vecs.reshape((1,self.FL_vecs.size),order='F')
        fig = plt.figure()
        fig.set_size_inches(8.3,3.9)
        ax = fig.gca()
        fig.suptitle('Projection over Floquet-Kohn-Sham space')

        cax = ax.matshow(abs(_fl_vec), cmap='YlOrBr', norm=pltcol.LogNorm(vmin=1.e-10,vmax=1.))
        fig.colorbar(cax,orientation='horizontal')

        lof_labels=[]
        for band in range(self.FL_vecs.shape[1]):
            for i in range(self.FL_vecs.shape[0]):
                tic='        '+str(i-self.max_fl_mode)
                lof_labels.append(tic)
        ax.set_xticks([x-0.5 for x in range(_fl_vec.shape[1])])
        ax.set_xticklabels(lof_labels)
        ax.set_yticks([])
        ax.set_yticklabels([])

        ax.text(0,-2.3, f'{abs(self.FL_vecs[self.max_fl_mode+1,0]):.1e}', fontsize=15)
        ax.text(4,-2.3, f'{abs(self.FL_vecs[self.max_fl_mode,0]):.1e}', fontsize=15)
        ax.text(10,-2.3, f'{abs(self.FL_vecs[self.max_fl_mode-1,1]):.1e}', fontsize=15)
        ax.text(13.5,-2.3, f'{abs(self.FL_vecs[self.max_fl_mode+1,1]):.1e}', fontsize=15)
        if (labels == '+2'):
            ax.text(16.5,-2.3, f'{abs(self.FL_vecs[self.max_fl_mode+2,1]):.1e}', fontsize=15)
        plt.grid(which='major',lw=1.,color='black')

        os.system(f'if [ ! -d {self.dir} ]; then mkdir {self.dir};fi')
        plt.savefig(f'{self.dir}/fig-FKS_projection.pdf')
        plt.close()

    def output(self):
        os.system(f'if [ ! -d {self.dir} ]; then mkdir {self.dir};fi')
        with open(f'{self.dir}/output-{self.tag}.dat','w') as f:
            for k,v in self.__dict__.items():
                f.write('--------------------------\n')
                f.write(f'Attribute: {k}\n')
                f.write('Values:\n')
                f.write(str(v)+'\n')

    def output_for_fortran(self,band,kpt,NL_band_1,file='fortran_input.txt'):
        if kpt == 1:
            filemode='w'
        else:
            filemode='a'
        with open(f'{file}',filemode) as f:
            for state in range(self.FL_vecs.shape[1]):
                for mode in range(self.FL_vecs.shape[0]):
                    _state = state + NL_band_1
                    f.write(f'FL_V_bands({_state},{mode+1},{band},{kpt},1) = {self.FL_vecs[mode,state].real}_SP + cI * {self.FL_vecs[mode,state].imag}_SP\n')
            f.write(f'FL_QE({band},{kpt},1) = {self.FL_qe/ha2ev}_SP\n')

##################################################################################
### AUXILIARY FUNCTIONS ###
##################################################################################

def build_exp_matrix(listof_times=None,t0_fs=None,t_step=None,n_steps=None,max_fl_mode=None,freq=None,l_inv=True):
  """builds matrix of exponentials with times
     and size given by self.times and self.tot_fl_modes
  """
  if max_fl_mode is None:
    raise ValueError("Missing total number of FL modes on input to build_exp_matrix")
  if listof_times is None and (t0_fs is None or t_step is None or n_steps is None):
    raise ValueError("Missing times on input to build_exp_matrix")

  if listof_times is None:
    listof_times = []
    for step in range(-int(n_steps*0.06),n_steps):
      t = t0_fs + t_step * step
      listof_times.append(t)

  tot_fl_modes = 2 * max_fl_mode + 1
  matrix = np.zeros((len(listof_times),tot_fl_modes),dtype=complex)
  for i,t in enumerate(listof_times):  # rows - time
   for j in range(tot_fl_modes): # columns - eta
     eta = j - max_fl_mode
     matrix[i,j]=np.exp(-1j * eta * (freq/hbar_eVfs) * t)

  if l_inv:
    inverse = np.linalg.inv(matrix[:tot_fl_modes,:])
    return matrix, inverse
  else:
    arrof_time = np.array(listof_times)
    return matrix, arrof_time

def get_frequency(db_path):
   ds=Dataset(db_path+'/ndb.Nonlinear')
   freq = float(ds['Field_Freq_range_1'][0])
   return freq

def get_bands(kpt=None,band=None):
  ds=Dataset('SAVE/ns.db1')
  KSevalues=ds['EIGENVALUES'][:]
  vb,vb_kpt,vbM=0,0,-111
  for _band in range(KSevalues.shape[2]):
    for _kpt in range(KSevalues.shape[1]):
      for _spin in range(KSevalues.shape[0]):
        if KSevalues[_spin,_kpt,_band] < 0. and KSevalues[_spin,_kpt,_band] > vbM:
           vb,vb__kpt,vbM = _band,_kpt,KSevalues[_spin,_kpt,_band]

  if kpt is None and band is None:
    return (KSevalues - vbM)*ha2ev
  else:
    return (KSevalues[0,kpt-1,band-1] - vbM)*ha2ev
