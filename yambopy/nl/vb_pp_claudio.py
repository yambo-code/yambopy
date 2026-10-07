# Copyright (c) 2026, C. Attaccalite
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
# All units are in a.u. 
# 
class VbPP_cla():

    def __init__(self,vb_db,nl_db,kpt,band,vb_pars=None): 
        """vb_db: VBdb object (from YamboVbandsDB class)
           freq : float, field frequency
           vb_pars: list of floats
        """
        self.n_timesteps = vb_db.n_timesteps -1
        self.tvecs = self.read_tvecs(vb_db,kpt,band)
        times_= []
        for i in range(self.n_timesteps):
            times_.append(vb_db.times[i]) 
        self.times = times_ 
        self.t_step= nl_db.RT_step
        self.kpt   = kpt
        self.band  = band
        self.basis_index = vb_db.basis_index
        self.basis_size = vb_db.basis_size
        self.ks_ev = get_bands(kpt=kpt,band=band)
        freq =  get_frequency(vb_db.vb_path)
        self.freq  =freq
        self.period = 2.0*np.pi/freq #in a.u.
        self.qe_thrs,self.err_thrs,self.step,self.max_iter,self.max_fl_mode = self.set_conv_parameters(vb_db,vb_pars=vb_pars)
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
        s+="Field freq    : "+str(self.freq*ha2ev)+" [eV] \n"
        s+="Field period  : "+str(self.period/fs2aut)+" [fs] \n"
        s+="Time step     : "+str(self.t_step/fs2aut)+" [fs] \n"
        s+="\n =================================\n\n"
        s+="Max Floq mode : "+str(self.max_fl_mode)+"\n"
        s+="Tot Floq modes: "+str(self.tot_fl_modes)+"\n"
        s+="\n =================================\n\n"
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
            qe_thrs,err_thrs,step,max_iter,max_fl_mode = vb_pars
            max_fl_mode = min(max_fl_mode,vb_db.fl_order)
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

    def qe_from_ratio(self,qe_ref=None,component=None):
        """Floquet quasienergy directly from the ratio of the raw coefficients
           at times differing by one period.

           For a Floquet state c_j(t) = exp(-i*qe*t/hbar) * p_j(t) with p_j(t+T)=p_j(t),
           hence  c_j(t+T)/c_j(t) = exp(-i*qe*T/hbar)   for every j and t.
           All (t,j) pairs are combined in the |c|^2-weighted estimator
               z = sum_{t,j} conj(c_j(t)) * c_j(t+T)   ->   qe = -hbar*arg(z)/T
           (so tiny, noisy components do not spoil the result).
           The phase is defined only modulo 2*pi, i.e. qe is defined modulo hbar*omega:
           the branch closest to qe_ref (default: KS energy) is returned.

           component : None (use all basis states) or 0-based index in the basis,
                       to use the ratio of a single coefficient only.
           Returns (qe_ev, info) with info a dict of diagnostics:
               eps_t     : qe estimated at each starting time (should be flat)
               spread    : std of eps_t [eV]
               coherence : |z| / sum|c(t)||c(t+T)|, equal to 1 for a perfect Floquet state
        """
        if qe_ref is None:
            qe_ref = self.ks_ev
        t = np.asarray(self.times)
        c = np.asarray(self.tvecs)
        T = self.period
        dt = np.diff(t)
        if t[-1] - t[0] < T:
            raise ValueError("Simulation shorter than one field period: cannot use the ratio method")

        shift = T/dt[0]
        if np.allclose(dt,dt[0],rtol=1e-6) and abs(shift-round(shift)) < 1e-3:
            # period is an integer number of time steps: no interpolation needed
            s = int(round(shift))
            c0, c1 = c[:len(t)-s], c[s:]
        else:
            print("\n\n WARNING!!! Period is not an integer number of the time steps \n\n")
            print("Time step :",str(dt[0])," fs ")
            print("Period    :",str(T)," fs ")
            # period not commensurate with the time grid: interpolate c(t+T)
            from scipy.interpolate import CubicSpline
            mask = (t + T) <= (t[-1] + 1e-9)
            c0 = c[mask]
            c1 = CubicSpline(t,c,axis=0)(t[mask]+T)

        if component is not None:
            c0, c1 = c0[:,component:component+1], c1[:,component:component+1]

        z_t = np.sum(np.conj(c0)*c1,axis=1)
        z   = np.sum(z_t)
        # phase is defined modulo 2pi -> qe modulo hbar*omega: choose branch closest to reference
        def _branch(eps0):
            return eps0 + np.round((qe_ref-eps0)/self.freq)*self.freq
        qe   = _branch(-hbar_eVfs*np.angle(z)/T)
        eps_t = _branch(-hbar_eVfs*np.angle(z_t)/T)
        coherence = np.abs(z)/np.sum(np.abs(c0)*np.abs(c1))
        info = {'eps_t':eps_t,'spread':np.std(eps_t),'coherence':coherence,
                'n_pairs':len(z_t)}
        return qe, info

    def plot_raw_tvecs(self,band_to_plot=1,outdir=None):
        """Plot the raw coefficients c_j(t) of the propagated state on KS state j
           (before removing the Floquet phase factor exp(-i qe t/hbar)).
           band_to_plot is 1-based position in the basis, as in plot_realtime.
        """
        c = np.array(self.tvecs)[:,band_to_plot-1]
        t = np.array(self.times)
        outdir = outdir or f'figs-raw_k{self.kpt}_b{self.band}'
        os.makedirs(outdir,exist_ok=True)
        fig,axes = plt.subplots(3,sharex=True)
        fig.set_size_inches(8.3,7.5)
        fig.suptitle(f'Raw coefficient on KS state {band_to_plot} (kpt {self.kpt}, band {self.band})')
        axes[0].plot(t,c.real,marker='.',color='tab:blue'); axes[0].set_ylabel(f'Re[ c_{band_to_plot} ]')
        axes[1].plot(t,c.imag,marker='.',color='tab:blue'); axes[1].set_ylabel(f'Im[ c_{band_to_plot} ]')
        axes[2].plot(t,np.abs(c),marker='.',color='tab:blue'); axes[2].set_ylabel(f'|c_{band_to_plot}|')
        axes[2].set_xlabel('Time (fs)')
        plt.show()
        #plt.savefig(f'{outdir}/fig-raw_coefficient_KS_state_{band_to_plot}.pdf')
        plt.close()

    def find_qe(self,qe_ev=None,tag=None,method='optimize',component=None):
        """method = 'optimize' : secant solver over run_nl2fl, minimising the
                                 error between nl_in and nl_out (original behaviour)
                    'ratio'    : quasienergy from c(t+T)/c(t) (see qe_from_ratio),
                                 no iteration; run_nl2fl is called once to get the
                                 Floquet coefficients and the periodicity error
           qe_ev     : initial guess ('optimize') / branch reference ('ratio')
        """
        if qe_ev   is None:
            qe_ev = self.ks_ev

        print("Initial KS energies : \n")
        print(qe_ev)
        sys.exit(0)
        if method == 'ratio':
            return self._find_qe_ratio(qe_ev,tag=tag,component=component)
        if method != 'optimize':
            raise ValueError("method must be 'optimize' or 'ratio'")
        return self._find_qe_optimize(qe_ev,tag=tag)

    def _find_qe_ratio(self,qe_ev,tag=None,component=None):
        qe, info = self.qe_from_ratio(qe_ref=qe_ev,component=component)
        evecs = self.run_nl2fl(qe,tag=tag,iter_num=0)
        evecs.nr_it  = 0
        evecs.nr_acc = info['spread']
        evecs.ratio_info = info
        return evecs

    def _find_qe_optimize(self,qe_ev,tag=None):
        """ Secant solver to iterate over executions of run_nl2fl
            and minimize the error between nl_in and nl_out
        """
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

    def plot_realtime(self,band_to_plot=1,t_step=0.0025):
        n_steps = int(self.period*4/t_step)+10
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
        plt.show()
        #plt.savefig(f'{self.dir}/fig-real_time_projection_over_KS_state_{band_to_plot+1}.pdf')
        #plt.close()

    def plot_realtime_raw(self,band_to_plot=1,t_step=0.0025):
        """Same as plot_realtime but for the coefficient BEFORE removing the Floquet
           phase, c(t) = exp(-i qe t/hbar) * p(t). Data points: raw NL coefficients;
           line: Floquet reconstruction exp(-i qe t/hbar) * sum_eta f_eta exp(-i eta w t).
        """
        n_steps = int(self.period*4/t_step)+10
        X,Y = self.calc_all_times_pVecs(t_step=t_step,n_steps=n_steps)
        phase_X = np.exp(-1j*self.FL_qe*X/hbar_eVfs)
        phase_t = np.exp(-1j*self.FL_qe*self.times/hbar_eVfs)
        C_fl  = Y[:,band_to_plot]*phase_X
        C_in  = self.NL_in[:,band_to_plot]*phase_t     # = raw tvecs
        C_out = self.NL_out[:,band_to_plot]*phase_t

        fig,axes = plt.subplots(3,sharex=True)
        fig.set_size_inches(8.3,7.5)
        fig.suptitle(f'Raw coefficient (with Floquet phase) on Kohn-Sham state: {band_to_plot+1}')
        for ax,f in zip(axes,[np.real,np.imag,np.abs]):
            ax.plot(X,f(C_fl),label='FL_calculated',color='tab:blue')
            ax.plot(self.times,f(C_in) ,label='NL_in (raw data)',marker='o',linestyle='',color='tab:olive',ms=12)
            ax.plot(self.times,f(C_out),label='NL_out',marker='.',linestyle='',color='tab:blue',ms=11)
        axes[0].set_ylabel(f'Re[ c_{band_to_plot+1} ]')
        axes[1].set_ylabel(f'Im[ c_{band_to_plot+1} ]')
        axes[2].set_ylabel(f'|c_{band_to_plot+1}|')
        axes[2].set_xlabel('Time (fs)')
        axes[0].legend()
        plt.show()
 #       os.makedirs(self.dir,exist_ok=True)
#        plt.savefig(f'{self.dir}/fig-real_time_raw_coefficient_KS_state_{band_to_plot+1}.pdf')
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
   freq = float(ds['Field_Freq_1'][0])
   return freq

def get_bands(kpt=None,band=None):
  save_folder='./SAVE'
  yel=YamboElectronsDB.from_db_file(folder=save_folder)
  efermi=yel.setFermiFixed()  # for insulators only
  return yel.eigenvalues_ibz[0,kpt-1,band-1] #-efermi


##################################################################################
### FLOQUET BAND/PLOT FUNCTIONS ###
##################################################################################

def findallk_qe(vbdb,report_file=None,method='optimize'):
    from contextlib import nullcontext
    #
    ks_evk = np.zeros([vbdb.n_kpts,vbdb.n_vbands])
    fl_qek = np.zeros([vbdb.n_kpts,vbdb.n_vbands])
    fl_eig = np.zeros([vbdb.n_kpts,vbdb.n_vbands,vbdb.fl_order*2+1,vbdb.basis_size],dtype=complex) # COMPLEX CHANGE IT

    with open(report_file,'w') if report_file else nullcontext(sys.stdout) as f:
        for _bnd in range(vbdb.n_vbands):
            for _kpt in range(vbdb.n_kpts):
                vbs = VbPP(vbdb,_kpt+1,_bnd+1)
                opt_evecs = vbs.find_qe(method=method)
                ks_evk[_kpt,_bnd] = vbs.ks_ev
                fl_qek[_kpt,_bnd] = opt_evecs.FL_qe
                fl_eig[_kpt,_bnd,:,:] = opt_evecs.FL_vecs[:,:]
                f.write(f'Bnd {_bnd+1:2} @Kpt {_kpt+1:2}: KS eval [eV] = {ks_evk[_kpt,_bnd]:.6f} --- FL qe [eV] = {fl_qek[_kpt,_bnd]:.6f} -- Diff [eV] = {vbs.ks_ev-opt_evecs.FL_qe:.2e} -- It.: {opt_evecs.nr_it:2} -- err in period = {opt_evecs.err:.2e}\n')
    return fl_eig, fl_qek, ks_evk  

def calc_rho(vbdb,fl_eig,order,harmonic,kpt,bnd_1,bnd_2):
    '''
    vbdb     :: databases with info on vbands as output by YanboVbandsDB
    fl_eig   :: Floquet eigenvectors as output by findallk_qe
    order    :: integer, perturbation order of the field, Gamma (positive) 
    harmonic :: integer, harmonic order, gamma (non negative)
    kpt,band1,band2 :: integers, k point, and band indexes of the element of the density matrix

    RETURNS
    rho :: density matrix element (kpt,band1,band2) contributing to harmonic order gamma 
           for perturbation order Gamma, for all possible combinations of harmonic orders eta, nu of the fl_eig 
    '''
    # dimension ckecks
    args =locals()
    val =list(args.values())
    key = list(args.keys())
    for i in range(len(val)):
        if not(any([str(key[i]) == 'vbdb', str(key[i]) == 'fl_eig'])): 
            if not(isinstance(val[i], int)):
                raise ValueError("value for "+str(key[i])+" must be integer")        
            if val[i] <0:
                raise ValueError("value for "+str(key[i])+" must be nonnegative")
    if not((order - harmonic)%2 == 0 and ((order - harmonic)>=0)):
        raise ValueError("Harmonic order gamma not compatible with perturbation order Gamma: Gamma = gamma + 2n, where n is 0,1,...")
    if kpt >= vbdb.n_kpts:
        raise ValueError("Value for kpt must be smaller than "+str(vbdb.n_kpts))
    if any([bnd_1 < vbdb.basis_index[0], bnd_1 > vbdb.basis_index[1], bnd_2 < vbdb.basis_index[0], bnd_2 > vbdb.basis_index[1] ]):
        raise ValueError("Basis indexes must be in range "+ str( vbdb.basis_index))
    # basis index 
    b_1 = bnd_1 -  vbdb.basis_index[0]
    b_2 = bnd_2 -  vbdb.basis_index[0]
    # determine possible eta,nu values 
    if harmonic == order:
        eta = list(range(-harmonic,1))
    elif harmonic < order:
        eta = [-int((harmonic+order)/2),int((order-harmonic)/2)]
    nu = list(int(x + harmonic) for x in eta)
    c_dim = len(eta)
    rho = np.zeros(c_dim,dtype=complex)
    path = []
    for ii in range(c_dim):
         for n in range(vbdb.n_vbands):
             rho[ii] = rho[ii] + np.conj(fl_eig[kpt,n,eta[ii],b_2])*fl_eig[kpt,n,nu[ii],b_1]
         path.append([eta[ii],nu[ii]])
    return rho, path

def get_qebands_path(lat,path,fl_qek,ks_evk): # see how to put kwargs into it...
    from yambopy.plot.bandstructure import YambopyBandStructure
    from yambopy.lattice import red_car
    from yambopy.kpoints import get_path

    bands_kpoints, bands_indexes, path_car = get_path(lat.car_kpoints,lat.rlat,lat.sym_car,path,debug=False) 

    ks_bs = YambopyBandStructure(ks_evk[bands_indexes],bands_kpoints,kpath=path_car) 
    fl_bs = YambopyBandStructure(fl_qek[bands_indexes],bands_kpoints,kpath=path_car)

    return ks_bs, fl_bs
    
def get_qebands_interpolate(lat,path,fl_qek,ks_evk,lpratio=5,fermie=0,nelect = 0,verbose=1): 
    
    cell = (lat.lat, lat.red_atomic_positions, lat.atomic_numbers)
    trev_for_interp = lat.time_rev
    symrel = [sym for sym,trev in zip(lat.sym_rec_red,lat.time_rev_list) if trev==False ]
    kpoints = lat.red_kpoints
    _, _, path_car = get_path(lat.car_kpoints,lat.rlat,lat.sym_car,path,debug=False)
    band_kpoints_rlu = path.get_klist()[:,:3]
    band_kpoints = red_car(band_kpoints_rlu,lat.rlat) # why different from above???

    fl_eigens = np.zeros([2,np.shape(fl_qek)[0], np.shape(fl_qek)[1]])
    ks_eigens = np.zeros([2,np.shape(ks_evk)[0], np.shape(ks_evk)[1]])
    for ik in range(np.shape(fl_qek)[0]):
        fl_eigens[0,ik,:] = sorted(fl_qek[ik,:])
        ks_eigens[0,ik,:] = sorted(ks_evk[ik,:])
    skw_fl = SkwInterpolator(lpratio,kpoints,fl_eigens,fermie,nelect,cell,symrel,trev_for_interp,verbose=verbose)
    skw_ks = SkwInterpolator(lpratio,kpoints,ks_eigens,fermie,nelect,cell,symrel,trev_for_interp,verbose=verbose)

    fl_eigens_kpath = skw_fl.interp_kpts(band_kpoints_rlu).eigens[0]
    fl_bs = YambopyBandStructure(fl_eigens_kpath,band_kpoints,kpath=path_car)
    ks_eigens_kpath = skw_ks.interp_kpts(band_kpoints_rlu).eigens[0]
    ks_bs = YambopyBandStructure(ks_eigens_kpath,band_kpoints,kpath=path_car)

    return ks_bs, fl_bs

@add_fig_kwargs
def plot_2D_kdist(data,lat,nspin=-1,plt_cbar=False,shift_BZ=True,**kwargs):
    """
    2D scatterplot in the k-BZ of any real quantity which is a function of only the k-grid
    
    data:: real quantity function of k-grid
        - if plt_cbar colorbar is shown
        - if shift_BZ adjacent BZs are also plotted (default)
        - kwargs example: marker='H', s=300, cmap='viridis', etc.
        
    NB: THIS IS A COPY OF PLOT_DIPOLE except it is out of the dipole class and some minor detail
    """
    from yambopy.plot.plotting import add_fig_kwargs,BZ_Wigner_Seitz

    kpts =lat.car_kpoints
    rlat = lat.rlat

    # Input check
    if len(data)!=len(kpts):
        raise ValueError('Something wrong in data dimensions (%d data vs %d kpts)'%(len(data),len(kpts)))

    # Global plot stuff
    fig, ax = plt.subplots(1, 1)
    ax.add_patch(BZ_Wigner_Seitz(lat))

    if plt_cbar:
        if 'cmap' in kwargs.keys(): color_map = plt.get_cmap(kwargs['cmap'])
        else:                       color_map = plt.get_cmap('viridis')
    lim = 1.05*np.linalg.norm(rlat[0])
    ax.set_xlim(-lim,lim)
    ax.set_ylim(-lim,lim)

    # Reproduce plot also in adjacent BZs
    if shift_BZ:
        BZs = shifted_grids_2D(kpts,rlat)
        for kpts_s in BZs: plot=ax.scatter(kpts_s[:,0],kpts_s[:,1],c=data,**kwargs)
    else:
        plot=ax.scatter(kpts[:,0],kpts[:,1],c=data,**kwargs)

    if plt_cbar: cbar = fig.colorbar(plot)

    plt.gca().set_aspect('equal')

    return ax

    #if plt_show: plt.show()
    #else: print_string = "Plot ready.\nYou can customise adding savefig, title, labels, text, show, etc..."
