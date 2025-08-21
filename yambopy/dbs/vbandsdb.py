# Copyright (c) 2025, Ignacio M Alliati & Myrta Gruening
# All rights reserved.
#
# This file is part of the yambopy project
#
from yambopy import *
from yambopy.plot import *
import numpy as np
import sys
import os
#
# This class reads data from the ndb.V_bands and the ndb.V_bands_K_section databases 
# 
# Time series is stored in IO_TIME_points variable
# 
class YamboVbandsDB():
    """
    Open the Vbands databases and store it in a VBandsDB class.
    """
    def get_n_timesteps(self):
        """gets the number of steps
        """
        ds=Dataset(self.vb_path+'/ndb.RT_V_bands')
        number_of_steps = int(ds['IO_TIME_steps_last_nsteps'][2])
        return number_of_steps
    
    def get_basis_size(self):
        """gets the size of the basis
        """
        ds=Dataset(self.vb_path+'/ndb.RT_V_bands_K_section')
        basis_size=int(ds.dimensions['RT_nbands'].size)
        n_vbands=int(ds.dimensions['nb_full'].size)
        n_kpts=int(ds.dimensions['QP_nk'].size)
        sp_pol = int(ds.dimensions['n_sp_pol'].size) # note all the rest is for spin unpolarised...
        if sp_pol != 1:
            raise Exception('Spin polarised case not yet implemented.')
        return basis_size, n_vbands, n_kpts#, sp_pol

    def get_basis_idx(self):
        """gets the basis index
        """
        basis_idx = [-1,-1]
        ds=Dataset(self.jobdir+'/ndb.RT_V_bands')
        basis_idx[0] = int(ds['RT_bands_kpts'][0])   
        basis_idx[1] = int(ds['RT_bands_kpts'][1])   
        return basis_idx
    
    def get_tvecs(self):
        """ 
        """
        ds=Dataset(self.vb_path+'/ndb.RT_V_bands_K_section')

        list_of_evecs = []
        for it in range(1,self.n_timesteps):
            evec=np.zeros((self.n_kpts,self.n_vbands,self.basis_size),dtype=complex)

            for kpt in range(self.n_kpts):
                for band in range(self.n_vbands):
                    for i in range(self.basis_size):
                        evec[kpt,band,i]=complex(
                        ds['V_bands'][it,0,kpt,band,i,0],ds['V_bands'][it,0,kpt,band,i,1])

            list_of_evecs.append(evec)

        return list_of_evecs

    def get_times(self):
        """gets times in fs as list
        """
        list_of_times = []
        ds=Dataset(self.vb_path+'/ndb.RT_V_bands')
        for i in range(1,self.n_timesteps):
            time_fs=float(ds['IO_TIME_points'][i]) #should I do it here?
            list_of_times.append(time_fs)
        return list_of_times

    def get_frequency(self):
        ds=Dataset(self.jobdir+'/ndb.Nonlinear')
        freq = float(ds['Field_Freq_range_1'][0])
        period = 2*np.pi / freq
        return freq,period

    def get_florder(self):
        """get the Floquet order
        """
        ds=Dataset(self.vb_path+'/ndb.RT_V_bands')
        floquet_order = int(ds['IO_Floquet_order'][0])
        return floquet_order

    def __init__(self,folder=".",calc="SAVE"):
        """folder: type string
           calc: type string
        """
        self.vb_path = '%s/%s/%s'%(folder,calc,nl_db)
        for ndb in ['/ndb.RT_V_bands','/ndb.RT_V_bands_K_section']:
            try:
                data_obs = Dataset(self.vb_path+ndb)
            except:
                raise ValueError("Error reading V_bands databases at %s"%self.vb_path)

        self.n_timesteps = self.get_ntimesteps()
        self.basis_size = self.get_basis_size()
        self.basis_index = self.get_basis_idx()
        self.tvecs = self.get_tvecs()
        self.times = self.get_times()
        self.freq, self.period = self.get_frequency()
        self.fl_order = self.get_florder()

    def __str__(self):
        """
        Print all info of the class
        """
        s="\n * * * ndb.V_bands dbs data * * * \n\n"
        s+="N timesteps   : "+str(self.n_timesteps)+"\n"
        s+="Basis size    : "+str(self.basis_size)+"\n"
        s+="Basis index   : "+str(self.basis_index)+"\n"
        s+="Number Kpt    : "+str(self.n_kpts)+"\n"
        s+="Number Vbands : "+str(self.n_vbands)+"\n"
        s+="Floquet order : "+str(self.fl_order)+"\n"
        # for i,v in enumerate(self.tvecs):
        #     s+="Time: "+str(self.times[i])+" au:\n"
        #     s+='idx  bnd   c.real             c.imag\n'
        #     for j in range(self.basis_size):
        #         idx = self.basis_index[0] + j
        #         s+=str(j)+"  "+str(idx)+"   "+str(v[j].real)+"   "+str(v[j].imag)+"\n"
        return s
