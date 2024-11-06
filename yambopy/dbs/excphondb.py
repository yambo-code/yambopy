# Copyright (c) 2018, Fulvio Paleari, Alejandro Molina-Sánchez, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
from yambopy import *
from netCDF4 import Dataset
from math import sqrt
import numpy as np
from yambopy.tools.string import marquee
import os
import matplotlib.pyplot as plt
from yambopy.units import ha2ev, I
from yambopy.plot.plotting import add_fig_kwargs,BZ_hexagon,shifted_grids_2D

class YamboExcitonPhononDB():
    """
    Python class to read the exciton-phonon matrix elements from yambo.

    By default it reads the full databases including fragments.

    - Input: YamboLatticeDB object
    - Input: paths of ndb.excph*
    - Input: if not read_all, read only header

    NB: So far we read G(exc_in_Q=0, exc_out_q, ph_q) 

    - Usage and main variables:

      :: yexcph = YamboExcitonPhononDB(ylat,save_excph=path1)

      :: yexcph.car_qpoints  #Momenta cc, iq
      :: yexcph.red_qpoints  #Momenta rlu, iq
      :: yexcph.exciton_in   #Exciton states in, iexc1
      :: yexcph.exciton_out  #Exciton states out (summed in self-energy), iexc2
      :: yexcph.phonon_modes #Phonon modes, il
      :: yexcph.excph        #Exc-ph matrix elements
      :: yexcph.excph_sq     #Exc-ph matrix elements squared

    Formats:
    - excph[iq][il][iexc1][iexc2]

    Plots provided:
    - Call function plot for scatterplot in the q-BZ of any quantity A(q)_{inu,ib1,ib2}
        -- if plt_show, show plot at runtime
        -- if plt_cbar, add colorbar

      Example, plot of |G(q)_{3,4,4}|:
           :: yexcph.plot_excph( np.abs(yexcph.excph[:,3,4,4]) )
    """
    def __init__(self,lattice,save_excph="./",read_all=True):
        
        # Find databases
        if os.path.isfile("%s/ndb.excph_gkkp"%save_excph): filename='%s/ndb.excph_gkkp'%save_excph
        else: raise FileNotFoundError("Databases not found.")
        self.frag_filename = filename + "_fragment_"
                    
        # Check if databases exist. Exit only if header is missing.
        try: database = Dataset(filename)
        except: raise FileNotFoundError("error opening %s in YamboElectronPhononDB"%filename)
        
        try: database_frag = Dataset("%s1"%self.frag_filename)
        except FileNotFoundError: print("[WARNING] Database fragment at q=0 not detected")
        else: database_frag.close()
                 
        #read dimensions of exciton phonon parameters
        self.nexc_i = database.variables['EXCITON_STATES'][1].astype(int)
        self.nexc_o = database.variables['EXCITON_SUM'][1].astype(int)
        self.nmodes = database.variables['PHONON_MODES'][0].astype(int)
        self.nqpoints = database.variables['HEAD_R_LATT'][3].astype(int)
        self.type_exc_i = database.variables['L_kind_in'][...].tostring().decode().strip()
        self.type_exc_o = database.variables['L_kind_out'][...].tostring().decode().strip()
        database.close()

        #Check how many databases are present
        self.nfrags = self.nqpoints
        for iq in range(self.nqpoints):
            if not os.path.isfile("%s%d"%(self.frag_filename,iq+1)): 
                self.nfrags = iq
                break
        
        # Necessary lattice information
        self.lattice = lattice
        self.alat    = lattice.alat
        self.rlat    = lattice.rlat
        
        # Keep reading
        if read_all: self.read_full_DB()

    def read_full_DB(self):
        """
        Read more variables in the ndb.excph_gkkp_fragments* dbs as attributes of this class
        """
        
        # qpoints
        self.read_qpoints()
        
        # Matrix elements
        self.read_excph()
        
    def read_qpoints(self):
        """
        Read q points and return cartesian and reduced coordinates
        """
        var_nm = "EXCPH_Q"
        
        self.qpoints = np.zeros([self.nfrags,3])

        for iq in range(self.nfrags):
            fil = self.frag_filename + "%d"%(iq+1)
            database = Dataset(fil)
            self.qpoints[iq] = database.variables['%s%d'%(var_nm,iq+1)][:].T        
    
        self.car_qpoints = np.array([ q/self.alat for q in self.qpoints ])

    def read_excph(self):
        """
        Read exciton-phonon matrix elements and their modulus squared
        
        NB: EXCPH_GKKP_Q is saved by yambo as (2,mode,exc_out,exc_in), but netCDF stores
            the *transpose* (exc_in,exc_out,mode,2).
            We want to change it to complex (iq,mode,exc_in,exc_out)
        """    
        var_nm    = "EXCITON_PH_GKKP_Q"
        var_sq_nm = "EXCITON_PH_GKKP_SQUARED_Q"
            
        # excph[q][mode][iexc1][iexc2]
        excph_full    = np.zeros([self.nfrags,self.nmodes,self.nexc_i,self.nexc_o],dtype=np.complex64)  
        excph_sq_full = np.zeros([self.nfrags,self.nmodes,self.nexc_i,self.nexc_o])  
        
        for iq in range(self.nfrags):
            fil = self.frag_filename + "%d"%(iq+1)
            database = Dataset(fil)
            excph = database.variables['%s%d'%(var_nm,iq+1)][:]
            excph_full[iq] = np.moveaxis( excph[:,:,:,0]+I*excph[:,:,:,1], -1,0 )
            #excph_full[iq] = np.swapaxes( np.swapaxes(excph[:,:,:,0] + I*excph[:,:,:,1],-1,0), -1,-2)
            
            excph_sq = database.variables['%s%d'%(var_sq_nm,iq+1)][:]
            excph_sq_full[iq] = np.moveaxis( excph_sq, -1,0)
            #excph_sq_full[iq] = np.swapaxes( np.swapaxes(excph_sq[:,:,:],-1,0), -1,-2)
            database.close()
        
        # Check integrity of elph values
        if np.isnan(excph_full).any(): print('[WARNING] NaN values detected in elph database.')

        self.excph    = excph_full
        self.excph_sq = excph_sq_full

    @add_fig_kwargs
    def plot_excph(self,data,plt_show=False,plt_cbar=False,**kwargs):
        """
        2D scatterplot in the q-BZ of the quantity A_{iq}(ib2,ib1,inu).
        
        Any real quantity which is a function of only the q-grid may be supplied.
        The indices iq,inu,ib1,ib2 are user-specified.
        
        - if plt_show plot is shown
        - if plt_cbar colorbar is shown
        - kwargs example: marker='H', s=300, cmap='viridis', etc.
        
        NB: So far requires a 2D system. 
            Can be improved to plot BZ planes at constant k_z for 3D systems.
        """
              
        qpts = self.car_qpoints
        
        # Input check
        if len(data)!=len(qpts): 
            raise ValueError('Something wrong in data dimensions (%d data vs %d kpts)'%(len(data),len(kpts)))
        
        # Global plot stuff
        self.fig, self.ax = plt.subplots(1, 1)
        self.ax.add_patch(BZ_Wigner_Seitz(self.lattice))
        
        if plt_cbar:
            if 'cmap' in kwargs.keys(): color_map = plt.get_cmap(kwargs['cmap'])
            else:                       color_map = plt.get_cmap('viridis') 
        lim = 1.05*np.linalg.norm(self.rlat[0])
        self.ax.set_xlim(-lim,lim)
        self.ax.set_ylim(-lim,lim)

        # Reproduce plot also in adjacent BZs
        BZs = shifted_grids_2D(qpts,self.rlat)
        for qpts_s in BZs: plot=self.ax.scatter(qpts_s[:,0],qpts_s[:,1],c=data,**kwargs)
        
        if plt_cbar: self.cbar = self.fig.colorbar(plot)
        
        plt.gca().set_aspect('equal')

        if plt_show: plt.show()
        else: print("Plot ready.\nYou can customise adding savefig, title, labels, text, show, etc...")

    def __str__(self,verbose=False):

        lines = []; app = lines.append
        app(marquee(self.__class__.__name__))
            
        app('nqpoints: %d'%self.nqpoints)
        app('nmodes: %d'%self.nmodes)
        app('nexcitons in (L_%s): %d'%(self.type_exc_i,self.nexc_i))
        app('nexcitons out (L_%s): %d'%(self.type_exc_o,self.nexc_o))
        if self.nfrags == self.nqpoints: app('fragments: %d'%self.nfrags)
        else: app('fragments: %d [WARNING] nfrags < nqpoints'%self.nfrags)
        
        return "\n".join(lines)
