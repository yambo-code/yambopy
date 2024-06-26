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
from mpl_toolkits.mplot3d import Axes3D
from yambopy.units import ha2ev, ev2cm1, I
from yambopy.plot.plotting import add_fig_kwargs,BZ_Wigner_Seitz,shifted_grids_2D

class YamboElectronPhononDB():
    """
    Python class to read the electron-phonon matrix elements from yambo.
    
    By default it reads the full databases including fragments.
    
    - Input: YamboLatticeDB object, paths of ndb.elph* and ns.db1
    - Input: if not read_all, read only header
    
    - Usage and main variables: 
    
      :: yph = YamboElectronPhononDB(ylat,folder_gkkp=path1,save=path2)
    
      :: yph.ph_energies     #Phonon energies (eV)      
      :: yph.ph_eigenvectors #Phonon modes
      :: yph.gkkp            #El-ph matrix elements (by default normalised with ph. energies):
      :: yph.gkkp_sq         #Couplings (square) 

      Additional variables (Experimental stuff)
      :: yph.gkkp_bare
      :: yph.gkkp_bare_sq
      :: yph.gkkp_mixed      #Coupling (mixed bare-dressed)
   
    Formats:
    - modes[il][iat][ix]
    - gkkp[iq][ik][il][ib1][ib2]

    Plots provided:
    - Call function plot for scatterplot in the k-BZ of any quantity A(k)_{iq,ib1,ib2,inu} 
        -- if plt_show, show plot at runtime
        -- if plt_cbar, add colorbar

      Example, plot of |g(k)_{0,3,4,4}|:      
           :: yph.plot_elph( np.abs(yph.gkkp[0,:,3,4,4]) )              
    """
    def __init__(self,lattice,filename='ndb.elph_gkkp',folder_gkkp='SAVE',save='SAVE',read_all=True):
        
        self.lattice = lattice

        # Find correct database names
        if os.path.isfile("%s/ndb.elph_gkkp"%folder_gkkp): filename='%s/ndb.elph_gkkp'%folder_gkkp
        elif os.path.isfile("%s/ndb.elph_gkkp_expanded"%folder_gkkp): filename='%s/ndb.elph_gkkp_expanded'%folder_gkkp
        else: filename = "%s/%s"%(folder_gkkp,filename)
        self.frag_filename = filename + "_fragment_"
        self.are_bare_there = False
        
        # necessary lattice information
        self.alat        = lattice.alat
        self.rlat        = lattice.rlat
        self.car_kpoints = lattice.car_kpoints
            
        # Check if databases exist. Exit only if header is missing.
        try: database = Dataset(filename)
        except: raise FileNotFoundError("error opening %s in YamboElectronPhononDB"%filename)
        
        try: database_frag = Dataset("%s1"%self.frag_filename)
        except FileNotFoundError: print("[WARNING] Database fragment at q=0 not detected")
        else: 
            # Check if bare matrix elements are present
            try: 
                database_frag.variables['ELPH_GKKP_BARE_Q1']
            except KeyError: 
                database_frag.close()
            else:
                self.are_bare_there = True
                database_frag.close()
                 
        #read qpoints    
        self.qpoints = database.variables['PH_Q'][:].T
        self.car_qpoints = np.array([ q/self.alat for q in self.qpoints ])
        #read dimensions of electron phonon parameters
        self.nmodes, self.nqpoints, self.nkpoints, b_1, b_2 = database.variables['PARS'][:5].astype(int)
        if b_1>b_2: # Old database (no GkkpBands in PARS)
            self.nbands = b_1
            self.b_in, self.b_out = [0,self.nbands-1]
        else: # New database (PARS with GkkpBands)
            self.b_in, self.b_out = [b_1-1,b_2-1]
            self.nbands = b_2-b_1+1
        self.natoms = int(self.nmodes/3)
        # read IBZ k-points
        self.ibz_kpoints_elph = database.variables['HEAD_KPT'][:].T
        self.ibz_car_kpoints = np.array([ k/self.alat for k in self.ibz_kpoints_elph ])
        try: # Check if full K-point list is provided (upon expansion), otherwise use the one from ns.db1
            self.kpoints_elph = database.variables['PH_K'][:].T
            self.car_kpoints = np.array([ k/self.alat for k in self.kpoints_elph ])
            database.close()
        except KeyError:
            database.close()
        
        #Check how many databases are present
        self.nfrags = self.nqpoints
        for iq in range(self.nqpoints):
            if not os.path.isfile("%s%d"%(self.frag_filename,iq+1)): 
                self.nfrags = iq
                break
        
        # Keep reading
        if read_all: self.read_full_DB()
        
    def read_full_DB(self):
        """
        Read all variables in the ndb.elph_gkkp* dbs as attributes of this class
        """
        
        # Frequencies
        self.read_frequencies()
        
        # Eigenmodes
        self.read_eigenmodes()

        # <dVscf> matrix elements plus <dVbare> if they exist
        self.read_elph()
        if self.are_bare_there: self.read_elph(kind='bare')   
   
        # Get square matrix elements
        self.get_gkkp_sq()
        
        # Get the symmetrised dressed-bare coupling
        if self.are_bare_there: self.get_gkkp_mixed()

    def read_frequencies(self):
        """
        Read phonon frequencies in eV
        """
        self.ph_energies  = np.zeros([self.nfrags,self.nmodes])
        
        for iq in range(self.nfrags):
            fil = self.frag_filename + "%d"%(iq+1)
            database = Dataset(fil)
            self.ph_energies[iq] = np.sqrt(database.variables['PH_FREQS%d'%(iq+1)][:])*ha2ev
            database.close()
        
    def read_eigenmodes(self):
        """
        Read phonon eigenmodes
        """
        self.ph_eigenvectors = np.zeros([self.nfrags,self.nmodes,self.natoms,3],dtype=np.complex64)
        
        for iq in range(self.nfrags):
            fil = self.frag_filename + "%d"%(iq+1)
            database = Dataset(fil)
            #eigs_q[cartesian][atom][mode][complex]
            eigs_q = database.variables['POLARIZATION_VECTORS'][:].T
            self.ph_eigenvectors[iq] = eigs_q[0,:,:,:] + eigs_q[1,:,:,:]*I
            database.close()
             
    def read_elph(self,kind='dressed',scale_g_with_ph_energies=True):
        """
        Read electron-phonon matrix elements
        
        - kind is 'dressed' or 'bare'
        - var_nm is 'ELPH_GKKP_Q' or 'ELPH_GKKP_BARE_Q'
        - If scale_g_with_ph_energies they are divided by sqrt(2*ph_E)

        NB: ELPH_GKKP_Q is saved by yambo as (2,mode,bnd1,bnd2,k), but netCDF stores
            the *transpose* (k,bnd2,bnd1,mode,2).
            We want to change it to complex (k,mode,bnd1,bnd2), therefore we need to
            *swap* bnd2<->mode.
        """    
        if kind!='dressed' and kind!='bare': 
            raise ValueError("Wrong kind %s (can be 'dressed' [Default] or 'bare')"%kind) 
        if kind=='dressed': var_nm = 'ELPH_GKKP_Q'
        if kind=='bare':    var_nm = 'ELPH_GKKP_BARE_Q'       
        
        # gkkp[q][k][mode][bnd1][bnd2]
        gkkp_full = np.zeros([self.nfrags,self.nkpoints,self.nmodes,self.nbands,self.nbands],dtype=np.complex64)   
        
        for iq in range(self.nfrags):
            fil = self.frag_filename + "%d"%(iq+1)
            database = Dataset(fil)
            gkkp = database.variables['%s%d'%(var_nm,iq+1)][:]
            gkkp_full[iq] = np.swapaxes(gkkp[:,:,:,:,0] + I*gkkp[:,:,:,:,1],-1,1)
            database.close()
        
        # Check integrity of elph values
        if np.isnan(gkkp_full).any(): print('[WARNING] NaN values detected in elph database.')
        
        # Scaling with phonon energies
        if scale_g_with_ph_energies: gkkp_full = self.scale_g(gkkp_full)    
                
        if kind=='dressed': self.gkkp = gkkp_full
        if kind=='bare': self.gkkp_bare = gkkp_full          
    
    def scale_g(self,dvscf):
        """
        Normalise matrix elements by the phonon energy as: 
       
        g_qnu = dvscf_qnu/sqrt(2*w_qnu)
        """
        
        g = np.zeros([self.nfrags,self.nkpoints,self.nmodes,self.nbands,self.nbands],dtype=np.complex64)
        for iq in range(self.nfrags):
            for inu in range(self.nmodes): 
                if iq==0 and inu in [0,1,2]: 
                    g[iq,:,inu,:,:] = 0. # Remove acoustic branches
                else:
                    ph_E = self.ph_energies[iq,inu]/ha2ev # Put back the energies in Hartree units
                    g[iq,:,inu,:,:] = dvscf[iq,:,inu,:,:]/np.sqrt(2.*ph_E)
        return g

    def get_gkkp_sq(self,read_bare=False):
        """
        Return g^2
        """
        self.gkkp_sq = np.abs(self.gkkp)**2. 
        if self.are_bare_there: self.gkkp_bare_sq = np.abs(self.gkkp_bare)**2. 

    def get_gkkp_mixed(self):
        """
        Return the symmetrised dressed-bare coupling
        """
        self.gkkp_mixed = np.real(self.gkkp)*np.real(self.gkkp_bare)+np.imag(self.gkkp)*np.imag(self.gkkp_bare)
        
    @add_fig_kwargs
    def plot_elph(self,data,kcoords=None,plt_show=False,plt_cbar=False,**kwargs):
        """
        2D scatterplot in the BZ:

         (i)  in k-space of the quantity A_{k}(iq,inu,ib1,ib2).
         (ii) in q-space of the quantity A_{q}(ik,inu,ib1,ib2).       

        Any real quantity which is a function of only the k-grid or q-grid may be supplied.
        The indices iq/ik,inu,ib1,ib2 are user-specified.

        - kcoords refers to the k/q-grid in Cartesian coordinates (i.e., yelph.car_qpoints and similar).
          If None is specified, a k-space, fixed-q plot is assumed.
        
        - if plt_show plot is shown
        - if plt_cbar colorbar is shown
        - kwargs example: marker='H', s=300, cmap='viridis', etc.
        
        NB: So far requires a 2D system. 
            Can be improved to plot BZ planes at constant k_z for 3D systems.
        """        
        if kcoords is None: kpts = self.car_kpoints # Assume k-space plot
        else:               kpts = kcoords # Plot on momentum map supplied by user       

        # Input check
        if len(data)!=len(kpts): 
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
        BZs = shifted_grids_2D(kpts,self.rlat)
        for kpts_s in BZs: plot=self.ax.scatter(kpts_s[:,0],kpts_s[:,1],c=data,**kwargs)
        
        if plt_cbar: self.fig.colorbar(plot)
        
        plt.gca().set_aspect('equal')

        if plt_show: plt.show()
        else: print("Plot ready.\nYou can customise adding savefig, title, labels, text, show, etc...")
        
    def __str__(self,verbose=False):

        try: self.ph_energies
        except AttributeError: self.read_frequencies()

        try: self.ph_eigenvectors
        except AttributeError: self.read_eigenmodes()

        lines = []; app = lines.append
        app(marquee(self.__class__.__name__))
            
        app('nqpoints: %d'%self.nqpoints)
        app('nkpoints: %d'%self.nkpoints)
        app('nmodes: %d'%self.nmodes)
        app('natoms: %d'%self.natoms)
        app('nbands: %d (%d - %d)'%(self.nbands,self.b_in,self.b_out))
        if self.nfrags == self.nqpoints: app('fragments: %d'%self.nfrags)
        else: app('fragments: %d [WARNING] nfrags < nqpoints'%self.nfrags)
        if self.are_bare_there: app('bare couplings are present')
        if verbose:
            app('-----------------------------------')
            for iq in range(self.nfrags):
                app('nqpoint %d'%iq)
                for n,mode in enumerate(self.ph_eigenvectors[iq]):
                    app('mode %d freq: %lf meV'%(n,self.ph_energies[iq,n]*1000.))
                    for a in range(self.natoms):
                        app(("%12.8lf "*3)%tuple(mode[a].real))
            app('-----------------------------------')
        return "\n".join(lines)
