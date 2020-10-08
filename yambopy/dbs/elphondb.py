# Copyright (c) 2018, Fulvio Paleari, Alejandro Molina-Sánchez, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
from yambopy import *
from netCDF4 import Dataset
from math import sqrt
import numpy as np
import os
from yambopy.units import ha2ev, ev2cm1, I

class YamboElectronPhononDB():
    """
    Python class to read the electron-phonon matrix elements from yambo.
    
    By default it reads dimension of elph parameters and qpoints.
    
    - Call function read_frequencies for frequencies.
    - Call function read_eigenmodes for phonon modes.
    - Call function read_elph for electron-phonon matrix elements.
    - Call function read_DB to read everything.
    
    Plot(s) provided:
    - Scatterplot in the BZ of G_{nk} = 1/N_q * \sum_{q,nu} | elph_{qnu,knn} |^2         
    """
    def __init__(self,lattice,filename='ndb.elph_gkkp',folder_gkkp='SAVE',save='SAVE'):
        
        filename = "%s/%s"%(folder_gkkp,filename)
        self.frag_filename = filename + "_fragment_"
        
        # necessary lattice information
        self.alat        = lattice.alat
        self.rlat        = lattice.rlat
        self.car_kpoints = lattice.car_kpoints
        self.red_kpoints = lattice.red_kpoints
            
        # Check if databases exist. Exit only if header is missing.
        try: database = Dataset(filename)
        except: raise FileNotFoundError("error opening %s in YamboElectronPhononDB"%self.filename)
        
        try: database_frag = Dataset("%s1"%self.frag_filename)
        except FileNotFoundError: print("[WARNING] Database fragment at q=0 not detected")
        else: database_frag.close() 
        
        #read qpoints    
        self.qpoints = database.variables['PH_Q'][:].T
        self.car_qpoints = np.array([ q/self.alat for q in self.qpoints ])
        #read dimensions of electron phonon parameters
        self.nmodes, self.nqpoints, self.nkpoints, self.nbands = database.variables['PARS'][:4].astype(int)
        self.natoms = int(self.nmodes/3)
        
        database.close()
        
        #Check how many databases are present
        self.nfrags = self.nqpoints
        for iq in range(self.nqpoints):
            if not os.path.isfile("%s%d"%(self.frag_filename,iq+1)): 
                self.nfrags = iq
                break

    def read_frequencies(self):
        """
        Read phonon frequencies
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
    
    def read_elph(self,iq=-1,ik=-1,ib1=-1,ib2=-1,inu=-1):
        """
        Driver to read electron-phonon matrix elements:
        
        - If no options are specified, read all calling read_elph_full
        - If options (q) or (k,b1,b2) are specified, read the appropriate slice of the gkkp
        
        """
        
        # Read a single q
        if iq>-1: 
            if iq<self.nfrags: self.read_elph_q(iq)
            else: raise ValueError("The database corresponding to iq = %d is not present."%iq)
        
        # Read all phonon quantum numbers (q and modes) for specified electronic quantum numbers (k band1 band2)
        elif ( ib1>-1 and ib2>-1 and ik>-1 ):
            if ib1<self.nbands and ib2<self.nbands and ik<self.nkpoints: return self.read_elph_knm(ik,ib1,ib2)
            else: raise ValueError("The database corresponding to (k,n,m) = (%d, %d, %d) is not present."%(ik,ib1,ib2))
            
        # Read a single transition (band1,band2) for all k,q,modes
        elif ( ib1>-1 and ib2>-1 ):
            if ib1<self.nbands and ib2<self.nbands: return self.read_elph_nm(ib1,ib2)
            else: raise ValueError("The database corresponding to transition (n,m) = (%d, %d) is not present."%(ib1,ib2))            
        
                        
        else: self.read_elph_full()

    def read_elph_nm(self,i_n,i_m):
        """
        Read electron-phonon matrix element between fixed electronic states (n,m)
        """
        gkkp_nm = np.zeros([self.nfrags,self.nkpoints,self.nmodes],dtype=np.complex64)
        
        for iq in range(self.nfrags):
            fil = self.frag_filename + "%d"%(iq+1)
            database = Dataset(fil)
            gkkp_nm[iq] = database.variables['ELPH_GKKP_Q%d'%(iq+1)][:,i_n,i_m,:,0] +\
                       I* database.variables['ELPH_GKKP_Q%d'%(iq+1)][:,i_n,i_m,:,1]
            database.close()
        return gkkp_nm
    
    def read_elph_knm(self,ik,i_n,i_m):
        """
        Read electron-phonon matrix element with fixed electronic (k,n,m) and running phononic (q,nu) quantum numbers
        """
        gkkp_knm = np.zeros([self.nfrags,self.nmodes],dtype=np.complex64)
        
        for iq in range(self.nfrags):
            fil = self.frag_filename + "%d"%(iq+1)
            database = Dataset(fil)
            gkkp_knm[iq] = database.variables['ELPH_GKKP_Q%d'%(iq+1)][ik,i_n,i_m,:,0] +\
                        I* database.variables['ELPH_GKKP_Q%d'%(iq+1)][ik,i_n,i_m,:,1]
            database.close()
        return gkkp_knm
        
    def read_elph_q(self,iq):
        """
        Read electron-phonon matrix element at a user-specified q point
        """
        fil = self.frag_filename + "%d"%(iq+1)
        database = Dataset(fil)
        g_q = database.variables['ELPH_GKKP_Q%d'%(iq+1)][:]
        gkkp_q = (g_q[:,:,:,:,0] + I*g_q[:,:,:,:,1]).reshape([self.nkpoints,self.nmodes,self.nbands,self.nbands])
        database.close()
        
        return gkkp_q
        
    def read_elph_full(self):
        """
        Read electron-phonon matrix elements
        """
        # gkkp[q][k][mode][bnd1][bnd2]
        self.gkkp = np.zeros([self.nfrags,self.nkpoints,self.nmodes,self.nbands,self.nbands],dtype=np.complex64)
        
        for iq in range(self.nfrags):
            fil = self.frag_filename + "%d"%(iq+1)
            database = Dataset(fil)
            gkkp = database.variables['ELPH_GKKP_Q%d'%(iq+1)][:]
            self.gkkp[iq] = (gkkp[:,:,:,:,0] + I*gkkp[:,:,:,:,1]).reshape([self.nkpoints,self.nmodes,self.nbands,self.nbands])
            database.close()

    def read_DB(self,only_freqs=False):
        """ 
        Load all the database data to memory
        """
        self.ph_energies  = np.zeros([self.nfrags,self.nmodes])
        self.ph_eigenvectors = np.zeros([self.nfrags,self.nmodes,self.natoms,3],dtype=np.complex64)
        self.gkkp = np.zeros([self.nfrags,self.nkpoints,self.nmodes,self.nbands,self.nbands],dtype=np.complex64)
        
        for iq in range(self.nfrags):
            fil = self.frag_filename + "%d"%(iq+1)
            database = Dataset(fil)
            self.ph_energies[iq] = np.sqrt(database.variables['PH_FREQS%d'%(iq+1)][:])
            eigs_q = database.variables['POLARIZATION_VECTORS'][:].T
            self.ph_eigenvectors[iq] = eigs_q[0,:,:,:] + eigs_q[1,:,:,:]*I
            gkkp = database.variables['ELPH_GKKP_Q%d'%(iq+1)][:]
            self.gkkp[iq] = (gkkp[:,:,:,:,0] + I*gkkp[:,:,:,:,1]).reshape([self.nkpoints,self.nmodes,self.nbands,self.nbands])
            database.close()
    
    def plot_elph(self,ib=1,inu=-1,cmap='viridis',size=300):
        """
        Scatterplot in the BZ of the quantity G_{nk} = 1/N_q * \sum_{q,nu} | elph_{qnu,knn} |^2 .
    
        Band and kpoint indices are user-specified.
    
        It is possible to plot the contribution of a single phonon mode specifying inu.
        """
        # Global plot stuff
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure(figsize=(12,12))
        color_map = plt.get_cmap(cmap)
        
        # Check user-provided input
        if ib >= self.nbands:   raise ValueError("The BAND index %d is not present."%ib)
        if inu > self.nmodes:   raise ValueError("The MODE index %d is not present."%inu)
    
        # Get raw data to plot
        kx, ky, kz = self.car_kpoints.T
        data = self.read_elph(ib1=ib,ib2=ib)
        
        # Prepare function G_{nk}
        to_plot = np.zeros(self.nkpoints)
        for iq in range(self.nfrags): #[ATTENTION] the sum will be performed on the AVAILABLE qpts.
            if inu>-1: 
                to_plot += np.abs(data[iq,:,inu])**2.
            else:
                for i_mode in range(self.nmodes): to_plot += np.abs(data[iq,:,i_mode])**2.
        to_plot = to_plot/self.nqpoints
        norm_to_plot = to_plot/max(to_plot)
            
        # Plot format and layout
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlabel('kx')
        ax.set_ylabel('ky')
        ax.set_zlabel('kz')
        ax.set_aspect('equal')
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        xlim, ylim, zlim = max(kx), max(ky), max(kz)
        ax.set_xlim(-xlim,xlim)
        ax.set_ylim(-ylim,ylim)
        ax.set_ylim(-zlim,zlim)
        
        # Plot title
        if inu>-1: ax.set_title('G^l_nk = 1/N_q * \sum_q | elph^ql_knn |^2')
        else: ax.set_title('G_nk = 1/N_q * \sum_ql | elph^ql_knn |^2')
        
        # Reciprocal lattice vectors
        lx,ly,lz = [ np.linalg.norm(self.rlat[i]) for i in range(3) ]
        ax.quiver(0., 0., 0., self.rlat[0,0], self.rlat[0,1], self.rlat[0,2], length=lx, normalize=True, color='black')
        ax.quiver(0., 0., 0., self.rlat[1,0], self.rlat[1,1], self.rlat[1,2], length=ly, normalize=True, color='black')
        ax.quiver(0., 0., 0., self.rlat[2,0], self.rlat[2,1], self.rlat[2,2], length=lz, normalize=True, color='black')
        
        # Gamma point
        ax.scatter(0.,0.,0.,marker='*',s=2./3.*size,edgecolors='black',facecolors='black', zorder=-1)
        
        # Actual plot
        plot = ax.scatter(kx,ky,kz,marker='o',s=size*norm_to_plot,edgecolors='black',c=norm_to_plot,cmap=color_map,zorder=1)
        fig.colorbar(plot)
    
    def __str__(self):

        try: self.ph_energies
        except AttributeError: self.read_frequencies()

        try: self.ph_eigenvectors
        except AttributeError: self.read_eigenmodes()
            
        s = 'nqpoints: %d\n'%self.nqpoints
        s+= 'nkpoints: %d\n'%self.nkpoints
        s+= 'nmodes: %d\n'%self.nmodes
        s+= 'natoms: %d\n'%self.natoms
        s+= 'nbands: %d\n'%self.nbands
        if self.nfrags == self.nqpoints: s+= 'fragments: %d\n'%self.nfrags
        else: s+= 'fragments: %d [WARNING] nfrags < nqpoints\n'%self.nfrags
        s+= '-----------------------------------\n'
        for iq in range(self.nfrags):
            s+= 'nqpoint %d\n'%iq
            for n,mode in enumerate(self.ph_eigenvectors[iq]):
                s+= 'mode %d freq: %lf meV\n'%(n,self.ph_energies[iq,n]*1000.)
                for a in range(self.natoms):
                    s += ("%12.8lf "*3+'\n')%tuple(mode[a].real)
        s+= '-----------------------------------\n'
        return s
