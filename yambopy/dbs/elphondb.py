# Copyright (c) 2018, Fulvio Paleari, Alejandro Molina-Sánchez, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
from yambopy import *
from netCDF4 import Dataset
from math import sqrt
import numpy as np
from cmath import exp
from yambopy.units import ha2ev, ev2cm1, I

class YamboElectronPhononDB():
    """
    Python class to read the electron-phonon matrix elements from yambo.
    
    By default it reads dimension of elph parameters and qpoints.
    
    Call function read_frequencies for frequencies.
    Call function read_eigenmodes for phonon modes.
    Call function read_elph for electron-phonon matrix elements.
    Call function read_DB to read everything.
    
    Plotting functions (scatterplots) are also provided.        
    """
    def __init__(self,lattice,filename='ndb.elph_gkkp',folder_gkkp='SAVE',save='SAVE',single_q=None):
        #self.lattice = lattice
        
        #self.save = save
        self.filename = "%s/%s"%(folder_gkkp,filename)
        #self.ph_eigenvalues = None
        
        #self.car_kpoints = lattice.car_kpoints
        #self.red_kpoints = lattice.red_kpoints
        #self.rlat        = lattice.rlat
        # necessary lattice information
        self.alat        = lattice.alat
    
        try: database = Dataset(self.filename)
        except: raise FileNotFoundError("error opening %s in YamboElectronPhononDB"%self.filename)
        
        #read qpoints    
        self.qpoints = database.variables['PH_Q'][:].T
        self.car_qpoints = np.array([ q/self.alat for q in self.qpoints ])
        #read dimensions of electron phonon parameters
        self.nmodes, self.nqpoints, self.nkpoints, self.nbands = database.variables['PARS'][:4].astype(int)
        self.natoms = int(self.nmodes/3)
        
        if single_q is not None: self.qrange=[single_q,single_q+1]
        else: self.qrange=[0,self.nqpoints]
        
        database.close()

    def read_frequencies(self):
        """
        Read phonon frequencies
        """
        self.ph_energies  = np.zeros([self.nqpoints,self.nmodes])
        qrange = self.qrange
        
        for iq in range(qrange[0],qrange[1]):
            filnm = self.filename + "_fragment_%d"%(iq+1)
            database = Dataset(filnm)
            self.ph_energies[iq] = np.sqrt(database.variables['PH_FREQS%d'%(iq+1)][:])*ha2ev
            database.close()
        
    def read_eigenmodes(self):
        """
        Read phonon eigenmodes
        """
        self.ph_eigenvectors = np.zeros([self.nqpoints,self.nmodes,self.natoms,3],dtype=np.complex64)
        qrange = self.qrange
        
        for iq in range(qrange[0],qrange[1]):
            filnm = self.filename + "_fragment_%d"%(iq+1)
            database = Dataset(self.filename)
            #eigs_q[cartesian][atom][mode][complex]
            eigs_q = database.variables['POLARIZATION_VECTORS']         #[:].T
            self.ph_eigenvectors[iq] = eigs_q[:,:,:,0] + eigs_q[:,:,:,1]*I
            database.close()
        
    def read_elph(self):
        """
        Read electron-phonon matrix elements
        
        TODO: Add memory-saving option to read only across some dimensions (or just a specific matrix element)
        """
        # gkkp[q][k][mode][bnd1][bnd2]
        self.gkkp = np.zeros([self.nqpoints,self.nkpoints,self.nmodes,self.nbands,self.nbands],dtype=np.complex64)
        qrange = self.qrange
        
        for iq in range(qrange[0],qrange[1]):
            filnm = self.filename + "_fragment_%d"%(iq+1)
            database = Dataset(self.filename)
            gkkp = database.variables['ELPH_GKKP_Q%d'%(iq+1)][:]
            self.gkkp[iq] = (gkkp[:,0,:,:] + I*gkkp[:,1,:,:]).reshape([self.nkpoints,self.nmodes,self.nbands,self.nbands])
            database.close()

    def read_DB(self,only_freqs=False):
        """ 
        Load all the database data to memory
        """
        self.ph_energies  = np.zeros([self.nqpoints,self.nmodes])
        self.ph_eigenvectors = np.zeros([self.nqpoints,self.nmodes,self.nmodes/3,3],dtype=np.complex64)
        self.gkkp = np.zeros([self.nqpoints,self.nkpoints,self.nmodes,self.nbands,self.nbands],dtype=np.complex64)
        qrange = self.qrange
        
        for iq in range(qrange[0],qrange[1]):
            filnm = self.filename + "_fragment_%d"%(iq+1)
            database = Dataset(self.filename)
            self.ph_energies[iq] = np.sqrt(database.variables['PH_FREQS%d'%(iq+1)][:])
            eigs_q = database.variables['POLARIZATION_VECTORS'][:].T
            self.ph_eigenvectors[iq] = eigs_q[:,:,:,0] + eigs_q[:,:,:,1]*I
            gkkp = database.variables['ELPH_GKKP_Q%d'%(iq+1)][:]
            self.gkkp[iq] = (gkkp[:,0,:,:] + I*gkkp[:,1,:,:]).reshape([self.nkpoints,self.nmodes,self.nbands,self.nbands])
            database.close()
        
    def __str__(self):

        try: 
            self.ph_energies
        except AttributeError:
            self.read_frequencies()
            self.read_eigenmodes()
            
        s = 'nqpoints: %d\n'%self.nqpoints
        s+= 'nkpoints: %d\n'%self.nkpoints
        s+= 'nmodes: %d\n'%self.nmodes
        s+= 'natoms: %d\n'%self.natoms
        s+= 'nbands: %d\n'%self.nbands
        for nq in range(self.nqpoints):
            s+= 'nqpoint %d\n'%nq
            for n,mode in enumerate(self.ph_eigenvectors[nq]):
                s+= 'mode %d freq: %lf cm-1\n'%(n,self.ph_energies[nq][n])
                for a in range(self.natoms):
                    s += ("%12.8lf "*3+'\n')%tuple(mode[a].real)
        return s