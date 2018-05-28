# Copyright (c) 2017, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
from yambopy import *
from netCDF4 import Dataset
from math import sqrt
import numpy as np
from cmath import exp

I = complex(0,1)
ha2ev  = 27.211396132
ev2cm1 = 8065.54429
abs2 = lambda x: x.real**2 + x.imag**2

class YamboElectronPhononDB():
    """
    Python class to read the electron-phonon matrix elements from yambo
    """
    def __init__(self,lattice,filename='ndb.elph_gkkp',save='SAVE',only_freqs=False):
        self.lattice = lattice
        
        self.save = save
        self.filename = "%s/%s"%(save,filename)
        self.ph_eigenvalues = None
        
        self.car_kpoints = lattice.car_kpoints
        self.red_kpoints = lattice.red_kpoints
        
        #read dimensions of electron phonon parameters
        try:
            database = Dataset(self.filename)
        except:
            print "error opening %s in YamboElectronPhononDB"%self.filename
            exit()
            
        self.qpoints = database.variables['PH_Q'][:].T
        self.nmodes, self.nqpoints, self.nkpoints, self.nbands = database.variables['PARS'][:4].astype(int)
        self.natoms = self.nmodes/3
        database.close()
        
        self.readDB()

    def get_elphon(self,dir=0):
        if self.gkkp is None:
            self.get_elphon_databases()

        kpts, nks, nss = self.expand_kpts()
        gkkp = self.gkkp

        return gkkp, kpts

    def readDB(self,only_freqs=False):
        """ 
        Load all the gkkp databases to memory
        """

        self.ph_eigenvalues  = np.zeros([self.nqpoints,self.nmodes])
        self.ph_eigenvectors = np.zeros([self.nqpoints,self.nmodes,self.nmodes/3,3],dtype=np.complex64)
        if not only_freqs:
            self.gkkp = np.zeros([self.nqpoints,self.nkpoints,self.nmodes,self.nbands,self.nbands],dtype=np.complex64)
        
        for nq in xrange(self.nqpoints):
            filename = '%s_fragment_%d'%(self.filename,nq+1)

            database = Dataset(filename)

            self.ph_eigenvalues[nq] = np.sqrt(database.variables['PH_FREQS%d'%(nq+1)][:])

            p_re = database.variables['POLARIZATION_VECTORS_REAL'][:].T
            p_im = database.variables['POLARIZATION_VECTORS_IMAG'][:].T
            self.ph_eigenvectors[nq] = p_re + p_im*I
            
            if not only_freqs:
                gkkp = database.variables['ELPH_GKKP_Q%d'%(nq+1)][:]
                self.gkkp[nq] = (gkkp[:,0,:,:] + I*gkkp[:,1,:,:]).reshape([self.nkpoints,self.nmodes,self.nbands,self.nbands])
            
            database.close()

        if not only_freqs:
            return self.gkkp
        
    def __str__(self):
        if self.ph_eigenvalues is None:
            self.get_elphon_databases()
        s = 'nqpoints: %d\n'%self.nqpoints
        s+= 'nkpoints: %d\n'%self.nkpoints
        s+= 'nmodes: %d\n'%self.nmodes
        s+= 'natoms: %d\n'%self.natoms
        s+= 'nbands: %d\n'%self.nbands
        for nq in xrange(self.nqpoints):
            s+= 'nqpoint %d\n'%nq
            for n,mode in enumerate(self.ph_eigenvectors[nq]):
                s+= 'mode %d freq: %lf cm-1\n'%(n,self.ph_eigenvalues[nq][n]*ha2ev*ev2cm1)
                for a in xrange(self.natoms):
                    s += ("%12.8lf "*3+'\n')%tuple(mode[a].real)
        return s

if __name__ == '__main__':
    elph = ElectronPhononDB()
    print elph
    elph.get_databases()
