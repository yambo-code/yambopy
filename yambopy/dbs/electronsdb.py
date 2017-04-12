# Copyright (c) 2016, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
from netCDF4 import Dataset
import numpy as np
from itertools import product
import collections
ha2ev  = 27.211396132
max_exp = 50

def fermi(e):
    """ fermi dirac function
    """
    if e > max_exp:
        return 0
    elif e < -max_exp:
        return 1
    return 1/(np.exp(e)+1)

def fermi_array(e_array,ef,invsmear):
    """ 
    Fermi dirac function for an array
    """
    e_array = (e_array-ef)/invsmear
    return [ fermi(e) for e in e_array]


class YamboElectronsDB():
    def __init__(self,lattice,save='SAVE',filename='ns.db1'):
        self.lattice = lattice
        self.filename = '%s/%s'%(save,filename)
        self.efermi = None
        self.readDB()
        if self.nkpoints != self.lattice.nkpoints: #sanity check
            raise ValueError("The number of k-points in the lattice database and electrons database is different.")
        self.expandEigenvalues()
        
    def readDB(self):
        try:
            database = Dataset(self.filename)
        except:
            raise IOError("Error opening file %s in YamboElectronsDB"%self.filename)

        self.eigenvalues_ibz  = database.variables['EIGENVALUES'][0,:]*ha2ev
        self.iku_kpoints      = database.variables['K-POINTS'][:].T
        dimensions = database.variables['DIMENSIONS'][:]
        self.nbands      = dimensions[5]
        self.temperature = dimensions[13]
        self.nelectrons  = int(dimensions[14])
        self.nkpoints    = int(dimensions[6])
        self.nbands      = int(dimensions[5])
        self.spin = int(dimensions[11])
        self.time_rev = dimensions[9]
        database.close()
        
        #spin degeneracy if 2 components degen 1 else degen 2
        self.spin_degen = [0,2,1][int(self.spin)]
        
        #number of occupied bands
        self.nbandsv = self.nelectrons / self.spin_degen
        self.nbandsc = self.nbands-self.nbandsv

    def expandEigenvalues(self):
        """
        Expand eigenvalues to the full brillouin zone
        """
        
        self.eigenvalues = self.eigenvalues_ibz[self.lattice.kpoints_indexes]
        
        self.nkpoints_ibz = len(self.eigenvalues_ibz)
        self.weights_ibz = np.zeros([self.nkpoints_ibz],dtype=np.float32)
        self.nkpoints = len(self.eigenvalues)
        
        #counter counts the number of occurences of element in a list
        for nk_ibz,inv_weight in collections.Counter(self.lattice.kpoints_indexes).items():
            self.weights_ibz[nk_ibz] = float(inv_weight)/self.nkpoints
        
        #kpoints weights
        self.weights = np.full((self.nkpoints), 1.0/self.nkpoints,dtype=np.float32)
 
    def setFermi(self,fermi,invsmear):
        """
        Set the fermi energy of the system
        """
        self.invsmear = invsmear
        self.efermi = fermi
        
        #full brillouin zone
        self.eigenvalues     -= self.efermi
        self.occupations = np.zeros([self.nkpoints,self.nbands],dtype=np.float32)
        for nk in xrange(self.nkpoints):
            self.occupations[nk] = fermi_array(self.eigenvalues[nk,:],0,self.invsmear)
        
        #for the ibz
        self.eigenvalues_ibz -= self.efermi
        self.occupations_ibz = np.zeros([self.nkpoints_ibz,self.nbands],dtype=np.float32)
        for nk in xrange(self.nkpoints_ibz):
            self.occupations_ibz[nk] = fermi_array(self.eigenvalues_ibz[nk,:],0,self.invsmear)
        
        return self.efermi
        
    def setFermiFixed(self,broad=1e-5):
        """
        Get fermi level using fixed occupations method
        Useful for semi-conductors
        """
        eigenvalues = self.eigenvalues_ibz
        weights     = self.weights_ibz
        nkpoints    = self.nkpoints_ibz

        nbands = self.nelectrons/self.spin_degen
        #top of valence
        top = np.max(eigenvalues[:,nbands])
        #bottom of conduction
        bot = np.max(eigenvalues[:,nbands-1])
        self.efermi = (top+bot)/2
        self.setFermi(self.efermi,broad)
    
    def getFermi(self,invsmear,setfermi=True):
        """ 
        Determine the fermi energy
        """
        if self.efermi: return self.efermi
        from scipy.optimize import bisect
        
        eigenvalues = self.eigenvalues_ibz
        weights     = self.weights_ibz
        nkpoints    = self.nkpoints_ibz
        
        min_eival, max_eival = np.min(eigenvalues), np.max(eigenvalues)
        self.invsmear = invsmear
        
        def occupation_minus_ne(ef):
            """ 
            The total occupation minus the total number of electrons
            """
            return sum([sum(self.spin_degen*fermi_array(eigenvalues[nk],ef,self.invsmear))*weights[nk] for nk in xrange(nkpoints)])-self.nelectrons

        fermi = bisect(occupation_minus_ne,min_eival,max_eival)
        if setfermi: self.setFermi(fermi,invsmear)
        return self.efermi

    def __str__(self):
        s = ""
        s += "spin_degen: %d\n"%self.spin_degen
        s += "nelectrons: %d\n"%self.nelectrons
        s += "nbands:   %d\n"%self.nbands
        s += "nkpoints: %d"%self.nkpoints
        return s
