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

class YamboElectronsDB():
    def __init__(self,lattice,save='SAVE',filename='ns.db1',expand=True):
        self.lattice = lattice
        self.filename = '%s/%s'%(save,filename)
        self.efermi = None
        self.readDB()
        if self.nkpoints != self.lattice.nkpoints: #sanity check
            raise ValueError("The number of k-points in the lattice database and electrons database is different.")
        if expand: self.expandEigenvalues()
        
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

    def getEcMinusEv(self):
        """
        Calculate G=ec-ev
        """
        #calculate G
        G = np.zeros([self.nkpoints,self.nbandsv,self.nbandsc])
                
        #iterate over the kpoints:
        for nk in xrange(self.nkpoints):
            
            #get eigenvalues for this kpoint
            eivk = self.eigenvalues[nk]
            
            #iterate over valence and conduction bands
            for v in xrange(self.nbandsv):
                    
                #calculate ec-ev
                G[nk,v,:] = eivk[self.nbandsv:]-eivk[v]
                
        return G
                
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
        
        
        
