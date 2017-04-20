# Copyright (c) 2017, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
from builtins import range
from builtins import object
from yambopy import *
from yamboparser import *

class YamboQPDB(object):
    """
    Class to read yambo ndb.QP files
    
    These files describe the quasiparticle states calculated from yambo
    Includes the quasi-particl energies, the lifetimes and the Z factors
    """
    def __init__(self,filename='ndb.QP',folder='.'):
        """
        Read a QP file using the yamboparser 
        """
        self.yfile = YamboFile('%s/%s' %(folder,filename))

        qps = self.yfile.data
        self.qps   = qps
        self.nqps  = len(qps['E'])
        self.nkpoints = len(qps['Kpoint'])

        #get kpoints
        kpts=[]
        for nk in range(self.nkpoints):
            kpts.append(qps['Kpoint'][nk])
        self.kpoints = np.array(kpts)

        #get nbands
        min_band = int(qps['Band'][0])
        max_band = int(qps['Band'][0])
        for iqp in range(self.nqps):
            band  = int(qps['Band'][iqp])
            if min_band > band: min_band = band
            if max_band < band: max_band = band
        self.nbands = max_band-min_band+1

    def get_qps(self,):
        """
        Get quasiparticle energies in a list
        """
        #get dimensions
        nqps   = self.nqps
        nkpts  = self.nkpoints
        nbands = self.nbands

        qps  = self.qps
        kpts = self.kpoints
        
        #start arrays 
        eigenvalues_lda = np.zeros([nkpts,nbands])
        eigenvalues_qp  = np.zeros([nkpts,nbands])
        lifetimes       = np.zeros([nkpts,nbands])
        for iqp in range(nqps):
            kindx = int(qps['Kpoint_index'][iqp])
            e     = qps['E'][iqp]
            e0    = qps['Eo'][iqp]
            band  = int(qps['Band'][iqp])
            kpt   = ("%8.4lf "*3)%tuple(kpts[kindx-1])
            Z     = qps['Z'][iqp]
            eigenvalues_qp[kindx-1,band-1] = e.real
            eigenvalues_lda[kindx-1,band-1] = e0.real
            lifetimes[kindx-1,band-1] = e.imag

        return eigenvalues_qp, eigenvalues_lda, lifetimes

    def __str__(self):
        s = ""
        s += "nqps:     %d\n"%self.nqps
        s += "nkpoints: %d\n"%self.nkpoints
        s += "nbands:   %d\n"%self.nbands
        return s        
        
