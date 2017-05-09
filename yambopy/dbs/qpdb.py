# Copyright (c) 2017, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
from builtins import range
from builtins import object
from yambopy import *
from yamboparser import *
import os

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
        if os.path.isfile('%s/%s'%(folder,filename)):
            self.yfile = YamboFile(filename,folder)
        else:
            raise ValueError('File %s/%s not found'%(folder,filename))

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

        #read the database
        self.eigenvalues_qp, self.eigenvalues_dft, self.lifetimes = self.get_qps()

    def get_qps(self,):
        """
        Get quasiparticle energies in a list
        """
        #get dimensions
        nqps   = self.nqps
        nkpts  = self.nkpoints

        qps  = self.qps
        kpts = self.kpoints
        nbands = int(np.max(qps['Band'][:]))

        #start arrays
        eigenvalues_dft = np.zeros([nkpts,nbands])
        eigenvalues_qp  = np.zeros([nkpts,nbands])
        lifetimes       = np.zeros([nkpts,nbands])
        for iqp in range(nqps):
            kindx = int(qps['Kpoint_index'][iqp])
            e     = qps['E'][iqp]*ha2ev
            e0    = qps['Eo'][iqp]*ha2ev
            band  = int(qps['Band'][iqp])
            kpt   = ("%8.4lf "*3)%tuple(kpts[kindx-1])
            Z     = qps['Z'][iqp]
            eigenvalues_qp[kindx-1,band-1] = e.real
            eigenvalues_dft[kindx-1,band-1] = e0.real
            lifetimes[kindx-1,band-1] = e.imag

        return eigenvalues_qp, eigenvalues_dft, lifetimes

    def __str__(self):
        s = ""
        s += "nqps:     %d\n"%self.nqps
        s += "nkpoints: %d\n"%self.nkpoints
        s += "nbands:   %d\n"%self.nbands
        return s

