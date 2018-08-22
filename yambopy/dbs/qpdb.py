# Copyright (c) 2018, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
import os
import numpy as np
from yambopy.units import ha2ev
from yambopy.tools.string import marquee
from yambopy.plot.bandstructure import YamboBandStructure
from yambopy.plot.plotting import add_fig_kwargs
from yamboparser import YamboFile

class YamboQPDB():
    """
    Class to read yambo ndb.QP files

    These files describe the quasiparticle states calculated from yambo
    Includes the quasi-particl energies, the lifetimes and the Z factors
    """
    def __init__(self,qps):
        """
        Initialize the YamboQP class
        """
        self.qps          = qps
        self.kpoints      = np.array(qps['Kpoint'])
        self.kpoint_index = np.array(qps['Kpoint_index'],dtype=int)
        self.band_index   = np.array(qps['Band'],dtype=int)
        self.e0           = np.array(qps['Eo']).real*ha2ev
        self.e            = np.array(qps['E']).real*ha2ev
        self.linewidths   = np.array(qps['E']).imag*ha2ev
        self.eigenvalues_qp, self.eigenvalues_dft, self.lifetimes = self.get_qps()
    
    @classmethod
    def from_db(cls,filename='ndb.QP',folder='.'):
        """
        Create instance of this class from a ndb.QP file
        """
        if os.path.isfile('%s/%s'%(folder,filename)):
            yfile = YamboFile(filename,folder)
        else:
            raise IOError('File %s/%s not found'%(folder,filename))
        return cls(yfile.data)
    
    def get_qps(self):
        """
        Get quasiparticle energies in a list
        """
        #start arrays
        eigenvalues_dft = np.zeros([self.nkpoints,self.nbands])
        eigenvalues_qp  = np.zeros([self.nkpoints,self.nbands])
        linewidths      = np.zeros([self.nkpoints,self.nbands])
        for ei,e0i,li,ki,ni in zip(self.e,self.e0,self.linewidths,self.kpoint_index,self.band_index):
            nkpoint = ki-self.min_kpoint
            nband = ni-self.min_band
            eigenvalues_dft[nkpoint,nband] = e0i
            eigenvalues_qp[nkpoint,nband] = ei
            linewidths[nkpoint,nband] = li

        return eigenvalues_dft, eigenvalues_qp, linewidths

    def get_bs(self,bs=None):
        """
        Get YamboBandStructure object with the KS and GW bands
        """
        if bs is None: bs = YamboBandStructure()

        eigenvalues_dft, eigenvalues_qp, lifetimes = self.get_qps()
        
        #add bands
        bs.add_bands(eigenvalues_dft,label='KS')
        bs.add_bands(eigenvalues_qp, label='GW')

        return bs

    @add_fig_kwargs
    def plot_bs(self):
        """
        Get and plot QP bandstructure
        """
        bs = self.get_bs()
        return bs.plot(show=False)

    @property
    def nqps(self):
        return len(self.qps)

    @property
    def min_kpoint(self):
        return min(self.kpoint_index)

    @property
    def max_kpoint(self):
        return max(self.kpoint_index)

    @property
    def nbands(self):
        return self.max_band-self.min_band+1

    @property
    def min_band(self):
        return min(self.band_index)

    @property
    def max_band(self):
        return max(self.band_index)

    @property
    def nkpoints(self):
        return len(self.kpoints)
    
    def __str__(self):
        lines = []; app = lines.append
        app(marquee(self.__class__.__name__))
        app("nqps:     %d"%self.nqps)
        app("nkpoints: %d"%self.nkpoints)
        app("nbands:   %d"%self.nbands)
        return "\n".join(lines)

