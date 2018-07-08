# Copyright (c) 2018, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
import os
import numpy as np
from yambopy.units import ha2ev
from yamboparser import YamboFile

class YamboQPDB():
    """
    Class to read yambo ndb.QP files

    These files describe the quasiparticle states calculated from yambo
    Includes the quasi-particl energies, the lifetimes and the Z factors
    """
    def __init__(self,filename='ndb.QP',folder='.'):
        """
        Read a QP file using the yamboparser
        """
        self.folder = folder
        self.filename = filename
        if os.path.isfile('%s/%s'%(folder,filename)):
            self.yfile = YamboFile(filename,folder)
        else:
            raise IOError('File %s/%s not found'%(folder,filename))

        self.qps          = self.yfile.data
        self.kpoints      = np.array(self.qps['Kpoint'])
        self.kpoint_index = np.array(np.array(self.qps['Kpoint_index']),dtype=int)
        self.band_index   = np.array(np.array(self.qps['Band'],dtype=int))
        self.e0           = self.qps['Eo'].real
        self.e            = self.qps['E'].real
        self.linewidths   = self.qps['E'].imag

        #read the database
        self.eigenvalues_qp, self.eigenvalues_dft, self.lifetimes = self.get_qps()

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

    def get_bs():
        """
        Get bandstructure
        """

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
    
    def qp_bs(self,lattice,path,debug=False):
        """
        Calculate quasi-particle band-structure
        """
        #get full kmesh
        kpoints = lattice.red_kpoints
        path = np.array(path)

        reps = list(range(-1,2))
        kpoints_rep, kpoints_idx_rep = replicate_red_kmesh(kpoints,repx=reps,repy=reps,repz=reps)
        band_indexes = get_path(kpoints_rep,path)
        band_kpoints = kpoints_rep[band_indexes]
        band_indexes = kpoints_idx_rep[band_indexes]

        if debug:
            for i,k in zip(band_indexes,band_kpoints):
                x,y,z = k
                plt.text(x,y,i)
            plt.scatter(kpoints_rep[:,0],kpoints_rep[:,1])
            plt.plot(path[:,0],path[:,1],c='r')
            plt.scatter(band_kpoints[:,0],band_kpoints[:,1])
            plt.show()
            exit()

        #get eigenvalues along the path
        #expand eigenvalues to the bull brillouin zone
        energies_qp = self.eigenvalues_qp[lattice.kpoints_indexes]
        #energies_qp = self.eigenvalues_qp

        #expand the quasiparticle energies to the bull brillouin zone
        energies_dft = self.eigenvalues_dft[lattice.kpoints_indexes]
        #energies_dft = self.eigenvalues_dft

        energies_dft = energies_dft[band_indexes]
        energies_qp  = energies_qp[band_indexes]

        return np.array(band_kpoints), energies_dft, energies_qp

    def plot_qp_bs(self,ax,lattice,path,what='DFT,QP',debug=False,label=False,**args):
        """
        Calculate the quasiparticle band-structure
        """
        bands_kpoints, energies_dft, energies_qp = self.qp_bs(lattice, path, debug)

        #calculate distances
        bands_distances = calculate_distances(bands_kpoints)

        #make the plots
        for b in range(self.min_band-1,self.max_band-1):
            if 'DFT' in what: 
                ax.plot(bands_distances, energies_dft[:,b], **args)
            if 'QP' in what:
                ax.plot(bands_distances, energies_qp[:,b],  **args)

        if 'DFT' in what: 
            ax.plot(bands_distances, energies_dft[:,self.max_band-1], label=label, **args)
        if 'QP' in what: 
            ax.plot(bands_distances, energies_qp[:,self.max_band-1],  label=label, **args)

        #add high-symmetry k-points vertical bars
        kpath_car = red_car(path,lattice.rlat)
        #calculate distances for high-symmetry points
        kpath_distances = calculate_distances( path ) 
        for d in kpath_distances:
            ax.axvline(d,c='k')

        xmin = np.min(bands_distances)
        xmax = np.max(bands_distances)
        plt.xlim([xmin,xmax])
        return kpath_distances

    def __str__(self):
        lines = []; app = lines.append
        app("nqps:     %d"%self.nqps)
        app("nkpoints: %d"%self.nkpoints)
        app("nbands:   %d"%self.nbands)
        return "\n".join(lines)

