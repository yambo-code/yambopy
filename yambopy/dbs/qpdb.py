# Copyright (c) 2018, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
from yambopy import *
from yamboparser import *
import os

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
        self.min_band = min_band
        self.max_band = max_band
        self.nbands = max_band-min_band+1

        #read the database
        self.eigenvalues_qp, self.eigenvalues_dft, self.lifetimes = self.get_qps()

    def qp_bs(self,lattice,path,debug=False):
        """
        Calculate qusi-particle band-structure
        """
        #get full kmesh
        kpoints = lattice.red_kpoints
        path = np.array(path)

        kpoints_rep, kpoints_idx_rep = replicate_red_kmesh(kpoints,repx=list(range(-1,2)),repy=list(range(-1,2)),repz=list(range(-1,2)))
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

    def get_qps(self):
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

