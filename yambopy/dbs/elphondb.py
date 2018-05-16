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
from pylab import *
import matplotlib.pyplot as plt

I = complex(0,1)
ha2ev  = 27.211396132
ev2cm1 = 8065.54429
abs2 = lambda x: x.real**2 + x.imag**2

class YamboElectronPhononDB():
    """
    Python class to read the electron-phonon matrix elements from yambo
    """
    def __init__(self,lattice,filename='ndb.elph_gkkp',folder_gkkp='SAVE',save='SAVE',only_freqs=False):
        self.lattice = lattice
        
        self.save = save
        self.filename = "%s/%s"%(folder_gkkp,filename)
        self.ph_eigenvalues = None
        
        self.car_kpoints = lattice.car_kpoints
        self.red_kpoints = lattice.red_kpoints
        self.rlat        = lattice.rlat
        #read dimensions of electron phonon parameters
        try:
            database = Dataset(self.filename)
        except:
            print "error opening %s in YamboElectronPhononDB"%self.filename
            exit()
            
        self.qpoints = database.variables['PH_Q'][:].T
        self.car_qpoints = np.array([ q/self.lattice.alat for q in self.qpoints ])

        self.nmodes, self.nqpoints, self.nkpoints, self.nbands = database.variables['PARS'][:4].astype(int)
        self.natoms = self.nmodes/3
        database.close()
        
        self.readDB_n_np(ib1=2,ib2=3,ik1=3)
        #self.readDB()
        #print self.gkkp

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

    def readDB_n_np(self,ib1=1,ib2=1,ik1=1):
        # Read gkkps for a given n,n' and k
        # The structure of the gkkps in Yambo is
        # GKKP(q)[k,complex,nmodes,nbands*nbands]

        iband = (ib1-1)*self.nbands + (ib2-1)
        
        self.gkkp_n_np_kn = np.zeros([self.nqpoints,self.nmodes],dtype=np.complex64)

        print 'iband', iband

        for nq in xrange(self.nqpoints):
            filename = '%s_fragment_%d'%(self.filename,nq+1)

            database = Dataset(filename)

            #self.ph_eigenvalues[nq] = np.sqrt(database.variables['PH_FREQS%d'%(nq+1)][:])

            #p_re = database.variables['POLARIZATION_VECTORS_REAL'][:].T
            #p_im = database.variables['POLARIZATION_VECTORS_IMAG'][:].T
            #self.ph_eigenvectors[nq] = p_re + p_im*I
            
            #if not only_freqs:
            self.gkkp_n_np_kn[nq] = database.variables['ELPH_GKKP_Q%d'%(nq+1)][ik1-1,0,:,iband] + I* database.variables['ELPH_GKKP_Q%d'%(nq+1)][ik1-1,1,:,iband]
            #self.gkkp_n_np_kn[nq] = (gkkp[:,0,:,:] + I*gkkp[:,1,:,:]).reshape([self.nkpoints,self.nmodes,self.nbands,self.nbands])
            
            database.close()

        return self.gkkp_n_np_kn

    def plot_map(self,fig,ib1=1,ib2=1,ik1=1,all_phonons=True,cmap='viridis',size=60,lim=0.15):
        """
        Alejandro Molina-Sanchez
        Plot the gkkp in a scatter plot (1st version developed by A. Molina-Sanchez)
        Options:
        cmap : colormap. Default viridis 
        log_scale : Logarithmic scale for the intensity (True or False) Do we put that?
        set_maximum : All plots are normalized 
        Further development: Option for the colorbar
 
        """
        import matplotlib.pyplot as plt
        import matplotlib.colors as colors
        #size=20,marker='H',set_origin=0.0,lim=0.2,cmap='viridis',log_scale=False,set_maximum=1.0

        # GKKP(k+q)[n,n',k]
        data=self.readDB_n_np(ib1,ib2,ik1)


        color_map = plt.get_cmap(cmap)
        
        kx_aux, ky_aux = self.car_qpoints[:,0], self.car_qpoints[:,1]

        kx = concatenate([kx_aux,kx_aux+self.rlat[0,0],kx_aux-self.rlat[0,0],kx_aux+self.rlat[1,0],kx_aux-self.rlat[1,0],kx_aux+self.rlat[0,0]-self.rlat[1,0],kx_aux-self.rlat[0,0]+self.rlat[1,0]])
        ky = concatenate([ky_aux,ky_aux+self.rlat[0,1],ky_aux-self.rlat[0,1],ky_aux+self.rlat[1,1],ky_aux-self.rlat[1,1],ky_aux+self.rlat[0,1]-self.rlat[1,1],ky_aux-self.rlat[0,1]+self.rlat[1,1]])
        
        """ 
        all_phonons options
        True:  Sum over all phonon modes
        False: Plot all gkkp from each phonon mode
        """
        if all_phonons:
            ax = fig.add_subplot(111)
            ax.set_aspect('equal')
            ax.axes.get_xaxis().set_visible(False)
            ax.axes.get_yaxis().set_visible(False)
            ax.set_xlim(-lim,lim)
            ax.set_ylim(-lim,lim)
            for ip in range(self.nmodes):
              gkkp_aux += abs(data[:,ip])
            max_gkkp = max(gkkp_aux)
            gkkp = concatenate(7*[gkkp_aux/max_gkkp])   
            ax.scatter( kx,ky,s=size,marker='H',c=gkkp,cmap=color_map)
        else:
            for ip in range(self.nmodes):
                square_size = 0.25
                x = 0.05 + (square_size+0.05)*(ip-ip/3*3)
                y = 0.75 - (square_size+0.05)*(ip/3)
                ax = fig.add_axes( [ x, y, square_size, square_size ])
                ax.set_aspect('equal')
                ax.axes.get_xaxis().set_visible(False)
                ax.axes.get_yaxis().set_visible(False)
                ax.set_xlim(-lim,lim)
                ax.set_ylim(-lim,lim)
                ax.set_facecolor(color_map(0.0))

                gkkp_aux = abs(data[:,ip])
                max_gkkp = max(gkkp_aux)
                gkkp = concatenate(7*[gkkp_aux/max_gkkp])   
                ax.scatter( kx,ky,s=size,marker='H',c=gkkp,cmap=color_map)

    def plot_modulus(self,ib1=1,ib2=1,ik1=1,all_phonons=True):
        # GKKP(k+q)[n,n',k]
        data=self.readDB_n_np(ib1,ib2,ik1)

        q_modulus = zeros(self.nqpoints)

        # Modulus of q-point
        for iq in range(self.nqpoints):
            q_modulus[iq] = sqrt(np.dot(self.car_qpoints[iq],self.car_qpoints[iq]))

        """ 
        all_phonons options
        True:  Sum over all phonon modes
        False: Plot all gkkp from each phonon mode
        """
        if all_phonons:
            gkkp     = zeros([self.nqpoints])
            gkkp_aux = zeros([self.nqpoints])
            for ip in range(self.nmodes):
                gkkp_aux[:] += abs(data[:,ip])
            gkkp[:] = gkkp_aux[:]#/max(gkkp_aux)
        else:
            gkkp = zeros([self.nqpoints,self.nmodes])
            for ip in range(self.nmodes):
                gkkp[:,ip] = abs(data[:,ip])
        print gkkp.shape
        # q_modulus : array dimension: nqpoints
        # gkkp      : matrix dimension: (nqpoints x nphonons) or (nqpoints)

        return q_modulus,gkkp

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
