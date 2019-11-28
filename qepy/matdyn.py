# Copyright (C) 2018 Henrique Pereira Coutada Miranda, Alejandro Molina Sanchez
# All rights reserved.
#
# This file is part of yambopy
#
from __future__ import print_function, division
import os
import re
from math import sqrt
import numpy as np
from .lattice import *
from qepy.auxiliary import *

eVtocm1 = 8065.54429
cm1toeV = 1.0/eVtocm1
ha2ev  = 27.211396132
ev2ha  = 1.0/ha2ev
Thz2cm1 = 33.35641
cm12Thz = 1.0/33.35641

__all__ = [
"Matdyn",
]

class Matdyn(object):
    """
    Class to read and plot the data from matdyn.modes files 
    """

    def __init__(self,qpoints,eig,eiv):
        self.qpoints  = np.array(qpoints)
        self.eig      = np.array(eig)
        self.eiv      = np.array(eiv)

    @classmethod
    def from_modes_file(cls,folder='.',filename='matdyn.modes'):
        """
        read the modes file from the hard drive 
        """
        with open(os.path.join(folder, filename),'r') as f:
            data_phon = f.readlines()

        #detect dimensions of the file
        #qpoints
        nqpoints = 0
        for line in data_phon:
            if 'q =' in line:
                nqpoints+=1
        nqpoints = nqpoints

        #modes
        nline = 0
        while 'freq' not in data_phon[nline]:
            nline += 1
        start_line = nline
        nline +=1
        while 'freq' not in data_phon[nline]:
            nline += 1
        end_line = nline
        natoms = end_line-start_line-1
        nmodes = natoms*3

        #empty stuff
        eig = []
        eiv = []
        qpoints = []

        #read qpoints, modes and energies
        for j in range(nqpoints):
            frec, v_frec = [], []
            k=2 + j*(nmodes*(natoms+1)+5)
            qpoints.append( float_from_string(data_phon[k]) )
            for i in range(nmodes):
                k=4 + j*(nmodes*(natoms+1)+5) + i*(natoms+1)
                y = float_from_string(data_phon[k])
                v_mode = []
                for ii in range(1,natoms+1):
                    z      = float_from_string(data_phon[k+ii])
                    v_atom = [complex(z[0],z[1]),complex(z[2],z[3]),complex(z[4],z[5])]
                    v_mode.append(v_atom)
                v_frec.append(v_mode)
                frec.append(y[1])
            eig.append(frec)
            eiv.append(v_frec)

        #store info
        natoms  = natoms
        nmodes  = nmodes
        nqpoints= nqpoints
        qpoints = np.array(qpoints)
        eig     = np.array(eig)
        eiv     = np.array(eiv).reshape(nqpoints,nmodes,nmodes)
        return cls(qpoints,eig,eiv)

    @property
    def modes(self):
        return self.eiv.reshape(self.nqpoints,self.nmodes,self.natoms,3)

    @property
    def nmodes(self):
        nqpoints,nmodes,_ = self.eiv.shape 
        return nmodes
 
    @property
    def natoms(self):
        nqpoints,nmodes,_ = self.eiv.shape 
        return int(nmodes/3)

    @property
    def nqpoints(self):
        nqpoints,nmodes,_ = self.eiv.shape 
        return nqpoints

    @property
    def write_modes(self,filename=None):
        """
        save the phonon modes in a file
        """
        s = " matrix written with qepy\n\n"
        for nq in range(self.nqpoints):
            s += ("q =  "+"%12.6lf "*3+"\n")%tuple(self.qpoints[nq])
            s += "*"*81+"\n"
            for n,mode in enumerate(self.eiv[nq]):
                phfreqmev = self.get_phonon_freq(nq,n+1,unit='THz')
                phfreqcm1 = self.get_phonon_freq(nq,n+1,unit='cm-1')
                s += '    freq ( %4d) = %12.6lf [ThZ] = %12.6lf [cm-1]\n'%(n+1,phfreqmev,phfreqcm1)
                for a in range(self.natoms):
                    xr,yr,zr = mode[a*3:(a+1)*3].real
                    xi,yi,zi = mode[a*3:(a+1)*3].imag
                    s += ("( "+"%12.6lf "*6+')\n')%(xr,xi,yr,yi,zr,zi)
            s += "*"*81+"\n"
            s += "\n\n"
        
        if filename:
            f = open(filename,'w')
            f.write(s)
            f.close()
        else:
            print(s)

    def rotate_phonons(self,eps=1e-5,debug=False):
        """
        Rotate the degenerate states to align then with the x,y axis
        Arguments:
        eps -> threshold for finding degeneracies
        """
        eig = self.eig
        eiv = self.eiv
        
        #the eigenvalues are probably sorted but just in case...
        eig, eiv = zip(*sorted(zip(eig,eiv), key=lambda x: x[0]))
        eig = np.array(eig)
        eiv = np.array(eiv)
        
        #iterate over qpoints
        for nq,(eigq,eivq) in enumerate(zip(eig,eiv)):
            #detect the degeneracies
            keys = [group.mean() for group in np.split(eigq, np.where(np.diff(eigq) > eps)[0]+1)]
            unique_frequencies = dict([(key,0) for key in keys])
            #count degeneracies
            for eig in eigq:
                for key in keys:
                    if np.isclose(key,eig,atol=eps):
                        unique_frequencies[key] += 1
            for val,deg in sorted(unique_frequencies.items()):
                if debug: print("frequency: %12.8lf [Thz] %12.8lf [cm-1] degeneracy: %2d" % (val,val,deg))
                if deg > 1:
                    #get the indexes of the modes with the same frequency
                    indexes = [n for n,f in enumerate(eigq) if np.isclose(f,val,atol=eps)]
                    r = np.array(eivq[indexes])
                    if debug:
                        print("input basis:")
                        for n in range(deg):
                            print("mode: %3d"%n)
                            for i in range(self.natoms):
                               print("atom %3d"%i+("%12.8lf"*3)%tuple(r[n,i*3:(i+1)*3].real))
                    #we make sure the first column vector the matrix r in non zero
                    rows,cols = r.shape
                    n = 0
                    while np.isclose(r[:,n],np.zeros([len(r)])).all():
                        n += 1
                    q, a = np.linalg.qr(r[:,n:])
                    r[:,n:] = a
                    if debug:
                        print("canonical basis:")
                        for n in range(deg):
                            print("mode: %3d"%n)
                            for i in range(self.natoms):
                                print("atom %3d"%i+("%12.8lf"*3)%tuple(r[n,i*3:(i+1)*3].real))
                    eivq[indexes] = r
                self.eiv[nq] = eivq

    def plot_eigen(self,path=[]):
        """ plot the phonon frequencies using matplotlib
        """
        import matplotlib.pyplot as plt

        if path:
            if isinstance(path,Path):
                path = path.get_indexes()
            plt.xticks( *list(zip(*path)) )
        plt.ylabel('\\omega (cm$^{-1}$)')

        #plot vertical line
        for point in path:
            x, label = point
            plt.axvline(x)

        #plot bands
        eig = np.array(self.eig)
        for ib in range(self.nmodes):
           plt.plot(range(self.nqpoints),eig[:,ib], 'r-', lw=2)
        plt.show()

    def get_phonon_freq(self,nq,n,unit="eV"):
        """
        Get the value of the phonon frequency
        nq -> q-point from where to get the frequency from
        n  -> mode of the phonon
        """
        if   unit == "eV":
            factor = cm1toeV
        elif unit == "Ha":
            factor = cm1toeV*eV2ha
        elif unit == "THz":
            factor = cm12Thz
        elif unit == "cm-1":
            factor = 1
        else:
            raise ValueError('Unit %s not known'%unit)

        return self.eig[nq][n-1]*factor

    def normalize(self):
        """
        Normalize the displacements u^n_{ai} according to:
        
        sum_ai ( u^n_{ai} )**2 = 1
        """
       
        for nq in range(self.nqpoints):
            for n in range(self.nmodes):
                print(np.linalg.norm(self.eiv[nq,n]))
                self.eiv[nq,n] /= np.linalg.norm(self.eiv[nq,n])

    def normalize_with_masses(self,masses): 
        """
        Normalize the displacements u^n_{ai} according to:
        
        sum_{ai} M_a u^n_{ai} u^m_{ai} = delta_{nm}
        
        u -> displacement
        n -> phonon mode
        a -> atom index
        i -> direction
        M -> mass
        """

        masses = np.array(masses)
        ref_mass = max(masses)
        masses = masses/ref_mass

        #divide by masses
        if self.check_orthogonality():
            for nq in range(self.nqpoints):
                for n in range(self.nmodes):
                    for a in range(self.natoms):
                        self.eiv[nq,n,a*3:(a+1)*3] *= 1.0/sqrt(masses[a])
        else:
            print("These eigenvectors are non-orthogonal, probably they are already scaled by the masses so I won't do it")

        #enforce delta_nm
        for nq in range(self.nqpoints):
            for n in range(self.nmodes):
                s = 0
                for a in range(self.natoms):
                    e = self.eiv[nq,n,a*3:(a+1)*3]
                    #get normalization constant
                    s += masses[a]*np.vdot(e,e).real
                self.eiv[nq,n] *= 1/sqrt(s)
    
    def check_orthogonality(self,atol=1e-5):
        """
        Check if the eigenvectors are orthogonal
        """

        orth = np.zeros([self.nmodes,self.nmodes])
        for nq in range(self.nqpoints):
            for n in range(self.nmodes):
                e1 = self.eiv[nq,n]
                for m in range(self.nmodes):
                    e2 = self.eiv[nq,m]
                    orth[n,m] = np.vdot(e1,e2).real
        
        return np.isclose(orth,np.eye(self.nmodes),atol=atol).all()
 
    def check_normalization(self,masses,atol=1e-5):
        """
        Check if the displacements are normalized according to:
        
        sum_{ai} M_a u^n_{ai} u^m_{ai} = delta_{nm}
        """
        
        masses = np.array(masses)
        ref_mass = max(masses)
        masses = masses/ref_mass
        
        #check normalization
        norm = np.zeros([self.nmodes])
        for nq in range(self.nqpoints):
            for n in range(self.nmodes):
                s = 0
                for a in range(self.natoms):
                    e = self.eiv[nq,n,a*3:(a+1)*3]
                    #get normalization constant
                    s += masses[a]*np.vdot(e,e).real
                norm[n] = s

        return np.isclose(norm,np.ones(self.nmodes),atol=atol).all()
 
    def __str__(self):
        s = ""
        for nq in range(self.nqpoints):
            for n,mode in enumerate(self.eiv[nq]):
                phfreqmev = self.get_phonon_freq(nq,n+1,unit='eV')*1000
                phfreqcm1 = self.get_phonon_freq(nq,n+1,unit='cm-1')
                s+= 'mode: %d freq: %8.2lf meV %8.2lf cm-1\n'%(n+1,phfreqmev,phfreqcm1)
                for a in range(self.natoms):
                    s += ("%12.8lf "*3+'\n')%tuple(mode[a*3:(a+1)*3].real)
        return s

