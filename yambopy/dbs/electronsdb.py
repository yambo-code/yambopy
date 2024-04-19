#
# License-Identifier: GPL
#
# Copyright (C) 2024 The Yambo Team
#
# Authors: HPC, FP
#
# This file is part of the yambopy project
#
from netCDF4 import Dataset
import numpy as np
from yambopy.tools.string import marquee
from itertools import product
import collections
ha2ev  = 27.211396132
max_exp = 50
min_exp =-100.

def abs2(x):
    return x.real**2 + x.imag**2

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

def histogram_eiv(eiv,weights,emin=-5.0,emax=5.0,step=0.01,sigma=0.05,ctype='lorentzian'):
    """
    Histogram of eigenvalues
    """
    eiv = np.array(eiv)
    #sigma = 0.005
    x = np.arange(emin,emax,step,dtype=np.float32)
    y = np.zeros([len(x)],dtype=np.float32)

    if ctype == 'gaussian':
        from math import sqrt
        c =  1.0/(sigma*sqrt(2))
        a = -1.0/(2*sigma)

    else:
        #lorentzian stuff
        s2 = (.5*sigma)**2
        c = (.5*sigma)

    eiv     = eiv.flatten()
    weights = weights.flatten()

    weights = weights[emin < eiv]
    eiv     = eiv[emin < eiv]
    weights = weights[eiv < emax]
    eiv     = eiv[eiv < emax]

    if ctype == 'gaussian':
        for e,w in zip(eiv,weights):
            x1 = (x-e)**2
            #add gaussian
            y += w*c*np.exp(a*x1)
    else:
        #lorentzian stuff
        for e,w in zip(eiv,weights):
            x1 = (x-e)**2
            #add lorentzian
            y += w*c/(x1+s2)
    return x, y

class YamboElectronsDB():
    """
    Class to read information about the electrons from the ``ns.db1`` produced by yambo

    Arguments:

        ``lattice``: instance of YamboLatticeDB or YamboSaveDB

        ``filename``: netcdf database to read from (default:ns.db1)

    NB: - spin polarized calculations not yet supported
        - spin-polarized eigenvalues are read and expanded for compatibility with DipolesDB
    """
    def __init__(self,lattice,save='SAVE',filename='ns.db1'):
        self.lattice = lattice
        self.filename = '%s/%s'%(save,filename)
        self.efermi = None
        self.readDB()
        self.expandEigenvalues()
        if self.nkpoints != lattice.nkpoints: #sanity check
            raise ValueError("The number of k-points in the lattice database"
                             " %d and electrons database %d is different."%(lattice.nkpoints,self.nkpoints))

    def readDB(self):
        try:
            database = Dataset(self.filename)
        except:
            raise IOError("Error opening file %s in YamboElectronsDB"%self.filename)

        self.eigenvalues_ibz  = database.variables['EIGENVALUES'][0,:]*ha2ev
        self.eigenvalues_ibz_sp_pol  = database.variables['EIGENVALUES'][:,:]*ha2ev
        self.iku_kpoints      = database.variables['K-POINTS'][:].T
        dimensions = database.variables['DIMENSIONS'][:]
        self.nbands      = dimensions[5]
        self.temperature = dimensions[13]
        self.nelectrons  = int(dimensions[14])
        self.nkpoints    = int(dimensions[6])
        self.nbands      = int(dimensions[5])
        self.spin = int(dimensions[12])
        self.spinor_components = int(dimensions[11]) 
        self.time_rev = dimensions[9] 
        database.close()

        #spin degeneracy if 2 components degen 1 else degen 2
        self.spin_degen = [0,2,1][int(self.spin)]
        #number of occupied bands
        # NB: in the spin-polarised case, nbands contains the total number
        #     of bands PER spin polarisation, i.e. half of the total number.
        #     Therefore, nbandsv and nbandsc are also given per
        #     per spin polarisation: this fact is used by DipolesDB
        self.nbandsv = int(self.nelectrons/2)
        self.nbandsc = int(self.nbands-self.nbandsv)
        if (self.spinor_components==2): self.nbandsv = int(self.nelectrons) 
        else:                           self.nbandsv = int(self.nelectrons/2) 
        if self.spin==2:
            self.nbands_tot  = self.nbands*self.spin
            self.nbandsv_tot = int(self.nelectrons/self.spin_degen)
            self.nbandsc_tot = int(self.nbands_tot-self.nbandsv_tot)
        
    def expandEigenvalues(self):
        """
        Expand eigenvalues to the full brillouin zone
        """

        self.eigenvalues = self.eigenvalues_ibz[self.lattice.kpoints_indexes]
        self.eigenvalues_sp_pol = self.eigenvalues_ibz_sp_pol[:,self.lattice.kpoints_indexes]

        self.nkpoints_ibz = len(self.eigenvalues_ibz)
        self.weights_ibz = np.zeros([self.nkpoints_ibz],dtype=np.float32)
        self.nkpoints = len(self.eigenvalues)

        #counter counts the number of occurences of element in a list
        for nk_ibz,inv_weight in list(collections.Counter(self.lattice.kpoints_indexes).items()):
            self.weights_ibz[nk_ibz] = float(inv_weight)/self.nkpoints

        #kpoints weights
        self.weights = np.full((self.nkpoints), 1.0/self.nkpoints,dtype=np.float32)

    def getDOS(self,broad=0.1,emin=-10,emax=10,step=0.01,ctype="lorentzian"):
        """
        Calculate the density of states.
        Should work for metals as well but untested for that case

        It can be used with QP values:
           - Instance yamboQPDB
           - set self.eigenvalues_ibz = yqp.eigenvalues_qp
           - TODO: provide QP broadening QP_broad=yqp.lifetimes
        """
        eigenvalues = self.eigenvalues_ibz
        weights = self.weights_ibz
        nkpoints = self.nkpoints_ibz

        na = np.newaxis
        weights_bands = np.ones(eigenvalues.shape,dtype=np.float32)*weights[:,na]
        energies, self.dos = histogram_eiv(eigenvalues,weights_bands,emin=emin,emax=emax,step=step,sigma=broad,ctype=ctype)

        return energies, self.dos

    def get_transitions(self):
        """
        Calculate transition energies
        """
        eigenvalues = self.eigenvalues_ibz
        nvalence  = self.nbandsv
        nconduction  = self.nbandsc
        nkpoints = self.nkpoints_ibz

        transitions = np.zeros([nkpoints,nvalence*nconduction])
        for k,v,c in product(range(nkpoints),range(nvalence),range(nconduction)):
            vc = v*nconduction+c
            transitions[k,vc] = eigenvalues[k,c+nvalence]-eigenvalues[k,v]
        self.transitions = transitions

        return self.transitions

    def getJDOS(self,broad=0.1,emin=0,emax=10,step=0.01,ctype="lorentzian"):
        """
        Calculate the joint density of states
        """
        transitions = self.get_transitions()
        weights = self.weights_ibz

        na = np.newaxis
        weights_transitions = np.ones(transitions.shape,dtype=np.float32)*weights[:,na]
        energies, self.jdos = histogram_eiv(self.transitions,weights_transitions,emin=emin,emax=emax,step=step,sigma=broad,ctype=ctype)

        return energies, self.jdos

    def setLifetimes(self,broad=0.1):
        """
        Add electronic lifetimes using the DOS
        """
        self.lifetimes_ibz = np.ones(self.eigenvalues_ibz.shape,dtype=np.float32)*broad
        self.lifetimes     = np.ones(self.eigenvalues.shape,dtype=np.float32)*broad

    def setLifetimesDOS(self,broad=0.1,debug=False):
        """
        Approximate the electronic lifetimes using the DOS
        """
        eigenvalues = self.eigenvalues_ibz
        weights = self.weights_ibz
        nkpoints = self.nkpoints_ibz

        #get dos
        emin = np.min(eigenvalues)-broad
        emax = np.max(eigenvalues)+broad
        energies, dos = self.getDOS(emin=emin, emax=emax, step=0.1, broad=broad)

        #normalize dos to broad
        dos = dos/np.max(dos)*broad

        #create a interpolation function to get the lifetimes for all the values
        from scipy.interpolate import interp1d
        f = interp1d(energies, dos, kind='cubic')

        if debug:
            """
            plot the calculated values for the DOS and the interpolated values
            """
            import matplotlib.pyplot as plt
            x = np.arange(emin+d,emax-d,0.001)
            plt.plot(energies,dos,'o')
            plt.plot(x,f(x))
            plt.show()
            exit()

        #add imaginary part to the energies proportional to the DOS
        self.lifetimes_ibz = np.array([ [f(eig) for eig in eigk] for eigk in self.eigenvalues_ibz],dtype=np.float32)
        self.lifetimes     = np.array([ [f(eig) for eig in eigk] for eigk in self.eigenvalues],dtype=np.float32)

    def setFermi(self,fermi,invsmear):
        """
        Set the fermi energy of the system
        """
        self.invsmear = invsmear
        self.efermi = fermi

        #full brillouin zone
        self.eigenvalues     -= self.efermi
        self.occupations = np.zeros([self.nkpoints,self.nbands],dtype=np.float32)
        for nk in range(self.nkpoints):
            self.occupations[nk] = fermi_array(self.eigenvalues[nk,:],0,self.invsmear)

        #for the ibz
        self.eigenvalues_ibz -= self.efermi
        self.occupations_ibz = np.zeros([self.nkpoints_ibz,self.nbands],dtype=np.float32)
        for nk in range(self.nkpoints_ibz):
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

        nbands = int(self.nelectrons/self.spin_degen)
        #top of valence
        top = np.max(eigenvalues[:,nbands])
        #bottom of conduction
        bot = np.max(eigenvalues[:,nbands-1])
        self.efermi = (top+bot)/2
        self.setFermi(self.efermi,broad)

    def energy_gaps(self,GWshift=0.):
        """
        Calculate the enegy of the gap and apply custom rigid shift
        """
        eiv = self.eigenvalues
        nv  = self.nbandsv
        nc  = self.nbandsc

        # First apply shift if there is one
        eiv[:,nv:]+=GWshift

        # Then compute gaps
        homo = np.max(eiv[:,nv-1])
        lumo = np.min(eiv[:,nv])
        Egap = lumo-homo
        for k in eiv:
            if k[nv-1]==homo:
                lumo_dir=k[nv]
        Edir = lumo_dir-homo

        print('DFT Energy gap: %s eV'%Egap)
        print('DFT Direct gap: %s eV'%Edir)
        print('GW shift:       %s eV'%GWshift)

        return np.copy(eiv)

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
            return sum([sum(self.spin_degen*fermi_array(eigenvalues[nk],ef,self.invsmear))*weights[nk] for nk in range(nkpoints)])-self.nelectrons

        fermi = bisect(occupation_minus_ne,min_eival,max_eival)
        if setfermi: self.setFermi(fermi,invsmear)
        return self.efermi

    def __str__(self):
        lines = []; app = lines.append
        app(marquee(self.__class__.__name__))
        app(f"spin polarizations: {self.spin}")
        app(f"spinor components:  {self.spinor_components}")
        app(f"nelectrons: {self.nelectrons}")
        app(f"nbands:     {self.nbands}")
        app(f"nkpoints:   {self.nkpoints}")
        return "\n".join(lines)
