#
# License-Identifier: GPL
#
# Copyright (C) 2024 The Yambo Team
#
# Authors: FP, HPC, AMS
#
# This file is part of the yambopy project
#
import os
import numpy as np
from netCDF4 import Dataset
from yambopy.tools.string import marquee
from yambopy.tools.funcs import fermi, fermi_array
from yambopy.lattice import car_red, rec_lat, vol_lat
from yambopy.kpoints import expand_kpoints, get_path
from yambopy.plot.spectra import get_spectra
from yambopy.units import ha2ev

class YamboElectronsDB():
    """
    Class to read information about the electrons from the ``ns.db1`` produced by yambo

    Arguments:

        ``filename``: netcdf database to read from (default:ns.db1)

    NB: - spin polarized calculations not yet supported
        - spin-polarized eigenvalues are read and expanded for compatibility with DipolesDB
    """

    def __init__(self,atomic_numbers,car_atomic_positions,eigenvalues_ibz,sym_car,
                 iku_kpoints,nbands,lat,alat,temperature,nelectrons,nkpoints_ibz,spin,time_rev,spinor_components,Expand):
        
        self.atomic_numbers         = atomic_numbers   
        self.car_atomic_positions   = car_atomic_positions
        self.eigenvalues_ibz        = eigenvalues_ibz     
        self.sym_car                = sym_car        
        self.iku_kpoints            = iku_kpoints
        self.nbands                 = nbands       
        self.lat                    = lat             
        self.alat                   = alat            
        self.temperature            = temperature     
        self.nelectrons             = nelectrons
        self.nkpoints_ibz           = nkpoints_ibz       
        self.spin                   = spin            
        self.time_rev               = time_rev
        self.spinor_components      = spinor_components
        self.EXPAND                 = Expand # expand eigenvalues and kpoints to full BZ
        self.expanded               = False  # flag to check if the kpoints were already expanded (starts from False)

        if self.EXPAND: self.expandEigenvalues()

    @classmethod
    def from_db_file(cls,folder='.',filename='ns.db1',Expand=False):
        """
        Read the ns.db1 database
        """
        path_filename = os.path.join(folder,filename)
        try:    database = Dataset(path_filename)
        except: raise IOError("Error opening file {path_filename} in YamboElectronsDB")

        with Dataset(path_filename) as database:
            
            dimensions            = database.variables['DIMENSIONS'][:]
            
            natoms_a = database.variables['N_ATOMS'][:].astype(int).T
            tmp_an = database.variables['atomic_numbers'][:].astype(int)
            tmp_apos = database.variables['ATOM_POS'][:,:]

            flatten = lambda l: [item for sublist in l for item in sublist]
            atomic_numbers = flatten([[tmp_an[n]]*na for n,na in enumerate(natoms_a)])
            atomic_positions = np.vstack([[tmp_apos[n,ia] for ia in range(na)] for n,na in enumerate(natoms_a) ])

            args = dict( atomic_numbers         = atomic_numbers,
                         car_atomic_positions   = atomic_positions,
                         eigenvalues_ibz        = database.variables['EIGENVALUES'][:,:]*ha2ev,
                         sym_car                = np.transpose( database.variables['SYMMETRY'][:], (0,2,1) ), # transpose leaving first axis as symm index
                         iku_kpoints            = database.variables['K-POINTS'][:].T,
                         nbands                 = int(dimensions[5]),
                         lat                    = database.variables['LATTICE_VECTORS'][:].T,
                         alat                   = database.variables['LATTICE_PARAMETER'][:].T,
                         temperature            = dimensions[13],
                         nelectrons             = int(dimensions[14]),
                         nkpoints_ibz           = int(dimensions[6]),
                         spin                   = int(dimensions[12]),
                         time_rev               = dimensions[9],
                         spinor_components      = int(dimensions[11]),
                         Expand                 = Expand 
                         )

        return cls(**args)

    @property
    def red_atomic_positions(self):
        return car_red(self.car_atomic_positions,self.lat)

    @property
    def spin_degen(self):
        """spin degeneracy if 2 components degen 1 else degen 2"""
        return [0,2,1][int(self.spin)]

    @property
    def min_eival(self):
        return np.min(self.eigenvalues) 
    
    @property
    def max_eival(self):
        return np.max(self.eigenvalues)

    @property
    def car_kpoints(self):
        """convert form internal yambo units to cartesian lattice units"""
        return np.array([ k/self.alat for k in self.iku_kpoints ])

    @property
    def red_kpoints(self):
        """convert from cartesian coordinates to reduced coordinates"""
        if not hasattr(self,"_red_kpoints"):
            self._red_kpoints = car_red(self.car_kpoints,self.rlat)
        return self._red_kpoints

    @property
    def rlat(self):
        """caclulate the reciprocal lattice"""
        return rec_lat(self.lat)

    @property
    def rlat_vol(self):
        return (2*np.pi)**3 * vol_lat(self.rlat)

    @property
    def lat_vol(self):
        return vol_lat(self.lat)

    @property
    def natoms(self):
        return len(self.atomic_positions)

    @property
    def nbandsv(self):
        """
        number of occupied bands
        """
        if (self.spinor_components==2): return int(self.nelectrons) 
        else:                           return int(self.nelectrons/2)

    @property
    def nbandsc(self):
        return int(self.nbands-self.nbandsv)

    @property
    def nbands_tot(self):
        """
        NB: in the spin-polarised case, nbands contains the total number
            of bands PER spin polarisation, i.e. half of the total number.
            Therefore, nbandsv and nbandsc are also given per
            spin polarisation: this fact is used by DipolesDB    
        """
        if self.spin==2: return self.nbands*self.spin
        else           : return self.nbands
    
    @property
    def nbandsv_tot(self):
        if self.spin==2: return int(self.nelectrons/self.spin_degen)
        else:            return self.nbandsv
    
    @property
    def nbandsc_tot(self):
        if self.spin==2: int(self.nbands_tot-self.nbandsv_tot)
        else:            return self.nbandsc

    @property
    def time_rev_list(self):
        """get a list of symmetries with time reversal"""
        time_rev_list = [False]*self.nsym
        for i in range(self.nsym):
            time_rev_list[i] = ( i >= self.nsym/(self.time_rev+1) )
        return time_rev_list

    @property
    def sym_rlu(self):
        """convert cartesian transformations to reduced transformations """
        sym_rlu = np.zeros([self.nsym,3,3])
        for n,s in enumerate(self.sym_car):
            a = np.dot(s.T,np.linalg.inv(self.rlat))
            sym_rlu[n] = np.dot(np.linalg.inv(self.lat.T),a)
        return sym_rlu

    @property
    def nsym(self):
        return len(self.sym_car)

    @property
    def sym_red(self):
        """Convert cartesian transformations to reduced transformations"""
        if not hasattr(self,"_sym_red"):
            sym_red = np.zeros([self.nsym,3,3],dtype=int)
            for n,s in enumerate(self.sym_car):
                sym_red[n] = np.round(np.dot(np.dot(self.lat,s.T),np.linalg.inv(self.lat)))
            self._sym_red = sym_red
        return self._sym_red

    @property
    def sym_rec_red(self):
        """Convert reduced transformations to reduced reciprocal transformations"""
        if not hasattr(self,"_sym_rec_red"):
            sym_rec_red = np.zeros([self.nsym,3,3],dtype=int)
            for n,s in enumerate(self.sym_red):
                sym_rec_red[n] = np.linalg.inv(s).T
            self._sym_rec_red = sym_rec_red
        return self._sym_rec_red

    @property
    def sym_rec(self):
        """Convert cartesian transformations to reciprocal transformations"""
        sym_rec = np.zeros([self.nsym,3,3])
        for n,s in enumerate(self.sym_car):
            sym_rec[n] = np.linalg.inv(s).T
        return sym_rec

    @property
    def efermi(self):
        if not hasattr(self,"_efermi") and self.Expand:
            # AMS: To be tested
            self._efermi = self.GetFermi()
        return self._efermi

    def GetFermi(self,inv_smear=0.001,verbose=0,setfermi=True):
        """ Determine the fermi energy
        """
        from scipy.optimize import bisect

        kpts, nks, nss = self.expand_kpoints()

        def fermi_array(e_array,ef):
            """ Fermi dirac function for an array
            """
            e_array = (e_array-ef)/inv_smear
            return [ fermi(e) for e in e_array]

        def occupation_minus_ne(ef):
            """ The total occupation minus the total number of electrons
            """
            if self.spinor == 1:
               return sum([sum(self.spin_degen*fermi_array(self.eigenvalues[0,nk],ef))*self.weights[nk] for nk in range(self.nkpoints)])-self.electrons
            elif self.spinor == 2:
               sum_up = sum([sum(self.spin_degen*fermi_array(self.eigenvalues[0,nk],ef))*self.weights[nk] for nk in range(self.nkpoints)]) 
               sum_dw = sum([sum(self.spin_degen*fermi_array(self.eigenvalues[1,nk],ef))*self.weights[nk] for nk in range(self.nkpoints)]) 
               return sum_up + sum_dw -self.electrons
     
        efermi = bisect(occupation_minus_ne,self.min_eival,self.max_eival)

        if verbose: print("fermi: %lf eV"%efermi)

        if setfermi: return self.setFermi(efermi, inv_smear)
        else:        return efermi

    def setFermi(self,fermi,invsmear):
        """
        Shift bands and get occupations
        """
        self.invsmear = invsmear
        self.efermi = fermi

        #full brillouin zone
        self.eigenvalues     -= self.efermi
        self.occupations = np.zeros([self.spin,self.nkpoints,self.nbands],dtype=np.float32)
        for nspin in range(self.spin):
            for nk in range(self.nkpoints):
                self.occupations[nspin,nk] = fermi_array(self.eigenvalues[nspin,nk,:self.nbands],0)

        #for the ibz
        self.eigenvalues_ibz -= self.efermi
        self.occupations_ibz = np.zeros([self.spin,self.nkpoints_ibz,self.nbands],dtype=np.float32)
        for nspin in range(self.spin):
            for nk in range(self.nkpoints_ibz):
                self.occupations_ibz[nk] = fermi_array(self.eigenvalues_ibz[nspin,nk,:],0,self.invsmear)

        return self.efermi

    def setFermiFixed(self,broad=1e-5):
        """
        Set Fermi level using fixed occupations method
        Useful for semi-conductors
        """
        eigenvalues = self.eigenvalues_ibz

        #top of valence
        top = np.max(eigenvalues[:,:,self.nbandsv-1])
        #bottom of conduction
        bot = np.max(eigenvalues[:,:,self.nbandsv])
        efermi = (top+bot)/2.
        self.setFermi(efermi,broad)

    def expandEigenvalues(self):
        """
        Expand eigenvalues to the full brillouin zone
        """
        
        self.expand_kpoints()
        self.eigenvalues = self.eigenvalues_ibz[:,self.kpoints_indexes]
        self.nkpoints = len(self.eigenvalues)

    def expand_kpoints(self,verbose=1,atol=1.e-4):
        """ 
        Wrapper for expand_kpoints in kpoints module

        Take a list of qpoints and symmetry operations and return the full brillouin zone
        with the corresponding index in the irreducible brillouin zone
        """

        #check if the kpoints were already exapnded
        if self.expanded == True: return self.kpoints_full, self.kpoints_indexes, self.symmetry_indexes

        weights, kpoints_indexes, symmetry_indexes, kpoints_full = expand_kpoints(self.car_kpoints,self.sym_car,self.rlat,atol=atol)

        #set the variables
        self.expanded = True
        self.kpoints_full     = kpoints_full
        self.kpoints_indexes  = kpoints_indexes
        self.symmetry_indexes = symmetry_indexes
        self.nkpoints         = len(kpoints_full)
        self.weights_ibz      = weights
        self.weights          = np.full((self.nkpoints), 1.0/self.nkpoints,dtype=np.float32)

        if verbose: print("%d kpoints expanded to %d"%(len(self.car_kpoints),len(kpoints_full)))

        return self.kpoints_full, self.kpoints_indexes, self.symmetry_indexes

    def energy_gaps(self,eigenvalues=None,GWshift=0.):
        """
        Calculate the energy of the gap and apply custom rigid shift

        The eigenvalues array [nk,nb] can be given as an argument (useful for QP energies or spin-polarized case),
        otherwise the eigenvalues are taken from the object (spin up channel).
        """
        if eigenvalues is None: eiv = self.eigenvalues_ibz[0]
        else:                   eiv = eigenvalues
        nv  = self.nbandsv

        # First apply shift if there is one
        eiv[:,nv:]+=GWshift

        # Then compute gaps
        Egap = np.min(eiv[:,nv]) - np.max(eiv[:,nv-1])
        Edir = np.min(eiv[:,nv]  -        eiv[:,nv-1])

        print('DFT Energy gap: %s eV'%Egap)
        print('DFT Direct gap: %s eV'%Edir)
        print('GW shift:       %s eV'%GWshift)

        return eiv

    ##
    ## Simple band plot
    ##
    def plot_bs_ax(self,ax,path,bandmin=None,bandmax=None,add_indexes=False,**kwargs):
        """
        Plot this bandstructure on Matpltolib ax
        """
        bands_kpoints, bands_indexes, path_car = get_path(self.car_kpoints,self.rlat,path)
        bands_highsym_qpts = path_car.kpoints
        self.get_fermi()

        #calculate distances
        bands_distances = [0]
        distance = 0
        for nk in range(1,len(bands_kpoints)):
            distance += np.linalg.norm(bands_kpoints[nk]-bands_kpoints[nk-1])
            bands_distances.append(distance)

        #plot high-symmetry qpoints
        distance = 0
        bands_highsym_qpts_distances = [0]
        for nk in range(1,len(bands_highsym_qpts)):
            ax.axvline(distance,color='k')
            distance += np.linalg.norm(bands_highsym_qpts[nk]-bands_highsym_qpts[nk-1])
            bands_highsym_qpts_distances.append(distance)
        ax.axvline(distance,color='k')

        #plot bands
        if self.spin == 2:
           color = kwargs.pop('c','red')
           ax.plot(bands_distances,self.eigenvalues[0,bands_indexes,bandmin:bandmax],c=color,**kwargs)
           color = kwargs.pop('c','blue')
           ax.plot(bands_distances,self.eigenvalues[1,bands_indexes,bandmin:bandmax],c=color,**kwargs)
        else:
            color = kwargs.pop('c','red')
            ax.plot(bands_distances,self.eigenvalues[0,bands_indexes,bandmin:bandmax],c=color,**kwargs)
        
        ax.set_xlim(0,max(bands_distances))

        if add_indexes:
            ax.set_xticks(bands_distances)
            ax.set_xticklabels(np.array(bands_indexes)+1)
            for d in bands_distances:
                ax.axvline(d,color='k',alpha=0.5)

        return ax

    ##
    ## DOS and JDOS
    ##
    def get_transitions(self,eigenvalues=None,nvalence=None,nconduction=None):
        """
        Calculate transition energies
        """
        if eigenvalues is None: eigenvalues = self.eigenvalues_ibz[0] #assume non spin polarized if not given
        if nvalence is None:    nvalence    = self.nbandsv
        if nconduction is None: nconduction = self.nbandsc

        eigs_valence    = eigenvalues[:, :nvalence] 
        eigs_conduction = eigenvalues[:, nvalence:nvalence + nconduction]
        transitions = eigs_conduction[:, None, :] - eigs_valence[:, :, None]
        transitions = transitions.reshape(self.nkpoints_ibz, nvalence * nconduction)

        return transitions

    def getDOS(self,eigenvalues=None,broad=0.1,emin=-10,emax=10,estep=0.01):
        """
        Retrieve the Density of States (DOS)
        (Should work for metals as well but untested)
        """
        if eigenvalues is None: eigenvalues = self.eigenvalues_ibz[0] # assume non spin polarized if not given
        if not self.EXPAND: self.expandEigenvalues()                  # we need this to get the weights
        weights = self.weights_ibz[:self.nkpoints_ibz]
        
        w, dos = get_spectra(eigenvalues,weights=weights,broadening=broad,emin=emin,emax=emax,estep=estep)

        return w, dos

    def getJDOS(self,eigenvalues=None,nconduction=None,broad=0.1,emin=0,emax=10,estep=0.01):
        """
        Calculate the joint density of states
        """
        transitions = self.get_transitions(eigenvalues=eigenvalues,nconduction=nconduction)
        if not self.EXPAND: self.expandEigenvalues() # we need this to get the weights
        weights = self.weights_ibz[:self.nkpoints_ibz]

        w, jdos = get_spectra(transitions,weights=weights,broadening=broad,emin=emin,emax=emax,estep=estep)

        return w, jdos
    
    def setLifetimes(self,broad=0.1):
        """
        Add electronic lifetimes using the DOS
        """
        self.lifetimes_ibz = np.ones(self.eigenvalues_ibz.shape,dtype=np.float32)*broad
        if self.Expand: self.lifetimes = np.ones(self.eigenvalues.shape,dtype=np.float32)*broad

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
            d=0.001
            x = np.arange(emin-d,emax+d,d)
            plt.plot(energies,dos,'o')
            plt.plot(x,f(x))
            plt.show()
            exit()

        #add imaginary part to the energies proportional to the DOS
        self.lifetimes_ibz = np.array([ [f(eig) for eig in eigk] for eigk in self.eigenvalues_ibz],dtype=np.float32)
        if self.Expand: self.lifetimes = np.array([ [f(eig) for eig in eigk] for eigk in self.eigenvalues],dtype=np.float32)
    
    def __str__(self):
        lines = []; app = lines.append
        app(marquee(self.__class__.__name__))
        app(f"spin polarizations: {self.spin}")
        app(f"spinor components:  {self.spinor_components}")
        app(f"nelectrons: {self.nelectrons}")
        app(f"nbands:     {self.nbands}")
        app(f"nkpoints (IBZ):   {self.nkpoints_ibz}")
        return "\n".join(lines)
