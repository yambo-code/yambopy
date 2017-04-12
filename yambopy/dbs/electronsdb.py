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
            y += c*np.exp(a*x1)
    else:
        #lorentzian stuff
        for e,w in zip(eiv,weights):
            x1 = (x-e)**2
            #add lorentzian
            y += w*c/(x1+s2)
    return x, y

def energy_gaps(eiv,nv,nc,GWshift=0.):
    HOMO = np.max(eiv[:,nv-1])
    LUMO = np.min(eiv[:,nv])
    Egap = LUMO-HOMO
    for k in eiv: 
        if k[nv-1]==HOMO: LUMO_dir=k[nv]
    Edir = LUMO_dir-HOMO
    print('DFT Energy gap: %s eV'%Egap)
    print('DFT Direct gap: %s eV'%Edir)
    eiv[:,nv:]+=GWshift
    print('GW shift:       %s eV'%GWshift)
    return eiv

def lorentzian(x,x0,g):
    height=1./(np.pi*g)
    return height*(g*g)/((x-x0)*(x-x0)+g*g)

def gaussian(x,x0,s):
    height=1./(np.sqrt(2.*np.pi)*s)
    argument=-0.5*((x-x0)/s)**2
    #Avoiding undeflow errors...
    np.place(argument,argument<min_exp,min_exp)
    return height*np.exp(argument)

def IP_eps2(eiv,nv,nc,weights,dipoles,ntot_dip=-1,GWshift=0.,broad=0.01,broadtype='l',nbnds=[-1,-1],Emin=0.,Emax=10.,Esteps=50):
    """Compute independent-particle absorption
    eiv -> eigenvalues_ibz from YamboElectronsDB.
    nv,nc -> nbandsv,nbandsc from YamboElectronsDB (i.e. from QE nscf calculation).
    weights -> weights_ibz from YamboElectronsDB.
    dipoles -> dipoles from YamboDipolesDB.
    ntot_dip -> nbands from YamboDipolesDB (i.e. included in the YAMBO calculation).
    GWshift -> rigid GW shift in eV.
    broad -> broadening of peaks.
    broadtype -> 'l' is lorentzian, 'g' is gaussian.
    nbnds -> number of [valence, conduction] bands included starting from Fermi level. Default means all are included.
    Emin,Emax,Esteps -> frequency range for the plot.
    """
    freq = np.linspace(Emin,Emax,Esteps)
    #Cut bands to the maximum number used for the dipoles
    if ntot_dip>0: 
        eiv = eiv[:,:ntot_dip]
        nc=ntot_dip-nv
    #Print band gap values and apply GW_shift
    eiv = energy_gaps(eiv,nv,nc,GWshift)
    #Check bands to include in the calculation
    if nbnds[0]<0: nbnds[0]=nv
    if nbnds[1]<0: nbnds[1]=nc
    iv = nv-nbnds[0] #first valence
    lc = nv+nbnds[1] #last conduction
    Ecv=[]
    DIP2=[]
    #Transition energies and matrix elements squared
    for c,v in product(range(nv,lc),range(iv,nv)):
        Ecv.append(eiv[:,c]-eiv[:,v])
        DIP2.append(abs2(dipoles[:,0,c-nv,v])+abs2(dipoles[:,1,c-nv,v])+abs2(dipoles[:,2,c-nv,v])) #make more efficient
    Ecv =np.array(Ecv)
    DIP2=np.array(DIP2)
    Ecv =np.swapaxes(Ecv,0,1)
    DIP2=np.swapaxes(DIP2,0,1)
    #Epsilon2
    EPS2= []
    for w in freq:
        EPS2_w=0.
        ind=0
        for c,v in product(range(nv,lc),range(iv,nv)):
            if broadtype=='g': EPS2_w += DIP2[:,ind]*weights*gaussian(w,Ecv[:,ind],broad)
            else:              EPS2_w += DIP2[:,ind]*weights*lorentzian(w,Ecv[:,ind],broad)
            ind=ind+1
        EPS2.append(np.sum(EPS2_w))
    EPS2 = np.array(EPS2)
    #Light polarization directions
    Pdir = 3
    #Spin polarization and kpoint_ibz weights already sum to one
    EPS2 = EPS2/float(Pdir)
    return freq,EPS2

class YamboElectronsDB():
    def __init__(self,lattice,save='SAVE',filename='ns.db1'):
        self.lattice = lattice
        self.filename = '%s/%s'%(save,filename)
        self.efermi = None
        self.readDB()
        if self.nkpoints != self.lattice.nkpoints: #sanity check
            raise ValueError("The number of k-points in the lattice database and electrons database is different.")
        self.expandEigenvalues()
        
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

    def getDOS(self,broad=0.1,emin=-10,emax=10,step=0.01):
        """
        Calculate the density of states.
        Should work for metals as well but untested for that case
        """
        eigenvalues = self.eigenvalues_ibz
        weights = self.weights_ibz
        nkpoints = self.nkpoints_ibz

        na = np.newaxis
        weights_bands = np.ones(eigenvalues.shape,dtype=np.float32)*weights[:,na]
        energies, self.dos = histogram_eiv(eigenvalues,weights_bands,emin=emin,emax=emax,step=step,sigma=broad)

        return energies, self.dos

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
        for nk in xrange(self.nkpoints):
            self.occupations[nk] = fermi_array(self.eigenvalues[nk,:],0,self.invsmear)
        
        #for the ibz
        self.eigenvalues_ibz -= self.efermi
        self.occupations_ibz = np.zeros([self.nkpoints_ibz,self.nbands],dtype=np.float32)
        for nk in xrange(self.nkpoints_ibz):
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

        nbands = self.nelectrons/self.spin_degen
        #top of valence
        top = np.max(eigenvalues[:,nbands])
        #bottom of conduction
        bot = np.max(eigenvalues[:,nbands-1])
        self.efermi = (top+bot)/2
        self.setFermi(self.efermi,broad)
    
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
            return sum([sum(self.spin_degen*fermi_array(eigenvalues[nk],ef,self.invsmear))*weights[nk] for nk in xrange(nkpoints)])-self.nelectrons

        fermi = bisect(occupation_minus_ne,min_eival,max_eival)
        if setfermi: self.setFermi(fermi,invsmear)
        return self.efermi

    def __str__(self):
        s = ""
        s += "spin_degen: %d\n"%self.spin_degen
        s += "nelectrons: %d\n"%self.nelectrons
        s += "nbands:   %d\n"%self.nbands
        s += "nkpoints: %d"%self.nkpoints
        return s
