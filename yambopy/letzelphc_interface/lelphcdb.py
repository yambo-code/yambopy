import numpy as np
from netCDF4 import Dataset
from yambopy.tools.string import marquee
from yambopy.units import ha2ev

class LetzElphElectronPhononDB():
    """
    Python class to read the electron-phonon matrix elements from LetzElPhC.

    About LetzElPhC: https://github.com/muralidhar-nalabothula/LetzElPhC
    
    By default it reads the full database g(k,q,m,s,b1,b2) including phonon energies.
    
    - Input: path of ndb.elph
    - Input: read_all (default True), read ph. eigenvectors and el-ph matrix elements
    - Input: div_by_energies (default True), divide el-ph mat. el. by sqrt(2* ph. energies)

    - Usage and main variables: 
    
      :: lph = LetzElphElectronPhononDB(path_ndb_elph)
    
      :: lph.kpoints         #kpoints in cryst. coords. (BZ)
      :: lph.qpoints         #qpoints in crist. coords. (BZ)
      :: lph.ph_energies     #Phonon energies (eV)      
      :: lph.ph_eigenvectors #Phonon modes
      :: lph.gkkp            #El-ph matrix elements (by default normalised with ph. energies):
      :: lph.gkkp_sq         #Couplings (square)

       
   
    Formats:
    - modes[iq][il][iat][ix]
    - gkkp[iq][ik][il][is][ib1][ib2]              
    """

    def __init__(self,filename,read_all=True,div_by_energies=True):

        # Open database
        try: database = Dataset(filename)
        except: raise FileNotFoundError("error opening %s in LetzElphElectronPhononDB"%filename)

        # Read DB dimensions
        self.nb1 = database.dimensions['initial_band'].size
        self.nb2 = database.dimensions['final_band_PH_abs'].size
        self.nm = database.dimensions['nmodes'].size
        self.nat = database.dimensions['atom'].size
        self.nk = database.dimensions['nk'].size
        self.nq = database.dimensions['nq'].size
        self.ns = database.dimensions['nspin'].size
        self.nsym = database.dimensions['nsym_ph'].size

        # Read DB
        self.kpoints = database.variables['kpoints'][:]
        self.qpoints = database.variables['qpoints'][:]
        self.bands   = database.variables['bands'][:]
        
        self.ph_energies = database.variables['FREQ'][:]*(ha2ev/2.) # Energy units are in Rydberg
        self.check_energies()

        if read_all: 
 
            self.read_eigenmodes(database)
            self.read_elph(database,scale_g_with_ph_energies=div_by_energies)

        for var in database.variables.values():
            if var.name=='elph_mat': self.ncfloat_type=var.dtype 

        database.close()

    def check_energies(self):
        """
        Inform the user about unexpected negative frequencies and set them to positive
        """
        indices = np.where(self.ph_energies < 0.)
        warn = False
        for Q in indices[0]:
            for M in indices[1]:
                if Q==0 and M in [0,1,2]: 
                    self.ph_energies[Q,M]=0.
                else:
                    warn = True
                    self.ph_energies[Q,M]=np.abs(self.ph_energies[Q,M])
        if warn: print("[WARNING] Found unexpected negative phonon energies, be careful!")

    def read_eigenmodes(self,database):
        """
        Read phonon eigenmodes
        """

        self.ph_eigenvectors = np.zeros([self.nq,self.nm,self.nat,3],dtype=np.complex64)
        #eivs_tmp[qpt][mode][atom][coord][cmplx]
        eivs_tmp = database.variables['POLARIZATION_VECTORS'][:]
        self.ph_eigenvectors = eivs_tmp[:,:,:,:,0] + 1j*eivs_tmp[:,:,:,:,1]

    def read_elph(self,database,scale_g_with_ph_energies=True):
        """
        Read electron-phonon matrix elements
        
        - If scale_g_with_ph_energies they are divided by sqrt(2*ph_E)
        """    
        gkkp_full = np.zeros([self.nq,self.nk,self.nm,self.ns,self.nb1,self.nb2],dtype=np.complex64)
        gkkp_tmp  = database.variables['elph_mat'][:]
        gkkp_full = gkkp_tmp[:,:,:,:,:,:,0]+1j*gkkp_tmp[:,:,:,:,:,:,1]
        
        # Check integrity of elph values
        if np.isnan(gkkp_full).any(): print('[WARNING] NaN values detected in elph database.')
        
        # Scaling with phonon energies
        if scale_g_with_ph_energies: gkkp_full = self.scale_g(gkkp_full) 

        self.gkkp = gkkp_full
        self.gkkp_sq = np.abs(gkkp_full)**2. 

    def scale_g(self,dvscf):
        """
        Normalise matrix elements by the phonon energy as: 
       
        g_qnu = dvscf_qnu/sqrt(2*w_qnu)
        """
        
        g = np.zeros([self.nq,self.nk,self.nm,self.ns,self.nb1,self.nb2],dtype=np.complex64)
        for iq in range(self.nq):
            for inu in range(self.nm): 
                if iq==0 and inu in [0,1,2]: 
                    g[iq,:,inu,:,:,:] = 0. # Remove acoustic branches
                else:
                    ph_E = self.ph_energies[iq,inu]/(ha2ev/2.) # Put back the energies in Rydberg units
                    g[iq,:,inu,:,:,:] = dvscf[iq,:,inu,:,:,:]/np.sqrt(2.*ph_E)
        return dvscf

    def __str__(self,verbose=False):

        lines = []; app = lines.append
        app(marquee(self.__class__.__name__))
            
        app('nqpoints: %d'%self.nq)
        app('nkpoints: %d'%self.nk)
        app('nmodes: %d'%self.nm)
        app('natoms: %d'%self.nat)
        app('nbands: %d %d'%(self.nb1,self.nb2))
 
        if verbose:

            if hasattr(self, 'ph_eigenvectors'):                 
                app('-----------------------------------')
                for iq in range(self.nq):
                    app('nqpoint %d'%iq)
                    for n,mode in enumerate(self.ph_eigenvectors[iq]):
                        app('mode %d freq: %lf meV'%(n,self.ph_energies[iq,n]*1000.))
                        for a in range(self.nat):
                            app(("%12.8lf "*3)%tuple(mode[a].real))
                    app('-----------------------------------')
            else:
                app('-----------------------------------')
                app('Eigenvectors and matrix elements not loaded')
                app('-----------------------------------')

        return "\n".join(lines)