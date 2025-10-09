import numpy as np
from netCDF4 import Dataset
from yambopy.tools.string import marquee
from yambopy.units import ha2ev
from yambopy.kpoints import build_ktree, find_kpt

class LetzElphElectronPhononDB():
    """
    Python class to read the electron-phonon matrix elements from LetzElPhC.

    About LetzElPhC: https://gitlab.com/lumen-code/LetzElPhC
    
    By default it reads the full database g(q,k,m,s,b1,b2) including phonon energies.
    
    - Input: path of ndb.elph
    - Input: read_all (default True), read ph. eigenvectors and el-ph matrix elements
    - Input: div_by_energies (default True), divide el-ph mat. el. by sqrt(2* ph. energies)

    - Usage and main variables: 
    
      :: lph = LetzElphElectronPhononDB(path_ndb_elph)
    
      :: lph.kpoints         #kpoints in cryst. coords. (BZ)
      :: lph.qpoints         #qpoints in crist. coords. (BZ)
      :: lph.ph_energies     #Phonon energies (eV), energies in LetzElPhCode [Ry]      
      :: lph.ph_eigenvectors #Phonon modes
      :: lph.gkkp            #El-ph matrix elements (by default normalised with ph. energies) [!!!! RYDBERG UNITS !!!!]:
      :: lph.gkkp_sq         #Couplings (square)

    Formats:
    - modes[iq][il][iat][ix]
    - gkkp[iq][ik][il][is][ib1][ib2]              
    """

    def __init__(self,filename,read_all=True,div_by_energies=True,verbose=False):

        # Open database
        try: database = Dataset(filename)
        except: raise FileNotFoundError("error opening %s in LetzElphElectronPhononDB"%filename)

        self.filename = filename
        # Read DB dimensions
        self.nb1 = database.dimensions['initial_band'].size
        self.nb2 = database.dimensions['final_band_PH_abs'].size
        self.nm = database.dimensions['nmodes'].size
        self.nat = database.dimensions['atom'].size
        self.nk = database.dimensions['nk'].size
        self.nq = database.dimensions['nq'].size
        self.ns = database.dimensions['nspin'].size
        self.nsym = database.dimensions['nsym_ph'].size
        self.div_by_energies = div_by_energies # if true, the elph store are normalized with 1/(2*w_ph)
        #
        #
        conv = database['convention'][...].data
        if isinstance(conv, np.ndarray):
            if conv.dtype.kind == 'S':  # Byte strings (C chars)
                conv = conv.tobytes().decode('utf-8').strip()
            else:
                conv = str(conv)  # Fallback for non-string arrays
        elif isinstance(conv, bytes):
            conv = conv.decode('utf-8').strip()
        else:
            conv = str(conv).strip()
        conv = conv.strip().replace('\0', '')
        #
        #
        if conv == 'standard':
            print("Convention used in Letzelphc : k -> k+q (standard)")
        else:
            print("Convention used in Letzelphc : k-q -> k (yambo)")
        self.convention = conv
        #
        # Read DB
        #For developes: Do not use database.variables['x'][:]
        #prefer to use : databaase.variables['x'][...].data to preserve precision
        self.kpoints = database.variables['kpoints'][...].data
        self.qpoints = database.variables['qpoints'][...].data
        self.bands   = database.variables['bands'][...].data
        self.ktree   = build_ktree(self.kpoints)
        self.qtree   = build_ktree(self.qpoints)
        self.kmap = database.variables['kmap'][...].data
        
        self.ph_energies = database.variables['FREQ'][...].data*(ha2ev/2.) # From [Ry] to [eV]
        self.check_energies()

        if read_all: 
 
            self.read_eigenmodes(database)
            self.read_elph(database,scale_g_with_ph_energies=div_by_energies)

        for var in database.variables.values():
            if var.name=='elph_mat': self.ncfloat_type=var.dtype 

        database.close()
        
        self.verbose = verbose

    def check_energies(self):
        """
        Inform the user about unexpected negative frequencies and set them to positive
        """
        indices = np.where(self.ph_energies < 0.)
        warn = False
        for Q in indices[0]:
            for M in indices[1]:
                if Q==0 and M in [0,1,2]: 
                    print('Acoustic modes have been set to zero')
                    self.ph_energies[Q,M]= 0 #np.abs(self.ph_energies[Q,M])#self.ph_energies[Q,M]=0.
                else:
                    warn = True
                    self.ph_energies[Q,M]=np.abs(self.ph_energies[Q,M])
        if warn: print("[WARNING] Found unexpected negative phonon energies, be careful!")

    def read_eigenmodes(self,database):
        """
        Read phonon eigenmodes
        """

        eivs_tmp = database.variables['POLARIZATION_VECTORS'][...].data
        self.ph_eigenvectors = np.zeros([self.nq,self.nm,self.nat,3],dtype=eivs_tmp.dtype)
        self.ph_eigenvectors = eivs_tmp[:,:,:,:,0] + 1j*eivs_tmp[:,:,:,:,1]

    def read_elph(self,database,scale_g_with_ph_energies=True):
        """
        Read electron-phonon matrix elements
        
        - If scale_g_with_ph_energies they are divided by sqrt(2*ph_E)
        """    
        gkkp_tmp  = database.variables['elph_mat'][...].data
        gkkp_full = np.zeros([self.nq,self.nk,self.nm,self.ns,self.nb1,self.nb2],dtype=gkkp_tmp.dtype)
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
        
        g = np.zeros([self.nq,self.nk,self.nm,self.ns,self.nb1,self.nb2],dtype=self.ph_eigenvectors.dtype)
        for iq in range(self.nq):
            for inu in range(self.nm): 
                if iq==0 and inu in [0,1,2]: 
                    g[iq,:,inu,:,:,:] = 0. # Remove acoustic branches
                else:
                    ph_E = self.ph_energies[iq,inu]/(ha2ev/2.) # Put back the energies in Rydberg units
                    g[iq,:,inu,:,:,:] = dvscf[iq,:,inu,:,:,:]/np.sqrt(2.*ph_E)
        return g


    def read_iq(self,iq, bands_range=[], database=None, convention='yambo'):
        """
        !!!! This method works in Ry !!!!
        Reads the electron-phonon matrix elements and phonon eigenvectors for a specific q-point index.

        If the data is already loaded in memory, it returns the corresponding array slice. Otherwise,
        it reads from the database without storing the data in memory.
        
        This function reads data for a single q-point instead of the entire dataset, which is useful
        for handling large databases efficiently.

        Parameters
        ----------
        iq : int
            Index of the q-point.
        bands_range : list, optional
            Specifies the range of bands to read. The start index follows Python indexing (starting from 0),
            and the end index is excluded. If not provided, it defaults to the minimum and maximum bands available.
        database : Dataset, optional
            If provided, the function will use this open dataset instead of opening the file again.
        convention : str, optional
            Defines the convention used for electron-phonon matrix elements.
            - 'yambo': Outputs \<k|dV|k-q>.
            - Any other value: Outputs \<k+q|dV|k>.

        Returns
        -------
        tuple
            A tuple containing:
            - ph_eigenvectors : ndarray
                The phonon eigenvectors.
            - ph_elph_me : ndarray
                The electron-phonon matrix elements with the specified convention [QE convention Ry].
        """
        #
        if len(bands_range) == 0:
            bands_range = [min(self.bands)-1,max(self.bands)]
        min_bnd = min(bands_range)
        max_bnd = max(bands_range)
        nbnds = max_bnd - min_bnd
        assert (min_bnd >= min(self.bands)-1)
        assert (max_bnd <= max(self.bands))
        start_bnd_idx = 1+min_bnd - min(self.bands)
        end_bnd = start_bnd_idx + nbnds
        
        if hasattr(self, 'ph_eigenvectors'):
            ph_eigs = self.ph_eigenvectors[iq]
            eph_mat = self.gkkp[iq, :, :, :, start_bnd_idx:end_bnd, start_bnd_idx:end_bnd ]
        else :
            ## else we load from the file
            close_file = False
            if not database :
                close_file = True
                database = Dataset(self.filename,'r')
            eph_mat = database['elph_mat'][iq, :, :, :, start_bnd_idx:end_bnd, start_bnd_idx:end_bnd, :].data
            # ( nk, nm, nspin, initial bnd, final bnd)
            ph_eigs = database['POLARIZATION_VECTORS'][iq,...].data
            eph_mat = eph_mat[...,0] + 1j*eph_mat[...,1]
            ph_eigs = ph_eigs[...,0] + 1j*ph_eigs[...,1]
            if self.div_by_energies:
                ph_freq_iq = np.sqrt(2.0*np.abs(self.ph_energies[iq]/(ha2ev/2.)))
                if iq > 0:
                    ph_freq_iq = 1.0/ph_freq_iq
                    eph_mat[:, 3:] *= ph_freq_iq[None,:,None,None,None]
            if close_file :database.close()
        return [ph_eigs, self.change_convention(self.qpoints[iq],eph_mat, convention).astype(eph_mat.dtype)]
        ## output elph matrix elements unit (Ry if div_by_energies else Ry^1.5)
        # ( nk, nm, nspin, initial bnd, final bnd)
    def change_convention(self, qpt, elph_iq, convention='yambo'):
        """
        Adjusts the convention of the electron-phonon matrix elements.
    
        Parameters
        ----------
        qpt : ndarray
            The q-point in crystal coordinates.
        elph_iq : ndarray
            The electron-phonon matrix elements.
        convention : str, optional
            Defines the output format:
            - 'yambo': Outputs \<k|dV|k-q>.
            - Any other value: Outputs \<k+q|dV|k>.

        Returns
        -------
        ndarray
            The electron-phonon matrix elements in the desired convention. The returned array is a view, not a copy.
        """
        if convention.strip() != 'yambo': convention = 'standard'
        if self.convention == convention: return elph_iq
        if convention == 'standard': factor = 1.0
        else: factor = -1.0
        idx_q = find_kpt(self.ktree, factor*qpt[None, :] + self.kpoints)
        return elph_iq[idx_q,:, ...]

    def __str__(self):

        lines = []; app = lines.append
        app(marquee(self.__class__.__name__))
            
        app('nqpoints: %d'%self.nq)
        app('nkpoints: %d'%self.nk)
        app('nmodes: %d'%self.nm)
        app('natoms: %d'%self.nat)
        app('nbands: %d %d'%(self.nb1,self.nb2))
        app('convention: %s'%self.convention)

        if self.verbose:

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
