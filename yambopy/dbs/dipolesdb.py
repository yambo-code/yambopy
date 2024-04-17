# Authors: HPCM, FP
#
# This file is part of the yambopy project
#
from yambopy.units import *
from yambopy import *
from math import sqrt
from time import time
import matplotlib.pyplot as plt
from yambopy.tools.string import marquee
from yambopy.tools.funcs import abs2,lorentzian, gaussian
from yambopy.plot.plotting import add_fig_kwargs,BZ_Wigner_Seitz,shifted_grids_2D

class YamboDipolesDB():
    """
    Class to read the dipoles databases from the ``ndb.dipoles`` files
    
    Can be used to plot, for example, the imaginary part of the dielectric
    function which corresponds to the optical absorption, or directly the matrix elements in kspace.

    Dipole matrix elements <ck|vec{r}|vk> are stored in self.dipoles with indices [k,r_i,c,v]. If the calculation is spin-polarised (nk->nks), then they are stored with indices [s,k,r_i,c,v]
    """
    def __init__(self,lattice,save='SAVE',filename='ndb.dip_iR_and_P',dip_type='iR',field_dir=[1,0,0],field_dir3=[0,0,1]):
        self.lattice = lattice
        self.filename = "%s/%s"%(save,filename)
        
        #read dipoles
        try:
            database = Dataset(self.filename, 'r')
        except:
            raise IOError("Error opening %s in YamboDipolesDB"%self.filename)
            
        self.nq_ibz, self.nq_bz, self.nk_ibz, self.nk_bz = database.variables['HEAD_R_LATT'][:].astype(int)
        self.spin = database.variables['SPIN_VARS'][0].astype(int)

        # indexv is the maximum partially occupied band
        # indexc is the minimum partially empty band
        self.min_band, self.max_band, self.indexv, self.indexc = database.variables['PARS'][:4].astype(int)
        database.close()

        # determine the number of bands
        self.nbands  = self.max_band-self.min_band+1
        self.nbandsv = self.indexv-self.min_band+1
        self.nbandsc = self.max_band-self.indexc+1
        self.indexv = self.indexv-1
        self.indexc = self.indexc-1
        self.index_firstv = self.min_band-1
        self.open_shell = False
        d_n_el = self.nbandsv + self.nbandsc - self.nbands
        if d_n_el != 0:
            # We assume the excess electron(s)/hole(s) are in the spin=0 channel
            print("[WARNING] You may be considering an open-shell system. Careful: optical absorption UNTESTED.")
            self.open_shell = True
            self.n_exc_el = d_n_el # this can be negative for excess holes
            # Spin-dependent band indices: these are taken from PARS so be careful
            self.indexv_os  = [self.indexv, self.indexv-d_n_el]
            self.indexc_os  = [self.indexc+d_n_el, self.indexc]
            self.nbandsv_os = [self.nbandsv, self.nbandsv-d_n_el ] 
            self.nbandsc_os = [self.nbandsc-d_n_el, self.nbandsc ]

        #read the database
        self.dipoles = self.readDB(dip_type)

        #expand the dipoles to the full brillouin zone
        if self.spin==1: self.expandDipoles(self.dipoles)
        if self.spin==2:
            dip_up, dip_dn = self.dipoles[0], self.dipoles[1]
            exp_dip      = self.expandDipoles
            self.dipoles = np.stack((exp_dip(dip_up,spin=0)[0], exp_dip(dip_dn,spin=1)[0]),axis=0) 

    def normalize(self,electrons):
        """ 
        Use the electrons to normalize the dipole matrix elements
        """
        # We take the eivs with the added spin dimensions even in non-spin pol case
        eiv = electrons.eigenvalues_sp_pol
        nkpoints, nbands = eiv[0].shape

        dipoles = self.dipoles
        if self.spin==1: dipoles = np.expand_dims(dipoles,axis=0)
        print(dipoles.shape)
        for nk in range(nkpoints):
            for ns in range(self.spin):

                eiv_sk = eiv[ns,nk]
            
                #create eigenvalues differences arrays
                norm = np.array([ [ec-ev for ev in eiv_sk] for ec in eiv_sk ])
            
                #normalize
                for i,j in product(list(range(nbands)),repeat=2):

                    if norm[i,j] == 0: dipoles[ns,nk,:,i,j] = 0.
                    else: dipoles[ns,nk,:,i,j] = dipoles[ns,nk,:,i,j]/norm[i,j]

        if self.spin==1: dipoles=np.squeeze(dipoles)
        self.dipoles = dipoles

    def readDB(self,dip_type):
        """
        The dipole matrix has the following indexes:
        [nspin, nkpoints, cartesian directions, nbands conduction, nbands valence]
        """
        #check if output is in the old format
        fragmentname = "%s_fragment_1"%(self.filename)
        if os.path.isfile(fragmentname): return self.readDB_oldformat(dip_type)

        self.dip_type = dip_type
        if self.spin==1: 
            dipoles = np.zeros([self.nk_ibz,3,self.nbandsc,self.nbandsv],dtype=np.complex64)
        if self.spin==2:
            dipoles = np.zeros([self.spin,self.nk_ibz,3,self.nbandsc,self.nbandsv],dtype=np.complex64)
        
        database = Dataset(self.filename)
        dip = database.variables['DIP_%s'%(dip_type)]
        if self.spin==1:
            dip = np.squeeze(dip)
            dip = (dip[:,:,:,:,0]+1j*dip[:,:,:,:,1]) # Read as nk,nv,nc,ir
        if self.spin==2:
            dip = (dip[:,:,:,:,:,0]+1j*dip[:,:,:,:,:,1]) # Read as ns,nk,nv,nc,ir
        dipoles = np.swapaxes(dip,self.spin,self.spin+2) # Swap indices as mentioned in the docstring
        database.close()

        return dipoles

    def readDB_oldformat(self,dip_type):
        """
        Legacy function for compatibility

        The dipole matrix has the following indexes:
        [nkpoints, cartesian directions, nspin, nbands conduction, nbands valence]
        """
        self.dip_type = dip_type
        dipoles = np.zeros([self.nk_ibz,3,self.nbandsc,self.nbandsv],dtype=np.complex64)
   
        #check dipole db format
        filename = "%s_fragment_1"%(self.filename)
        database = Dataset(filename)
        tag1 = 'DIP_iR_k_0001_spin_0001'
        tag2 = 'DIP_iR_k_0001_xyz_0001_spin_0001'
        if tag1 in list(database.variables.keys()):
            dipoles_format = 1
        elif tag2 in list(database.variables.keys()):
            dipoles_format = 2
        database.close()
        
        for nk in range(self.nk_ibz):

            #open database for each k-point
            filename = "%s_fragment_%d"%(self.filename,nk+1)
            database = Dataset(filename)

            if dipoles_format == 1:
                dip = database.variables['DIP_%s_k_%04d_spin_%04d'%(dip_type,nk+1,1)]
                dip = (dip[:,:,:,0]+1j*dip[:,:,:,1])
                for i in range(3):
                    dipoles[nk,i] = dip[:,:,i].T
            elif dipoles_format == 2:
                for i in range(3):
                    dip = database.variables['DIP_%s_k_%04d_xyz_%04d_spin_%04d'%(dip_type,nk+1,i+1,1)][:]
                    dipoles[nk,i] = dip[0].T+dip[1].T*1j

            #close database
            database.close()

        return dipoles
        
    def expandDipoles(self,dipoles=None,field_dir=[1,0,0],field_dir3=[0,0,1],spin=None):
        """
        Expand diples from the IBZ to the FBZ
        """
        if dipoles is None:
            dipoles = self.dipoles
            
        #check if we need to expand the dipoles to the full BZ
        lattice = self.lattice
        kpts = lattice.car_kpoints
        nks  = lattice.kpoints_indexes
        nss  = lattice.symmetry_indexes
        
        #normalize the fields
        field_dir  = np.array(field_dir)
        field_dir  = field_dir/np.linalg.norm(field_dir)
        field_dir3 = np.array(field_dir3)
        field_dir3 = field_dir3/np.linalg.norm(field_dir3)
        
        #calculate polarization directions
        field_dirx = field_dir
        field_diry = np.cross(field_dir3,field_dirx)
        field_dirz = field_dir3

        #get band indexes
        nkpoints = len(nks)
        indexv = self.index_firstv #self.min_band-1
        indexc = self.indexc       #self.indexc-1 
        nbandsv = self.nbandsv
        nbandsc = self.nbandsc
        nbands = self.min_band+self.nbands-1
        if self.open_shell: indexc  = self.indexc_os[spin]
        if self.open_shell: nbandsv = self.nbandsv_os[spin]
        if self.open_shell: nbandsc = self.nbandsc_os[spin]

        #Note that P is Hermitian and iR anti-hermitian.
        # [FP] Other possible dipole options (i.e., velocity gauge) to be checked. Treat them as not supported.
        if self.dip_type == 'P':
            factor =  1.0
        else:
            factor = -1.0
           
        #save dipoles in the ibz
        self.dipoles_ibz = dipoles 
        #get dipoles in the full Brillouin zone
        self.dipoles = np.zeros([nkpoints,3,nbands,nbands],dtype=np.complex64)
        for nk_fbz,nk_ibz,ns in zip(list(range(nkpoints)),nks,nss):
            
            #if time rev we conjugate
            if lattice.time_rev_list[ns]:
                dip = np.conjugate(dipoles[nk_ibz,:,:,:])
            else:
                dip = dipoles[nk_ibz,:,:,:]
            
            #get symmmetry operation
            sym = lattice.sym_car[ns].T
            #get projection operation
            pro = np.array([field_dirx,field_diry,field_dirz])
            #transformation
            tra = np.dot(pro,sym)
            
            #rotate dipoles
            for c,v in product(list(range(nbandsc)),list(range(nbandsv))):
                self.dipoles[nk_fbz,:,indexc+c,indexv+v] = np.dot(tra,dip[:,c,v])
        
            #make hermitian
            for c,v in product(list(range(nbandsc)),list(range(nbandsv))):
                self.dipoles[nk_fbz,:,indexv+v,indexc+c] = factor*np.conjugate(self.dipoles[nk_fbz,:,indexc+c,indexv+v])
                        
        self.field_dirx = field_dirx
        self.field_diry = field_diry
        self.field_dirz = field_dirz
        
        return self.dipoles, kpts
       
    def plot(self,ax,kpoint=0,dir=0,func=abs2):
        return ax.matshow(func(self.dipoles[kpoint,dir]))

    @add_fig_kwargs
    def plot_dipoles(self,data,nspin=-1,plt_show=False,plt_cbar=False,shift_BZ=True,**kwargs):
        """
        2D scatterplot in the k-BZ of the quantity A_{k}(is,ik,idir,ic,iv).
        TODO: this is the same function as plot_elph in elphondb. They should be merged.

        Any real quantity which is a function of only the k-grid may be supplied.
        The indices is,ik,idir,ic,iv are user-specified. 

        - if plt_show plot is shown
        - if plt_cbar colorbar is shown
        - if shift_BZ adjacent BZs are also plotted (default)
        - kwargs example: marker='H', s=300, cmap='viridis', etc.

        NB: So far requires a 2D system.
            Can be improved to plot BZ planes at constant k_z for 3D systems.
        """
        kpts = self.lattice.car_kpoints
        lattice = self.lattice
        rlat = self.lattice.rlat

        # Input check
        if len(data)!=len(kpts):
            raise ValueError('Something wrong in data dimensions (%d data vs %d kpts)'%(len(data),len(kpts)))

        # Global plot stuff
        self.fig, self.ax = plt.subplots(1, 1)
        self.ax.add_patch(BZ_Wigner_Seitz(lattice))

        if plt_cbar:
            if 'cmap' in kwargs.keys(): color_map = plt.get_cmap(kwargs['cmap'])
            else:                       color_map = plt.get_cmap('viridis')
        lim = 1.05*np.linalg.norm(rlat[0])
        self.ax.set_xlim(-lim,lim)
        self.ax.set_ylim(-lim,lim)

        # Reproduce plot also in adjacent BZs
        if shift_BZ:
            BZs = shifted_grids_2D(kpts,rlat)
            for kpts_s in BZs: plot=self.ax.scatter(kpts_s[:,0],kpts_s[:,1],c=data,**kwargs)
        else:
            plot=self.ax.scatter(kpts[:,0],kpts[:,1],c=data,**kwargs)

        if plt_cbar: self.cbar = self.fig.colorbar(plot)

        plt.gca().set_aspect('equal')

        if plt_show: plt.show()
        else: print("Plot ready.\nYou can customise adding savefig, title, labels, text, show, etc...")
        
    def ip_eps2(self,electrons,mode='imag',pol=1,ntot_dip=-1,nspin=-1,GWshift=0.,broad=0.1,broadtype='l',nbnds=[-1,-1],emin=0.,emax=10.,esteps=500,res_k=False):
        """
        Compute independent-particle absorption [interband transitions]

        electrons -> electrons YamboElectronsDB
        pol -> polarization direction(s). Can be integer or list of dirs to be summed over.
        ntot_dip -> if nbands_dip in ndb.dipoles < nbands_el in ns.db1, set ntot_dip=nbands_dip 
        nspin -> if -1 spin polarisations are summed (default)
                 if  0 only majority spin channel is considered
                 if  1 only minority spin channel is considered 
        GWshift -> rigid GW shift in eV
        broad -> broadening of peaks in eV
        broadtype -> 'l' is lorentzian, 'g' is gaussian
        nbnds -> number of [valence, conduction] bands included starting from Fermi level. Default means all are included
        emin,emax,esteps -> frequency range for the plot

        mode -> 'imag': Im[eps(w)] resonant case [DEFAULT] i.e. absorption spectrum / Fermi's golden rule
                'full': complex eps(w) including antiresonant case i.e. dielectric function / additional optical functions

        By R. Reho
        res_k -> if True, it returns an additional array epskres with IPA absorption for each k-point. 
                 In this way, we can plot it on the 2D-BZ (e.g. integrating over an energy range).

                This feature can be used like in the following example:

                :: code block ::
                    emin=0.
                    emax=3.5
                    step = int((emax-emin)/0.0025)
                    _, _, datakres = ydip.ip_eps2(yel,pol[0,1],ntot_dip=-1,broad=0.12,broadtype='l',emin=emin,emax=emax,nbnds=[2,2],esteps=step,res_k=True)
                    kres_int = np.sum(datakres,axis=0) #suitable integral over a frequency range
                    ydip.plot_dipoles(dataplot,marker='H',s=300,cmap='viridis')
               ::  end block ::
        """

        #get eigenvalues and weights of electrons
        eiv = electrons.eigenvalues_sp_pol
        weights = electrons.weights
        nv = electrons.nbandsv
        nc = electrons.nbandsc   
        nkpoints = len(eiv[0]) 

        #get dipoles
        dipoles = self.dipoles
        # (shape them according to spin polarization)
        sp_pol = self.spin # sum over all spin channels
        if self.spin==1: 
            dipoles = np.expand_dims(dipoles,axis=0)
        if self.spin==2 and (nspin==0 or nspin==1) : 
                dipoles = np.expand_dims(self.dipoles[nspin],axis=0)
                eiv = np.expand_dims(eiv[nspin],axis=0)
                sp_pol = self.spin-1 # sum over one spin channel only

        #get frequencies and im
        freq = np.linspace(emin,emax,esteps)
        if mode=='imag': eps = np.zeros([len(freq)])
        if mode=='full': eps = np.zeros([len(freq)],dtype=np.complex64)

        #Cut bands to the maximum number used for the dipoles
        if ntot_dip>0: 
            eiv = eiv[:,:,:ntot_dip]
            nc=ntot_dip-nv

        #Print band gap values and apply GW_shift
        electrons.energy_gaps(GWshift)

        #Check bands to include in the calculation
        if nbnds[0]<0: nbnds[0]=nv
        if nbnds[1]<0: nbnds[1]=nc
        iv = nv-nbnds[0] #first valence
        lc = nv+nbnds[1] #last conduction

        #choose broadening
        if mode=='imag' or res_k:
            if "l" in broadtype: broadening = lorentzian
            else:                broadening = gaussian

        #dimensional factors
        # [NB] This cofactor is not consistent with the yambo output:
        #      - In 3D there is a factor missing
        #      - In 2D there is a frequency dependence eps(w)->eps(w)/w missing (and a factor)
        #setting to 1. for now
        if self.spin == 1 : spin_deg=2
        if self.spin == 2 : spin_deg=1
        cofactor = spin_deg*8.*np.pi/self.lattice.rlat_vol
        cofactor = spin_deg*1.

        na = np.newaxis
        epskres = np.zeros([esteps,nkpoints])
        #calculate epsilon
        for c,v,s in product(range(nv,lc),range(iv,nv),range(sp_pol)):

                #get electron-hole energy and dipoles
                #(sum over pol directions if needed)
                eivs = eiv[s]
                ecv  = eivs[:,c]-eivs[:,v]

                dips = dipoles[s]
                dip2=0.
                try:
                    for p in pol: dip2 = dip2 + np.abs(dips[:,p,c,v])**2
                except TypeError: 
                    dip2 = np.abs(dips[:,pol,c,v])**2.

                #make dimensions match
                dip2a = dip2[na,:]
                ecva  = ecv[na,:]
                freqa = freq[:,na]
                wa    = weights[na,:]       

                if mode=='imag' or res_k: 
                    #calculate the lorentzians 
                    broadw = broadening(freqa,ecva,broad)
   
                    #scale broadening with dipoles and weights
                    epsk =  wa*dip2a*broadw

                    #k-resolved absorption
                    if res_k: epskres+=epsk

                    #integrate over kpoints
                    if mode=='imag': eps += np.sum(epsk,axis=1)

                if mode=='full':
                    #construct complex-valued response function
                    #including resonant and antiresonant components
                    G1 = -1./(freqa-ecva+broad*I)
                    G2 = -1./(-freqa-ecva-broad*I)

                    # oscillators
                    osc = wa*dip2a

                    # +=: sum over (c,v,s) ; np.sum(axis=1): sum over k                    
                    eps += np.sum(osc*(G1+G2),axis=1)/np.pi

        eps = eps*cofactor

        if res_k: return freq, eps, epskres
        else:     return freq, eps

    def add_drude(self,freq,eps,omegap,gammap):
        """
        Add 3D Drude term from semiclassical electron gas, i.e.,
        INTRABAND transitions to the dielectric function

        - freq, eps -> outputs of ip_eps2
                       - freq: energy window previously defined in ip_eps2
                       - eps: dielectric function previously computed with ip_eps2
                              * if float: assume it is Im[eps(w)] and add imaginary part of Drude term
                              * if cmplx: add real and imaginary parts of Drude term
        - omegap: plasma frequency
        - gammap:  Drude broadening

        Output: modified freq (removed E<=0 part) and eps with Drude added.
        """
        # Cut energy window
        zero_ind = np.where(freq==0.)[0]
        if zero_ind.size!=0:
            print('[WARNING] Values for energies <=0 are removed when adding the Drude term')
            freq = freq[zero_ind[0]+1:]
            eps  = eps[zero_ind[0]+1:]            

        # Drude correction, real and imaginary parts
        Drude_term  = omegap/(freq**2.+gammap**2.)
        Drude_real  = 1.-Drude_term
        Drude_imag  = Drude_term*gammap/freq
        Drude_cmplx = Drude_real+I*Drude_imag

        if np.issubdtype(eps.dtype, np.complexfloating): eps+=Drude_cmplx
        if np.issubdtype(eps.dtype, np.floating):        eps+=Drude_imag

        return freq,eps

    def __str__(self):
        lines = []; app = lines.append
        app(marquee(self.__class__.__name__))
        app("kpoints:")
        app("nk_ibz : %d"%self.nk_ibz)
        app("nk_bz  : %d"%self.nk_bz)
        app("bands:")
        app("nbands : %d" % self.nbands)
        app("nbandsv: %d" % self.nbandsv)
        app("nbandsc: %d" % self.nbandsc)
        app("indexv : %d" % (self.min_band-1))
        app("indexc : %d" % self.indexc)
        app("gauge  : %s" % (self.dip_type))
        app("spin:")
        app("spin pol         : %d" % (self.spin))
        if self.spin==2:
            app("open shell       : %s" % (self.open_shell))
            if self.open_shell: app("excess electrons : %d" % (self.n_exc_el))
        app("field_dirx: %10.6lf %10.6lf %10.6lf"%tuple(self.field_dirx))
        app("field_diry: %10.6lf %10.6lf %10.6lf"%tuple(self.field_diry))
        app("field_dirz: %10.6lf %10.6lf %10.6lf"%tuple(self.field_dirz))
        return "\n".join(lines)

if __name__ == "__main__":
    ddb = DipolesDB()
    ddb.get_databases()
    print(ddb)
