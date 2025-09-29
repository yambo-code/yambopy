#
# License-Identifier: GPL
#
# Copyright (C) 2024 The Yambo Team
#
# Authors: HPC,FP,RR
#
# This file is part of the yambopy project
#
import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from itertools import product 
import matplotlib.pyplot as plt
from yambopy import YamboLatticeDB
from yambopy.units import I
from yambopy.tools.types import CmplxType
from yambopy.tools.string import marquee
from yambopy.tools.funcs import abs2,lorentzian, gaussian
from yambopy.plot.plotting import add_fig_kwargs,BZ_Wigner_Seitz,shifted_grids_2D

class YamboDipolesDB():
    """
    Class to read the dipoles databases from the ``ndb.dipoles`` files
    
    Can be used to plot, for example, the imaginary part of the dielectric
    function which corresponds to the optical absorption, or directly the matrix elements in kspace.

    Dipole matrix elements <ck|vec{r}|vk> are stored in self.dipoles with indices [k,r_i,c,v]. 
    If the calculation is spin-polarised (nk->nks), then they are stored with indices [s,k,r_i,c,v]
    """
    def __init__(self,lattice,nq_ibz,nq_bz,nk_ibz,nk_bz,spin,min_band,max_band,indexv,indexc,bands_range,nbands,nbandsv,nbandsc,\
                      dip_bands_ordered,dipoles,dip_type,field_dir,project,polarization_mode,expand):
        """
        Initialize the YamboDipolesDB
        
        """        

        self.lattice   = lattice
        self.field_dir = field_dir
        self.nq_ibz    = nq_ibz
        self.nq_bz     = nq_bz
        self.nk_ibz    = nk_ibz
        self.nk_bz     = nk_bz
        self.spin      = spin
        self.min_band  = min_band # first band included in dipole calculation
        self.max_band  = max_band # last band included in dipole calculation
        self.indexv    = indexv   # index of maximum (partially) occupied band
        self.indexc    = indexc   # index of minimum (partially) empty band
        self.bands_range = bands_range # bands range in FORTRAN indexing 
        self.nbands    = nbands   # no. of bands in dipole calculation
        self.nbandsc   = nbandsc  # no. of conduction bands in dipole calculation
        self.nbandsv   = nbandsv  # no. of valence bands in dipole calculation
        self.dip_type  = dip_type
        # Standard dipoles: <c|e.r|v> between val and cond
        # Without dip_bands_ordered: <v|e.r|v'> and <c|e.r|c'> also present
        self.dip_bands_ordered = dip_bands_ordered
        # Additional operations to be done on dipoles
        self.project   = project
        self.expand    = expand
        # Explicitly select field polarzation for absorption plots (disables projection)
        self.polarization_mode = polarization_mode
        if self.polarization_mode is not None: self.project = False

        # This part is needed if dealing with open-shell systems
        self.index_firstv = self.min_band-1
        self.open_shell = False
        d_n_el = self.nbandsv + self.nbandsc - self.nbands
        if d_n_el != 0:
            # We assume the excess electron(s)/hole(s) are in the spin=0 channel
            print("[WARNING] You may be considering an open-shell system. Careful: bands_range option and optical absorption UNTESTED.")
            self.open_shell = True
            self.n_exc_el = d_n_el # this can be negative for excess holes
            # Spin-dependent band indices: these are taken from PARS so be careful
            self.indexv_os  = [self.indexv, self.indexv-d_n_el]
            self.indexc_os  = [self.indexc+d_n_el, self.indexc]
            self.nbandsv_os = [self.nbandsv, self.nbandsv-d_n_el ] 
            self.nbandsc_os = [self.nbandsc-d_n_el, self.nbandsc ]
        # End open-shell part

        # Dipoles matrix elements

        ## Note : Yambo stores dipoles are for light emission.
        ## In case of light absoption, conjugate it

        if not self.dip_bands_ordered:
            # use a slice of the full array
            self.dipoles_full = dipoles # Here keep the vv' and cc' transitions
            self.dipoles = dipoles[...,self.nbandsv:,:self.nbandsv] # cv only
        else:
            self.dipoles = dipoles

        # Expand the dipoles to the full brillouin zone 
        # and project them along field_dir
        if (self.expand):
            if self.spin==1: self.expandDipoles(self.dipoles)
            if self.spin==2:
                dip_up, dip_dn = self.dipoles[0], self.dipoles[1]
                exp_dip      = self.expandDipoles
                self.dipoles = np.stack((exp_dip(dip_up,spin=0,project=project)[0], exp_dip(dip_dn,spin=1)[0]),axis=0) 

    @classmethod
    def from_db_file(cls,lattice,filename='ndb.dipoles',dip_type='iR',field_dir=[1,1,1],bands_range=[],expand=True,project=True,polarization_mode=None,debug=False):
        """
        Initialize this class from a file

        The dipole matrix has the following indexes:
        [nspin, nkpoints, cartesian directions, nbands conduction, nbands valence]
      
        Parameters:
          -----------
        
        lattice: YamboLatticeDB
            Lattice information

        filename: str, optional
            Path and name of the database, default is 'ndb.dipoles'

        dip_type: str, optional
            Type of dipole matrix elements to read, can be 'iR', 'v', 'P', default is 'iR'

        field_dir: list of 3 floats, optional
            Direction of the electric field, default is [1,1,1]. Used if no polarization_mode is set. If project is True, applied already during k-expansion of dipoles.

        bands_range : array, optional
            select bands_range for computation of dipoles (Fortran indexing). Default: empty list []
            example: [7,10]  you consider from 7th (v) to the 10th (c) bands

        expand : k-expansion of dipoles (needed for eps). Default is True.

        project: bool, optional
            Whether to project the dipoles along the field direction during k-expansion, default is True

        polarization_mode: str, optional
            Polarization mode, can be 'linear', 'circular', 'dichroism', etc (check the docstring of related function; if set, `project` is turned to False and the chosen field polarization is applied for the epsilon calculation; field_dir is only used with the 'linear' option
            
        """
        if not isinstance(lattice,YamboLatticeDB):
            raise ValueError('Invalid type for lattice argument. It must be YamboLatticeDB')

        if not os.path.isfile(filename):
            raise FileNotFoundError("error opening %s in YamboLatticeDB"%filename)

        with Dataset(filename) as database:

            # General parameters
            nq_ibz, nq_bz, nk_ibz, nk_bz = database.variables['HEAD_R_LATT'][:].astype(int)
            spin = database.variables['SPIN_VARS'][0].astype(int)
            min_band, max_band, indexv, indexc = database.variables['PARS'][:4].astype(int)
            dip_bands_ordered = database.variables['PARS'][8].astype(int)
            
            # Determine the number of bands to read
            # We have four cases:
            # 1. full  range and     bands_ordered -> dipoles[0:Nv,0:Nc]
            # 2. full  range and not bands_ordered -> dipoles[1:(Nv+Nc),1:(Nv+Nc)]
            # 3. small range and     bands_ordered -> dipoles[i_v:Nv,1:f_c]
            # 4. small range and not bands_ordered -> dipoles[i_v:f_c,i_v:f_c]
            #
            # FP: I didn't want to change the previous class variable names for compatibility.
            #     However, it would be better to define unambiguously min_band and min_band_usr,
            #     nbands and nbands_usr and so on to distinguish between original and user bands

            # Cases 3. and 4.
            if len(bands_range) != 0:  # Custom selection of bands range
                if bands_range[0] not in range(min_band,indexv+1) or bands_range[1] not in range(indexc,max_band+1):
                    raise ValueError(f"[ERROR] invalid bands_range, db contains [{min_band},{max_band}]")
                
                min_band = min(bands_range)
                max_band = max(bands_range)
                nbands   = max_band-min_band+1  

                if dip_bands_ordered: # Standard case  
                    nbandsv = indexv-min_band+1
                    nbandsc = max_band-indexc+1
                    indexv = indexv-1
                    indexc = indexc-1 
                    nbands1, nbands2 = [nbandsv, nbandsc]
                    start_idx_v, start_idx_c = [bands_range[0]-1,0]
                    end_idx_v, end_idx_c = [indexv+1, nbandsc]

                if not dip_bands_ordered: # Yambo calculation with DipBandsALl
                    nbandsv = lattice.nbandsv-min_band+1
                    nbandsc = max_band-nbandsv
                    indexv  = nbandsv-1
                    indexc  = nbandsv
                    nbands1, nbands2 = [nbands, nbands]
                    start_idx_v, start_idx_c = [bands_range[0]-1,bands_range[0]-1]
                    end_idx_v, end_idx_c = [bands_range[1], bands_range[1]]

            # Cases 1. and 2.
            if len(bands_range) == 0:    # Read full database
                bands_range = [min_band,max_band]
                nbands      = max_band-min_band+1           
                
                if dip_bands_ordered: # Standard case
                    nbandsv = indexv-min_band+1
                    nbandsc = max_band-indexc+1
                    indexv = indexv-1
                    indexc = indexc-1
                    nbands1, nbands2 = [nbandsv, nbandsc]
                    start_idx_v, start_idx_c = [0,0]
                    end_idx_v, end_idx_c = [nbandsv, nbandsc]

                if not dip_bands_ordered: # Yambo calculation with DipBandsAll
                    nbandsv = lattice.nbandsv-min_band+1
                    nbandsc = max_band-nbandsv
                    indexv  = nbandsv-1
                    indexc  = nbandsv
                    nbands1, nbands2 = [nbands, nbands]
                    start_idx_v, start_idx_c = [0,0]
                    end_idx_v, end_idx_c = [nbands, nbands]

            if debug: # headache
                print(f"max_band: {max_band}")
                print(f"min_band: {min_band}")
                print(f"spin: {spin}")
                print(f"kpts: {nk_ibz}")
                print(f"bands range: {bands_range}")
                print(f"nbands: {nbands}")
                print(f"bandsv: {nbandsv}")
                print(f"bandsc: {nbandsc}")
                print(f"bands1: {nbands1}")
                print(f"bands2: {nbands2}")
                print(f"indexv: {indexv}")
                print(f"indexc: {indexc}")
                print(f"start_idx_v: {start_idx_v}")
                print(f"start_idx_c: {start_idx_c}")
                print(f"end_idx_v: {end_idx_v}")
                print(f"end_idx_c: {end_idx_c}")

            dipoles = database[f'DIP_{dip_type}'][:,:,start_idx_v:end_idx_v,start_idx_c:end_idx_c,:].data # Read as nk,nv,nc,ir
            dipoles = dipoles.view(dtype=CmplxType(dipoles)).reshape((spin,nk_ibz,nbands1,nbands2,3))

            if spin==1: dipoles = np.squeeze(dipoles,axis=0)
 
            dipoles = np.swapaxes(dipoles,spin,spin+2) # Swap indices as mentioned in the docstring

        return cls(lattice,nq_ibz,nq_bz,nk_ibz,nk_bz,spin,min_band,max_band,indexv,indexc,bands_range,nbands,nbandsv,nbandsc,dip_bands_ordered,\
                   dipoles,dip_type=dip_type,field_dir=field_dir,project=project,polarization_mode=polarization_mode,expand=expand)

    def normalize(self,electrons):
        """ 
        Use the electrons to normalize the dipole matrix elements
        """
        # We take the eivs with the added spin dimensions even in non-spin pol case
        eiv = electrons.eigenvalues
        nkpoints, nbands = eiv[0].shape

        dipoles = self.dipoles
        if self.spin==1: dipoles = np.expand_dims(dipoles,axis=0)
        
        for nk in range(nkpoints):
            for ns in range(self.spin):

                eiv_sk = eiv[ns,nk]
            
                #create eigenvalues differences arrays
                norm = np.array([ [ec-ev for ev in eiv_sk] for ec in eiv_sk ])
            
                #normalize
                for i,j in product(list(range(nbands)),repeat=2):

                    if norm[i,j] == 0: dipoles[ns,nk,:,i,j] = 0.
                    else: dipoles[ns,nk,:,i,j] = dipoles[ns,nk,:,i,j]/norm[i,j]

        if self.spin==1: dipoles=np.squeeze(dipoles,axis=0)
        self.dipoles = dipoles

    def expandDipoles(self,dipoles=None,spin=None):
        """
        Rotate dipoles from the IBZ to the FBZ
        and (if `project` is True) project them along field_dir
        [Equivalent to DIP_rotated and DIP_projected in Yambo]
        """
        if dipoles is None:
            dipoles = self.dipoles
            
        #check if we need to expand the dipoles to the full BZ
        lattice = self.lattice
        kpts = lattice.car_kpoints
        nks  = lattice.kpoints_indexes
        nss  = lattice.symmetry_indexes
        
        #normalize the fields
        self.field_dir  = np.array(self.field_dir)
        self.field_dir  = self.field_dir/np.linalg.norm(self.field_dir)
        
        #calculate polarization directions
        field_dirx = np.array([self.field_dir[0],0.,0.])
        field_diry = np.array([0.,self.field_dir[1],0.])
        field_dirz = np.array([0.,0.,self.field_dir[2]])

        #get band indexes
        nkpoints = len(nks)
        indexv = 0 # the dipoles array already starts from the min val band considered
        indexc = self.indexc - self.index_firstv # we have to start with conduction minus valence offset
        nbandsv = self.nbandsv
        nbandsc = self.nbandsc
        nbands = self.nbands
        if self.open_shell: 
            indexv = self.index_firstv
            indexc  = self.indexc_os[spin]
            nbandsv = self.nbandsv_os[spin]
            nbandsc = self.nbandsc_os[spin]

        assert self.dip_type in ['iR', 'v', 'P'], \
            "Dipole rotation is supported only for dip_type = iR, v, P."
        #
        # Note that P, v is Hermitian and iR anti-hermitian.
        if self.dip_type == 'iR': factor =  -1.0
        else:                     factor =  1.0
        
        ##
        #save dipoles in the ibz
        self.dipoles_ibz = dipoles 
        #get projection operation
        pro = np.array([field_dirx,field_diry,field_dirz])
        #get dipoles in the full Brillouin zone
        self.dipoles = np.zeros([nkpoints,3,self.nbnds_range,self.nbnds_range],dtype=dipoles.dtype)
        rot_mats = lattice.sym_car[nss, ...]
        if self.project: rot_mats = pro[None,:,:]@rot_mats
        # dipoles (nk, pol, c, v).
        dip_expanded = np.einsum('kij,kjcv->kicv',rot_mats,dipoles[nks])
        # Take care of time reversal
        trev = int(np.rint(lattice.time_rev))
        time_rev_s = (nss >= (len(lattice.sym_car) / (trev + 1)))
        dip_expanded[time_rev_s] = dip_expanded[time_rev_s].conj()
        #store them
        self.dipoles[:,:,self.nval_bands:self.nbnds_range,:self.nval_bands] = dip_expanded
        self.dipoles[:,:,:self.nval_bands,self.nval_bands:self.nbnds_range] = factor*dip_expanded.transpose(0,1,3,2).conj()
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
        
    def ip_eps2(self,electrons,mode='imag',ntot_dip=-1,nspin=-1,GWshift=0.,broad=0.1,broadtype='l',nbnds=[-1,-1],emin=0.,emax=10.,esteps=500,res_k=False,system_2D=False):
        """
        Compute independent-particle absorption [interband transitions]

        electrons -> electrons YamboElectronsDB over full BZ (Expand=True)
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
        
        2D_system -> if True, returns 2D polarizability instead of eps2
        res_k -> if True, it returns an additional array epskres with IPA absorption for each k-point. 
                 In this way, we can plot it on the 2D-BZ (e.g. integrating over an energy range).

                This feature can be used like in the following example:

                :: code block ::
                    emin=0.
                    emax=3.5
                    step = int((emax-emin)/0.0025)
                    _, _, datakres = ydip.ip_eps2(yel,ntot_dip=-1,broad=0.12,broadtype='l',emin=emin,emax=emax,nbnds=[2,2],esteps=step,res_k=True)
                    kres_int = np.sum(datakres,axis=0) #suitable integral over a frequency range
                    ydip.plot_dipoles(dataplot,marker='H',s=300,cmap='viridis')
               ::  end block ::
        """

        # Normalize field direction
        self.field_dir  = np.array(self.field_dir)
        self.field_dir  = self.field_dir/np.linalg.norm(self.field_dir)

        #get eigenvalues and weights of electrons
        if electrons.EXPAND == False:
            print("[WARNING] Expanding the electrons database")
            electrons.expandEigenvalues()
        eiv = electrons.eigenvalues[...,self.min_band-1:self.min_band-1+self.nbands]
        weights = electrons.weights
        nv = self.nbandsv # can be smaller than electrons.nbandsv
        nc = self.nbandsc # can be smaller than electrons.nbandsc
        nkpoints = len(eiv[0]) 

        #Print band gap values and apply GW_shift
        eiv[0]=electrons.energy_gaps(eiv[0],GWshift,nv=nv)

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
        if ntot_dip>0 and ntot_dip<self.nbands: 
            eiv = eiv[:,:,:ntot_dip]
            nc=ntot_dip-nv

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
        if self.spin == 1 : spin_deg=2
        if self.spin == 2 : spin_deg=1
        cofactor = spin_deg*8.*np.pi/(self.lattice.rlat_vol)

        na = np.newaxis
        epskres = np.zeros([esteps,nkpoints])
        #calculate epsilon
        for c,v,s in product(range(nv,lc),range(iv,nv),range(sp_pol)):

                #get electron-hole energy and dipoles
                #(sum over pol directions if needed)
                eivs = eiv[s]
                ecv  = eivs[:,c]-eivs[:,v]

                # these are the expanded+projected dipoles already
                dips = dipoles[s]
                
                # Initialize dip2 to accumulate contributions from all polarization vectors
                dip2 = np.zeros(nkpoints, dtype=np.float64)
                #RR: I believe the class should work without projecting but only via polarization_mode.
                #This series of if is there only to prevent back-compatibility problems.
                if self.project and self.polarization_mode is None:
                    dip2= np.abs( np.sum( dips[:,:,c,v], axis=1) )**2.
                
                elif not self.project and self.polarization_mode is None:
                    dip2 = np.abs( np.einsum('j,ij->i', self.field_dir , dips[:,:,c,v]) )**2

                if self.polarization_mode is not None: 
                    for e_vec, w in self._polarization_vectors():           # ⇐ NEW
                        #  ┌──── e_vec (3,) , dips_slice (nk,3)  ─────┐
                        proj = np.einsum('d,kd->k', e_vec, dips[:,:,c,v])   # shape (nk,)
                        dip2 += w*np.abs(proj)**2

                #make dimensions match
                dip2a = dip2[na,:]
                ecva  = ecv[na,:]
                freqa = freq[:,na]
                # rescale weight factors because we are in the expanded BZ
                wa    = weights[na,:]*self.nk_ibz/nkpoints       

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
       
        # Treat 2D case
        # THERE IS STILL A FACTOR MISSING
        #if system_2D:
        #    if mode=='imag':
        #        print("[2D system] Returning 2D polarizability instead of eps2 (bohr units)")
        #        # assuming supercell direction to be the largest one
        #        L = np.max(self.lattice.alat) 
        #        eps = eps#*L/4./np.pi
        #    if mode=='full':
        #        print("[2D system] Remember that for spectra, the physical quantity is eps.imag*L/(2*pi), with L aperiodic direction in bohr")

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

    def _polarization_vectors(self):
        """
        Return a list of (vector, weight) tuples based on the polarization mode.
        
        Each tuple contains:
        • vector  : complex ndarray(3,)   — The polarization direction in the Cartesian basis.
        • weight  : real scalar           — The weight of each contribution to the dielectric function.
        Polarization mode options:
        • 'linear'      : The polarization is defined by a user-specified direction.
        • 'unpolarized' : Averages over the three Cartesian directions: x, y, and z.
        • 'circular+'   : Circularly right polarized light in the xy/yz/xz planes. 
        • 'circular-'   : Circularly left polarized light in the xy/yz/xz planes.                   
        • 'dichroism'   : Difference between right and left circularly polarized light.                   
        """
        mode = self.polarization_mode.lower()

        if mode == 'linear':                   # ✱ ê set by user
            e = np.asarray(self.field_dir, dtype=np.complex64)
            e /= np.linalg.norm(e)
            return [(e, 1.0)]

        if mode == 'unpolarized':              # ✱ average over x,y,z
            return [ (np.array([1,0,0],np.complex64), 1/3),
                    (np.array([0,1,0],np.complex64), 1/3),
                    (np.array([0,0,1],np.complex64), 1/3) ]

        if mode == 'circularxy+':                # ✱ σ⁺  (propagation ‖ z)
            e = np.array([1, 1j, 0], np.complex64)/np.sqrt(2)
            return [(e, 1.0)]

        if mode == 'circularxy-':                # ✱ σ⁻
            e = np.array([1,-1j, 0], np.complex64)/np.sqrt(2)
            return [(e, 1.0)]
        if mode == 'circularxz+':                # ✱ σ⁺  (propagation ‖ z)
            e = np.array([1, 0, 1j], np.complex64)/np.sqrt(2)
            return [(e, 1.0)]

        if mode == 'circularxz-':                # ✱ σ⁻
            e = np.array([1,0, -1j], np.complex64)/np.sqrt(2)
            return [(e, 1.0)]
        if mode == 'circularyz+':                # ✱ σ⁺  (propagation ‖ z)
            e = np.array([0, 1, 1j], np.complex64)/np.sqrt(2)
            return [(e, 1.0)]

        if mode == 'circularyz-':                # ✱ σ⁻
            e = np.array([0, 1, -1j], np.complex64)/np.sqrt(2)
            return [(e, 1.0)]

        if mode == 'dichroism_xy':             # ✱ σ⁺ − σ⁻  (CD signal)
            e_p = np.array([1,  1j, 0], np.complex64)/np.sqrt(2)
            e_m = np.array([1, -1j, 0], np.complex64)/np.sqrt(2)
            return [(e_p,  1.0),              # add σ⁺
                    (e_m, -1.0)]              # subtract σ⁻
        if mode == 'dichroism_xz':             # ✱ σ⁺ − σ⁻  (CD signal)
            e_p = np.array([1,  0, 1j], np.complex64)/np.sqrt(2)
            e_m = np.array([1, 0, -1j], np.complex64)/np.sqrt(2)
            return [(e_p,  1.0),              # add σ⁺
                    (e_m, -1.0)]              # subtract σ⁻       
        if mode == 'dichroism_yz':             # ✱ σ⁺ − σ⁻  (CD signal)
            e_p = np.array([0, 1,  1j], np.complex64)/np.sqrt(2)
            e_m = np.array([0, 1, -1j], np.complex64)/np.sqrt(2)
            return [(e_p,  1.0),              # add σ⁺
                    (e_m, -1.0)]              # subtract σ⁻ 

        raise ValueError(f'Unknown polarization mode: {mode}')
    def __str__(self):
        lines = []; app = lines.append
        app(marquee(self.__class__.__name__))
        app("kpoints:")
        app("   nk_ibz : %d"%self.nk_ibz)
        app("   nk_bz  : %d"%self.nk_bz)
        app("bands:")
        app("   nbands : %d" % self.nbands)
        app("   nbandsv: %d" % self.nbandsv)
        app("   nbandsc: %d" % self.nbandsc)
        app("   indexv : %d" % (self.min_band-1))
        app("   indexc : %d" % self.indexc)
        app("gauge  : %s" % (self.dip_type))
        app("field_dir: %10.6lf %10.6lf %10.6lf"%tuple(self.field_dir))
        app("spin:")
        app("   spin pol         : %d" % (self.spin))
        if self.spin==2:
            app("   open shell       : %s" % (self.open_shell))
            if self.open_shell: app("   excess electrons : %d" % (self.n_exc_el))
        return "\n".join(lines)