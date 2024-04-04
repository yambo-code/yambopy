
from yambopy.units import *
from yambopy.plot.plotting import add_fig_kwargs,BZ_Wigner_Seitz
from yambopy.plot.bandstructure import *
from yambopy.lattice import replicate_red_kmesh, calculate_distances, get_path, car_red
from yambopy.tools.funcs import gaussian, lorentzian
from yambopy.tools.skw import SkwInterpolator
from yambopy.dbs.savedb import *
from yambopy.dbs.latticedb import *
from yambopy.dbs.electronsdb import *
from yambopy.dbs.dipolesdb import *
from yambopy.dbs.qpdb import *
from yambopy.tools.skw import SkwInterpolator

class ExcitonList():
    """
    Container class to perform operations on lists of excitons
    """
    def __init__(self,excitonlist):
        self.excitonlist = excitonlist

    def __str__(self):
        lines = []; app = lines.append
        for exciton in self.excitonlist:
            app( str(exciton.get_string(singleline=True)) )
        return "\n".join(lines)

class Exciton():
    """
    Basic container of data for a single exciton
    TODO: classify the excitons according to their symmetry
    """
    def __init__(self,energy,intensity,degeneracy,coeffs=None,wf=None):
        self.energy = energy
        self.intensity = intensity
        self.degeneracy = degeneracy
        self.coeffs = coeffs
        self.wf = wf

    def get_string(self,singleline=False):
        lines = []; app = lines.append
        app( 'energy:     %8.4'%self.energy )
        app( 'intensity:  %8.4'%self.intensity )
        app( 'degeneracy: %8d'%self.degeneracy )
        if singleline: return "".join(lines)
        return "\n".join(lines)
    
    def __str__(self):
        return self.get_string()

class YamboExcitonDB(YamboSaveDB):
    """ Read the excitonic states database from yambo

        Exciton eigenvectors are arranged as eigenvectors[i_exc, i_kvc]
        Transitions are unpacked in table[ i_k, i_v, i_c, i_s_c, i_s_v ] (last two are spin indices)
    """
    def __init__(
            self,lattice,Qpt,eigenvalues,l_residual,r_residual,spin_pol='no',car_qpoint=None,q_cutoff=None,table=None,eigenvectors=None, \
            pl_l_residual=None, pl_r_residual = None
            ):
        if not isinstance(lattice,YamboLatticeDB):
            raise ValueError('Invalid type for lattice argument. It must be YamboLatticeDB')

        self.Qpt = Qpt
        self.lattice = lattice
        self.eigenvalues = eigenvalues
        self.l_residual = l_residual
        self.r_residual = r_residual
        #optional
        self.car_qpoint = car_qpoint
        self.q_cutoff = q_cutoff
        self.table = table
        self.eigenvectors = eigenvectors
        self.spin_pol = spin_pol
        self.pl_l_residual = pl_l_residual
        self.pl_r_residual = pl_r_residual

    @classmethod
    def from_db_file(cls,lattice,filename='ndb.BS_diago_Q1',folder='.'):
        """ initialize this class from a file
        """
        path_filename = os.path.join(folder,filename)
        if not os.path.isfile(path_filename):
            raise FileNotFoundError("File %s not found in YamboExcitonDB"%path_filename)

        # Qpoint
        Qpt = filename.split("Q",1)[1]
        pl_l_residual = None
        pl_r_residual = None
        with Dataset(path_filename) as database:
            if 'BS_left_Residuals' in list(database.variables.keys()):
                #residuals
                rel,iml = database.variables['BS_left_Residuals'][:].T
                rer,imr = database.variables['BS_right_Residuals'][:].T
                l_residual = rel+iml*I
                r_residual = rer+imr*I
            if 'BS_Residuals' in list(database.variables.keys()):
                #residuals
                rel,iml,rer,imr = database.variables['BS_Residuals'][:].T
                l_residual = rel+iml*I
                r_residual = rer+imr*I

            if 'BS_PL_Residuals' in list(database.variables.keys()):
                #residuals
                replr,implr = database.variables['BS_PL_Residuals'][:].T
                pl_r_residual = replr+implr*I
                pl_l_residual = np.conjugate(pl_r_residual)

            car_qpoint = None
            if 'Q-point' in list(database.variables.keys()):
                # Finite momentum
                car_qpoint = database.variables['Q-point'][:]/lattice.alat
            if Qpt=="1": car_qpoint = np.zeros(3)

            #energies
            eig =  database.variables['BS_Energies'][:]*ha2ev
            eigenvalues = eig[:,0]+eig[:,1]*I
                
            #eigenvectors
            table = None
            eigenvectors = None
            if 'BS_EIGENSTATES' in database.variables:
                eiv = database.variables['BS_EIGENSTATES'][:]
                eiv = eiv[:,:,0] + eiv[:,:,1]*I
                eigenvectors = eiv
                table = database.variables['BS_TABLE'][:].T.astype(int)

            table = table
            eigenvectors = eigenvectors
            spin_vars = [int(database.variables['SPIN_VARS'][:][0]), int(database.variables['SPIN_VARS'][:][1])]
            if spin_vars[0] == 2 and spin_vars[1] == 1:
               spin_pol = 'pol'
            else:
               spin_pol = 'no'
        # Check if Coulomb cutoff is present
        path_cutoff = os.path.join(path_filename.split('ndb',1)[0],'ndb.cutoff')  
        q_cutoff = None
        if os.path.isfile(path_cutoff):
            with Dataset(path_cutoff) as database:
                bare_qpg = database.variables['CUT_BARE_QPG'][:]
                bare_qpg = bare_qpg[:,:,0]+bare_qpg[:,:,1]*I
                q_cutoff = np.abs(bare_qpg[0,int(Qpt)-1])

        return cls(lattice,Qpt,eigenvalues,l_residual,r_residual,spin_pol,q_cutoff=q_cutoff,car_qpoint=car_qpoint,table=table,eigenvectors=eigenvectors, 
                   pl_l_residual=pl_l_residual, pl_r_residual=pl_r_residual)

    @property
    def unique_vbands(self):
        return np.unique(self.table[:,1]-1)

    @property
    def unique_cbands(self):
        return np.unique(self.table[:,2]-1)

    @property
    def transitions_v_to_c(self):
        """Compute transitions from valence to conduction"""
        if hasattr(self,"_transitions_v_to_c"): return self._transitions_v_to_c
        uniq_v = self.unique_vbands
        uniq_c = self.unique_cbands
        transitions_v_to_c = dict([ ((v,c),[]) for v,c in product(uniq_v,uniq_c) ])

        #add elements to dictionary
        kidx = set()
        for eh,kvc in enumerate(self.table-1):
            k,v,c = kvc[0:3]
            kidx.add(k)
            transitions_v_to_c[(v,c)].append((k,eh))
        self.nkpoints = len(kidx)

        #make an array 
        for t,v in list(transitions_v_to_c.items()):
            if len(np.array(v)):
                transitions_v_to_c[t] = np.array(v)
            else:
                del transitions_v_to_c[t]

        self._transitions_v_to_c = transitions_v_to_c 
        return transitions_v_to_c

    @property
    def nkpoints(self): return max(self.table[:,0])

    @property
    def nvbands(self): return len(self.unique_vbands)

    @property
    def ncbands(self): return len(self.unique_cbands)

    @property
    def nbands(self): return self.ncbands+self.nvbands

    @property
    def mband(self): return max(self.unique_cbands)+1
 
    @property
    def ntransitions(self): return len(self.table)

    @property
    def nexcitons(self): return len(self.eigenvalues)
    
    @property
    def start_band(self): return min(self.unique_vbands)

    def write_sorted(self,prefix='yambo'):
        """
        Write the sorted energies and intensities to a file
        """
        #get intensities
        eig = self.eigenvalues.real
        intensities = self.get_intensities()

        #get sorted energies
        sort_e, sort_i = self.get_sorted()     

        #write excitons sorted by energy
        with open('%s_E.dat'%prefix, 'w') as f:
            for e,n in sort_e:
                f.write("%3d %12.8lf %12.8e\n"%(n+1,e,intensities[n])) 

        #write excitons sorted by intensities
        with open('%s_I.dat'%prefix,'w') as f:
            for i,n in sort_i:
                f.write("%3d %12.8lf %12.8e\n"%(n+1,eig[n],i)) 

    def get_nondegenerate(self,eps=1e-4):
        """
        get a list of non-degenerate excitons
        """
        non_deg_e   = [0]
        non_deg_idx = [] 

        #iterate over the energies
        for n,e in enumerate(self.eigenvalues):
            if not np.isclose(e,non_deg_e[-1],atol=eps):
                non_deg_e.append(e)
                non_deg_idx.append(n)

        return np.array(non_deg_e[1:]), np.array(non_deg_idx)

    def get_intensities(self):
        """
        get the intensities of the excitons
        """
        intensities = self.l_residual*self.r_residual
        intensities /= np.max(intensities)
        return intensities

    def get_sorted(self):
        """
        Return the excitonic weights sorted according to energy and intensity
        """
        #get intensities
        eig = self.eigenvalues.real
        intensities = self.get_intensities()

        #list ordered with energy
        sort_e = sorted(zip(eig, list(range(self.nexcitons))))

        #list ordered with intensity
        sort_i = sorted(zip(intensities, list(range(self.nexcitons))),reverse=True)

        return sort_e, sort_i 

    def get_degenerate(self,index,eps=1e-4):
        """
        Get degenerate excitons
        
        Args:
            eps: maximum energy difference to consider the two excitons degenerate in eV
        """
        energy = self.eigenvalues[index-1]
        excitons = [] 
        for n,e in enumerate(self.eigenvalues):
            if np.isclose(energy,e,atol=eps):
                excitons.append(n+1)
        return excitons

    def exciton_bs(self,energies,path,excitons=(0,),debug=False):
        """
        Calculate exciton band-structure
            
            Arguments:
            energies -> can be an instance of YamboSaveDB or YamboQBDB
            path     -> path in reduced coordinates in which to plot the band structure
            exciton  -> exciton index to plot
            spin     -> ??
        """
        if self.eigenvectors is None:
            raise ValueError('This database does not contain Excitonic states,'
                              'please re-run the yambo BSE calculation with the WRbsWF option in the input file.')
        if isinstance(excitons, int):
            excitons = (excitons,)
        #get full kmesh
        kpoints = self.lattice.red_kpoints
        path = np.array(path)

        rep = list(range(-1,2))
        kpoints_rep, kpoints_idx_rep = replicate_red_kmesh(kpoints,repx=rep,repy=rep,repz=rep)
        band_indexes = get_path(kpoints_rep,path,debug=debug)
        band_kpoints = kpoints_rep[band_indexes] 
        band_indexes = kpoints_idx_rep[band_indexes]

        if debug:
            import matplotlib.pyplot as plt
            for i,k in zip(band_indexes,band_kpoints):
                x,y,z = k
                plt.text(x,y,i) 
            plt.scatter(kpoints_rep[:,0],kpoints_rep[:,1])
            plt.plot(path[:,0],path[:,1],c='r')
            plt.scatter(band_kpoints[:,0],band_kpoints[:,1])
            plt.show()
            exit()

        #get eigenvalues along the path
        if isinstance(energies,(YamboSaveDB,YamboElectronsDB)):
            #expand eigenvalues to the full brillouin zone
            # SPIN-UP CHANNEL ONLY. Check with BSE WFs
            energies = energies.eigenvalues[0,self.lattice.kpoints_indexes]
            
        elif isinstance(energies,YamboQPDB):
            #expand the quasiparticle energies to the bull brillouin zone
            pad_energies = energies.eigenvalues_qp[self.lattice.kpoints_indexes]
            min_band = energies.min_band
            nkpoints, nbands = pad_energies.shape
            energies = np.zeros([nkpoints,energies.max_band])
            energies[:,min_band-1:] = pad_energies 
        else:
            raise ValueError("Energies argument must be an instance of YamboSaveDB,"
                             "YamboElectronsDB or YamboQPDB. Got %s"%(type(energies)))

        weights = self.get_exciton_weights(excitons)      
        energies = energies[band_indexes]
        weights  = weights[band_indexes]

        #make top valence band to be zero
        energies -= max(energies[:,max(self.unique_vbands)])
        
        return np.array(band_kpoints), energies, weights 

    def magnon_bs(self,energies,path,excitons=(0,),debug=False):
        """
        Calculate exciton band-structure
            
            Arguments:
            energies -> can be an instance of YamboSaveDB or YamboQBDB
            path     -> path in reduced coordinates in which to plot the band structure
            exciton  -> exciton index to plot
        """
        if self.eigenvectors is None:
            raise ValueError('This database does not contain Excitonic states,'
                              'please re-run the yambo BSE calculation with the WRbsWF option in the input file.')
        if isinstance(excitons, int):
            excitons = (excitons,)
        #get full kmesh
        kpoints = self.lattice.red_kpoints
        path = np.array(path)

        rep = list(range(-1,2))
        kpoints_rep, kpoints_idx_rep = replicate_red_kmesh(kpoints,repx=rep,repy=rep,repz=rep)
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
        print('self.nvbands')
        print(self.nvbands)
        print('self.ncbands')
        print(self.ncbands)
        exit()
        #get eigenvalues along the path
        if isinstance(energies,(YamboSaveDB,YamboElectronsDB)):
            #expand eigenvalues to the full brillouin zone
            # SPIN-UP Valence
            energies_v = energies.eigenvalues[0,self.lattice.kpoints_indexes,:self.nvbands]
            energies_c = energies.eigenvalues[1,self.lattice.kpoints_indexes,self.nvbands:(self.nvbands+self.ncbands)]


        elif isinstance(energies,YamboQPDB):
            #expand the quasiparticle energies to the bull brillouin zone
            pad_energies = energies.eigenvalues_qp[self.lattice.kpoints_indexes]
            min_band = energies.min_band
            nkpoints, nbands = pad_energies.shape
            energies = np.zeros([nkpoints,energies.max_band])
            energies[:,min_band-1:] = pad_energies 
        else:
            raise ValueError("Energies argument must be an instance of YamboSaveDB,"
                             "YamboElectronsDB or YamboQPDB. Got %s"%(type(energies)))

        weights = self.get_exciton_weights(excitons)      
        #print(energies.shape)
        #exit()
        energies = energies[band_indexes]
        weights  = weights[band_indexes]

        #make top valence band to be zero
        energies -= max(energies[:,max(self.unique_vbands)])
        
        return np.array(band_kpoints), energies, weights 


    def arpes_intensity(self,energies_db,path,excitons,ax):   #,size=1,space='bands',f=None,debug=False): later on
        size=1 # luego lo ponemos como input variable 
        n_excitons = len(excitons)
        #
        kpath   = path
        # kpoints IBZ
        kpoints = self.lattice.red_kpoints
        # array of high symmetry k-points
        path    = np.array(path.kpoints)

        # Expansion of IBZ kpoints to Path kpoints
        rep = list(range(-1,2))
        kpoints_rep, kpoints_idx_rep = replicate_red_kmesh(kpoints,repx=rep,repy=rep,repz=rep)
        band_indexes = get_path(kpoints_rep,path)
        band_kpoints = np.array(kpoints_rep[band_indexes])
        band_indexes = kpoints_idx_rep[band_indexes]

        # Eigenvalues Full BZ
        # Dimension nk_fbz x nbands
        energies = energies_db.eigenvalues[self.lattice.kpoints_indexes]

        # Calculate omega
        # omega_vk,lambda = e_(v,k-q) + omega_(lambda,q) only for q=0
        '''
        omega_vkl = np.zeros([self.nkpoints, self.nvbands,n_excitons])
        for i_l,exciton in enumerate(excitons):
            for i_k in range(self.nkpoints):
                for i_v in range(self.nvbands):
                    i_v2 = self.unique_vbands[i_v]
                    # omega_vk,lambda      = e_(v,k-q) + omega_(lambda,q)
                    omega_vkl[i_k,i_v,i_l] = energies[i_k,i_v2] + self.eigenvalues.real[exciton-1]

        '''
        omega_vkl = self.calculate_omega(energies,excitons)
        rho       = self.calculate_rho(excitons)
        # Calculate rho's
        # rho_vk = Sum_{c} |A_cvk|^2
#        rho = np.zeros([self.nkpoints, self.nvbands, n_excitons])


#        for i_exc, exciton in enumerate(excitons):
#            # get the eigenstate
#            eivec = self.eigenvectors[exciton-1]
#            for t,kvc in enumerate(self.table):
#                k,v,c = kvc[0:3]-1    # This is bug's source between yambo 4.4 and 5.0 check all this part of the class
#                i_v = v - self.nvbands                    # index de VB bands (start at 0)
#                i_c = c - self.ncbands - self.nvbands     # index de CB bands (start at 0)
#                rho[k,i_v,i_exc] += abs2(eivec[t])

        # Eigenvalues Path contains in Full BZ
        energies_path  = energies[band_indexes]
        rho_path       = rho[band_indexes]
        omega_vkl_path = omega_vkl[band_indexes]

        #make top valence band to be zero
        energies_path -= max(energies_path[:,max(self.unique_vbands)])

        plot_energies = energies_path[:,self.start_band:self.mband]
  
        # LDA or GW band structure
        ybs_bands = YambopyBandStructure(plot_energies, band_kpoints, kpath=kpath)


        # Intensity Plot
        print('shape energies_path')
        nkpoints_path=energies_path.shape[0]
        #exit()
        # Intensity histogram
        # I(k,omega_band)
        omega_band = np.arange(0.0,7.0,0.01)
        n_omegas = len(omega_band)
        Intensity = np.zeros([n_omegas,nkpoints_path]) 
        Im = 1.0j
           #for i_o in range(n_omegas):

        for i_o in range(n_omegas):
            for i_k in range(nkpoints_path):
                for i_v in range(self.nvbands):
                    for i_exc in range(n_excitons):
                        delta = 1.0/( omega_band[i_o] - omega_vkl_path[i_k,i_v,i_exc] + Im*0.2 )
                        Intensity[i_o,i_k] += rho_path[i_k,i_v,i_exc]*delta.imag

        distances = [0]
        distance = 0
        for nk in range(1,nkpoints_path):
            distance += np.linalg.norm(band_kpoints[nk]-band_kpoints[nk-1])
            distances.append(distance)
        distances = np.array(distances)
        X, Y = np.meshgrid(distances, omega_band)
        import matplotlib.pyplot as plt
        #plt.imshow(Intensity, interpolation='bilinear',cmap='viridis_r')
        plt.pcolor(X, Y, Intensity,cmap='viridis_r',shading='auto')
        # Excitonic Band Structure
        for i_v in range(self.nvbands):
            for i_exc in range(n_excitons):
                plt.plot(distances,omega_vkl_path[:,i_v,i_exc],color='w',lw=0.5) 
        # Electronic Band Structure
       
        for i_b in range(energies_db.nbands):
            plt.plot(distances,energies_path[:,i_b],lw=1.0,color='r')
        plt.xlim((distances[0],distances[-1]))
        plt.ylim((-5,10))
        plt.show()
        exit()

        # ARPES band structure
        ybs_omega = []
        for i_exc in range(n_excitons):
            plot_omega    = omega_vkl_path[:,:,i_exc]
            plot_rho      = rho_path[:,:,i_exc]
            size *= 1.0/np.max(plot_rho)
            ybs_omega.append( YambopyBandStructure(plot_omega, band_kpoints, weights=plot_rho, kpath=kpath, size=size) )

        # Plot bands
        ybs_bands.plot_ax(ax,color_bands='black',lw_label=2)

        for ybs in ybs_omega:
            ybs.plot_ax(ax,color_bands='black',lw_label=0.1)

        return rho

    def calculate_omega(self,energies,excitons):
        """ Calculate:
            omega_vk,lambda = e_(v,k-q) + omega_(lambda,q) only for q=0
        """

        n_excitons = len(excitons)
        omega_vkl = np.zeros([self.nkpoints, self.nvbands,n_excitons])
        for i_l,exciton in enumerate(excitons):
            for i_k in range(self.nkpoints):
                for i_v in range(self.nvbands):
                    i_v2 = self.unique_vbands[i_v]
                    # omega_vk,lambda      = e_(v,k-q) + omega_(lambda,q)
                    omega_vkl[i_k,i_v,i_l] = energies[i_k,i_v2] + self.eigenvalues.real[exciton-1]
         
        return omega_vkl

    def calculate_rho(self,excitons):
        """ Calculate:
            rho_vkl = Sum_{c} |A_cvk,l|^2
        """
        n_excitons = len(excitons)
        print('self.nkpoints, self.nvbands, n_excitons')
        print(self.nkpoints, self.nvbands, n_excitons)
        print('self.unique_vbands')
        print(self.unique_vbands)
        print('self.unique_cbands')
        print(self.unique_cbands)
        rho = np.zeros([self.nkpoints, self.nvbands, n_excitons])
        for i_exc, exciton in enumerate(excitons):
            # get the eigenstate
            eivec = self.eigenvectors[exciton-1]
            for t,kvc in enumerate(self.table):
                k,v,c = kvc[0:3]-1    # This is bug's source between yambo 4.4 and 5.0 check all this part of the class
                i_v = v - self.unique_vbands[0] # index de VB bands (start at 0)
                #i_c = c - self.unique_cbands[0] # index de CB bands (start at 0)
                rho[k,i_v,i_exc] += abs2(eivec[t])

        return rho

    #def arpes_interpolate(self,energies,path,excitons,lpratio=5,f=None,size=1,verbose=True,**kwargs):
    def arpes_intensity_interpolated(self,energies_db,path,excitons,lpratio=5,f=None,size=1,verbose=True,**kwargs):
        """ 
            Interpolate arpes bandstructure using SKW interpolation from Abipy (version 1)
            Change to the Fourier Transform Interpolation
            DFT energies == energies_db
            All is done internally. No use of the bandstructure class
            (something to change)
        """

        Im = 1.0j # Imaginary
        
        # Number of exciton states
        n_excitons = len(excitons)

        # Options kwargs

        # Alignment of the Bands Top Valence
        fermie      = kwargs.pop('fermie',0)
        # Parameters ARPES Intensity
        omega_width = kwargs.pop('omega_width',0)
        omega_1     = kwargs.pop('omega_1',0)
        omega_2     = kwargs.pop('omega_2',0)
        omega_step  = kwargs.pop('omega_step',0)
        omega_band  = np.arange(omega_1,omega_2,omega_step)
        n_omegas = len(omega_band)
        cmap_name   = kwargs.pop('cmap_name',0)
        scissor    = kwargs.pop('scissor',0)
       
        # Lattice and Symmetry Variables
        lattice = self.lattice
        cell = (lattice.lat, lattice.red_atomic_positions, lattice.atomic_numbers)

        symrel = [sym for sym,trev in zip(lattice.sym_rec_red,lattice.time_rev_list) if trev==False ]
        time_rev = True

        nelect = 0  # Why?

        # DFT Eigenvalues FBZ
        energies = energies_db.eigenvalues[0,self.lattice.kpoints_indexes] #SPIN-UP
        # Rho FBZ
        rho      = self.calculate_rho(excitons)
        if f: rho = f(rho)
        # Omega FBZ
        omega    = self.calculate_omega(energies,excitons)

        size *= 1.0/np.max(rho)

        ibz_nkpoints = max(lattice.kpoints_indexes)+1
        kpoints = lattice.red_kpoints

        #map from bz -> ibz:
        ibz_rho     = np.zeros([ibz_nkpoints,self.nvbands,n_excitons])
        ibz_kpoints = np.zeros([ibz_nkpoints,3])
        ibz_omega   = np.zeros([ibz_nkpoints,self.nvbands,n_excitons])
        for idx_bz,idx_ibz in enumerate(lattice.kpoints_indexes):
            ibz_rho[idx_ibz,:,:]   = rho[idx_bz,:,:] 
            ibz_kpoints[idx_ibz]   = lattice.red_kpoints[idx_bz]
            ibz_omega[idx_ibz,:,:] = omega[idx_bz,:,:] 

        #get DFT or GW eigenvalues
        if isinstance(energies_db,(YamboSaveDB,YamboElectronsDB)):
            ibz_energies = energies_db.eigenvalues[0,:,self.start_band:self.mband] #spin-up
        elif isinstance(energies_db,YamboQPDB):   # Check this works !!!!
            ibz_energies = energies_db.eigenvalues_qp
        else:
            raise ValueError("Energies argument must be an instance of YamboSaveDB,"
                             "YamboElectronsDB or YamboQPDB. Got %s"%(type(energies)))

        # set k-path
        kpoints_path = path.get_klist()[:,:3]
        distances = calculate_distances(kpoints_path)
        nkpoints_path = kpoints_path.shape[0]

        na = np.newaxis
        rho_path   = np.zeros([1, nkpoints_path, self.nvbands, n_excitons])
        omega_path = np.zeros([1, nkpoints_path, self.nvbands, n_excitons])

        for i_exc in range(n_excitons):

            # interpolate rho along the k-path
            skw_rho   = SkwInterpolator(lpratio,ibz_kpoints,ibz_rho[na,:,:,i_exc],fermie,nelect,cell,symrel,time_rev,verbose=verbose)
            rho_path[0,:,:,i_exc] = skw_rho.interp_kpts(kpoints_path).eigens

            # interpolate omega along the k-path
            skw_omega = SkwInterpolator(lpratio,ibz_kpoints,ibz_omega[na,:,:,i_exc],fermie,nelect,cell,symrel,time_rev,verbose=verbose)
            omega_path[0,:,:,i_exc] = skw_omega.interp_kpts(kpoints_path).eigens

        # interpolate energies
        skw_energie = SkwInterpolator(lpratio,ibz_kpoints,ibz_energies[na,:,:],fermie,nelect,cell,symrel,time_rev,verbose=verbose)
        energies_path = skw_energie.interp_kpts(kpoints_path).eigens

        top_valence_band = np.max(energies_path[0,:,0:self.nvbands])
        omega_path    = omega_path - top_valence_band
        energies_path = energies_path - top_valence_band

        import matplotlib.pyplot as plt

        # I(k,omega_band)
        Intensity = np.zeros([n_omegas,nkpoints_path]) 
         
        for i_exc in range(n_excitons):
            for i_o in range(n_omegas):
                for i_k in range(nkpoints_path):
                    for i_v in range(self.nvbands):
                        delta = 1.0/( omega_band[i_o] - omega_path[0, i_k, i_v, i_exc] + Im*omega_width ) # check this
                        Intensity[i_o,i_k] += rho_path[0, i_k, i_v, i_exc]*delta.imag

        X, Y = np.meshgrid(distances, omega_band)
        import matplotlib.pyplot as plt

        # Plot I(k,w)
        plt.pcolor(X, Y, Intensity,cmap=cmap_name,shading='auto')


        # Plot Excitonic Energies
        #for i_exc in range(n_excitons):
        #    for i_v in range(self.nvbands):
        #        plt.plot(distances,omega_path[0,:,i_v,i_exc],color='white',lw=0.5)

        # Plot Valence Band Energies
        #for i_b in range(energies_path.shape[2]):
        for i_b in range(self.nvbands):
            plt.plot(distances,energies_path[0,:,i_b],lw=1.0,color='white')
        for i_b in range(self.ncbands):
            plt.plot(distances,energies_path[0,:,i_b+self.nvbands]+scissor,lw=1.0,color='white')
            #plt.plot(distances,energies_path[0,:,i_b+9],lw=0.5,color='r')

        plt.xlim((distances[0],distances[-1]))
        plt.ylim((omega_1,omega_2-omega_width))

        #plt.axhline(np.max(energies_path[0,:,0:self.nvbands]),c='white')
        for kpoint, klabel, distance in path:
            plt.axvline(distance,c='w')
        plt.xticks(path.distances,path.klabels)
        plt.show()

        #create band-structure object
        #exc_bands = YambopyBandStructure(energies[0],kpoints_path,kpath=path,weights=exc_weights[0],size=size,**kwargs)
        #exc_bands.set_fermi(self.nvbands)
        #exit()
        #return exc_bands
        return 

    def get_exciton_weights(self,excitons):
        """get weight of state in each band"""
        weights = np.zeros([self.nkpoints,self.mband])
        for exciton in excitons:
            #get the eigenstate
            eivec = self.eigenvectors[exciton-1]

            #add weights
            sum_weights = 0
            for t,kcv in enumerate(self.table):
                k,c,v = kcv[0:3]-1
                this_weight = abs2(eivec[t])
                weights[k,c] += this_weight
                weights[k,v] += this_weight
                sum_weights += this_weight
            if abs(sum_weights - 1) > 1e-3: raise ValueError('Excitonic weights does not sum to 1 but to %lf.'%sum_weights)
 
        return weights
    
    def get_exciton_total_weights(self,excitons):
        """get weight of state in each band"""
        total_weights = np.zeros(self.nkpoints)
        for exciton in excitons:
            #get the eigenstate
            eivec = self.eigenvectors[exciton-1]
            #add weights
            sum_weights = 0
            for t,kcv in enumerate(self.table):
                k,c,v = kcv[0:3]
                total_weights[k-1] += abs2(eivec[t])
            if abs(sum(total_weights) - 1) > 1e-3: raise ValueError('Excitonic weights does not sum to 1 but to %lf.'%sum_weights)
 
        return total_weights

    def get_exciton_transitions(self,excitons):
        """get weight of state in each band"""
        # Double check the part of the array w_k_v_to_c
        # We should comment more this part
        #weights = np.zeros([self.nkpoints,self.mband])
        w_k_v_to_c = np.zeros([self.nkpoints,self.nvbands,self.ncbands])
        v_min = self.unique_vbands[0]
        c_min = self.unique_cbands[0]
        for exciton in excitons:
            #get the eigenstate
            eivec = self.eigenvectors[exciton-1]
            #add weights
            #sum_weights = 0
            for t,kcv in enumerate(self.table):
                k,c,v = kcv-1                             
                #k,v,c = kcv-1                 # bug?? Double-check
                this_weight = abs2(eivec[t])
                w_k_v_to_c[k,v-v_min,c-c_min] = this_weight   # new
            #if abs(sum_weights - 1) > 1e-3: raise ValueError('Excitonic weights does not sum to 1 but to %lf.'%sum_weights)
 
        #return weights, w_k_v_to_c
        return w_k_v_to_c

    def get_exciton_2D(self,excitons,f=None):
        """get data of the exciton in 2D"""
        weights = self.get_exciton_weights(excitons)
        #sum all the bands
        weights_bz_sum = np.sum(weights,axis=1)
        if f: weights_bz_sum = f(weights_bz_sum)

        kmesh_full, kmesh_idx = replicate_red_kmesh(self.lattice.red_kpoints,repx=range(-1,2),repy=range(-1,2))
        x,y = red_car(kmesh_full,self.lattice.rlat)[:,:2].T
        weights_bz_sum = weights_bz_sum[kmesh_idx]
        return x,y,weights_bz_sum

    def get_exciton_lifetimes(self, statelist = None, degen_step = 0.001, gauge = 'length', verbosity = 'False', pl_res = False):
        """get exciton lifetimes in 2D assuming the plane is on x,y
            TODO: T different than 0
        """
        if (statelist is None ):
            statelist = np.arange(1,self.nexcitons+1)
        Omega = self.lattice.lat_vol
        a1 = self.lattice.lat[0]
        a2 = self.lattice.lat[1]
        A = np.linalg.norm(np.cross(a1,a2))
        muS2 = 0 
        tmpE = np.ma.asarray(self.eigenvalues.real)
        if (pl_res):
            gauge = 'velocity' # force velocity gauge for PL
            tmpI = np.ma.asarray((self.pl_l_residual*self.pl_r_residual))
        else:
            tmpI = np.ma.asarray(self.l_residual * self.r_residual)
        excE = sorted(tmpE)
        excI = tmpI

        tau0_tot = np.zeros(len(statelist))
        merged_states = np.empty(len(statelist), dtype='object')
        # q0 norm factor 
        q0_norm = 1e-5

        for l, state in enumerate(statelist):
            state_internal = state 

            # find states within degen_step window from state
            mesk = np.logical_and(excE>=excE[state_internal]-degen_step,excE<=excE[state_internal]+degen_step)
            states_idx = np.where(mesk==True)[0]
            #compute gamma
            gamma0 = 0
            for i, st in enumerate(states_idx):
                ES = excE[st]/ha2ev
                if (gauge=='length'):
                    muS2 = excI[st]/(q0_norm**2) # if you inspect the Yambo code you might expect another 1/((2*np.pi)**3) but I think that d3k_factor/((2np.pi)**3) is actually 1/Omega
                elif (gauge == 'velocity'):
                    muS2 = excI[st]/(ES**2)
                gg = 4.*np.pi*ES*(muS2/self.nkpoints)/(A*speed_of_light)
                gamma0 += gg
            #compute tau
            tau0_tot[l] = AU2S/(gamma0.real)
            merged_states[l] = '{}<->{}'.format(min(states_idx)+1,max(states_idx)+1)
        
        if (verbosity):
            for l,state in enumerate(statelist):
                print(' Exciton radiative lifetime')
                print('#')
                #print(f'# Effective mass = {m_e}\n')
                print(f"# Global Gauge = {gauge}")
                print(f'Energy degeneration step = {degen_step} eV')
                print('#')
                print(f"# {'State':>5} {'Energy (eV)':>15} {'$tau_0$ (s)':>20} {'Merged states':>20}")
                print(f"{state:>5} {excE[state]:>15} {tau0_tot[l]:>20} {merged_states[l]:>20}")
        return excE, tau0_tot, merged_states
    
    def get_exciton_2D_spin_pol(self,excitons,f=None):
        """get data of the exciton in 2D for spin polarized calculations"""
        weights_up, weights_dw = self.get_exciton_weights_spin_pol(excitons)

        #sum all the bands
        weights_bz_sum_up = np.sum(weights_up,axis=1)
        weights_bz_sum_dw = np.sum(weights_dw,axis=1)

        if f: weights_bz_sum_up = f(weights_bz_sum_up)
        if f: weights_bz_sum_dw = f(weights_bz_sum_dw)

        kmesh_full, kmesh_idx = replicate_red_kmesh(self.lattice.red_kpoints,repx=range(-1,2),repy=range(-1,2))
        x,y = red_car(kmesh_full,self.lattice.rlat)[:,:2].T
        weights_bz_sum_up = weights_bz_sum_up[kmesh_idx]
        weights_bz_sum_dw = weights_bz_sum_dw[kmesh_idx]

        return x,y,weights_bz_sum_up,weights_bz_sum_dw
 
 
    def plot_exciton_2D_ax(self,ax,excitons,f=None,mode='hexagon',limfactor=0.8,spin_pol=None,**kwargs):
        """
        Plot the exciton weights in a 2D Brillouin zone
       
           Arguments:
            excitons -> list of exciton indexes to plot
            f -> function to apply to the exciton weights. Ex. f=log will compute the 
                 log of th weight to enhance the small contributions
            mode -> possible values are 'hexagon'/'square' to use hexagons/squares as markers for the 
                    weights plot and 'rbf' to interpolate the weights using radial basis functions.
            limfactor -> factor of the lattice parameter to choose the limits of the plot 
            scale -> size of the markers
        """
        if spin_pol is not None: print('Plotting exciton mad in 2D axis for spin polarization: %s' % spin_pol)

        if spin_pol is not None:
           x,y,weights_bz_sum_up,weights_bz_sum_dw = self.get_exciton_2D_spin_pol(excitons,f=f)
        else:
           x,y,weights_bz_sum = self.get_exciton_2D(excitons,f=f)

        weights_bz_sum=weights_bz_sum/np.max(weights_bz_sum)
        
        #filter points outside of area
        lim = np.max(self.lattice.rlat)*limfactor
        dlim = lim*1.1
        if spin_pol is not None:
           filtered_weights_up = [[xi,yi,di] for xi,yi,di in zip(x,y,weights_bz_sum_up) if -dlim<xi and xi<dlim and -dlim<yi and yi<dlim]
           filtered_weights_dw = [[xi,yi,di] for xi,yi,di in zip(x,y,weights_bz_sum_dw) if -dlim<xi and xi<dlim and -dlim<yi and yi<dlim]
           x,y,weights_bz_sum_up = np.array(filtered_weights_up).T
           x,y,weights_bz_sum_dw = np.array(filtered_weights_dw).T
        else:
           filtered_weights = [[xi,yi,di] for xi,yi,di in zip(x,y,weights_bz_sum) if -dlim<xi and xi<dlim and -dlim<yi and yi<dlim]
           x,y,weights_bz_sum = np.array(filtered_weights).T
        # Add contours of BZ
        ax.add_patch(BZ_Wigner_Seitz(self.lattice))

        #plotting
        if mode == 'hexagon': 
            scale = kwargs.pop('scale',1)
            if spin_pol=='up':
               s=ax.scatter(x,y,s=scale,marker='H',c=weights_bz_sum_up,rasterized=True,**kwargs)
            elif spin_pol=='dw':
               s=ax.scatter(x,y,s=scale,marker='H',c=weights_bz_sum_dw,rasterized=True,**kwargs)
            else:
               s=ax.scatter(x,y,s=scale,marker='H',c=weights_bz_sum,rasterized=True,**kwargs)
            ax.set_xlim(-lim,lim)
            ax.set_ylim(-lim,lim)
        elif mode == 'square': 
            scale = kwargs.pop('scale',1)
            if spin_pol=='up':
               s=ax.scatter(x,y,s=scale,marker='s',c=weights_bz_sum_up,rasterized=True,**kwargs)
            elif spin_pol=='dw':
               s=ax.scatter(x,y,s=scale,marker='s',c=weights_bz_sum_dw,rasterized=True,**kwargs)
            else:
               s=ax.scatter(x,y,s=scale,marker='s',c=weights_bz_sum,rasterized=True,**kwargs)
            ax.set_xlim(-lim,lim)
            ax.set_ylim(-lim,lim)
        elif mode == 'rbf':
            from scipy.interpolate import Rbf
            npts = kwargs.pop('npts',100)
            interp_method = kwargs.pop('interp_method','bicubic')
            if spin_pol=='up':
               rbfi = Rbf(x,y,weights_bz_sum_up,function='linear')
               x = y = np.linspace(-lim,lim,npts)
               weights_bz_sum_up = np.zeros([npts,npts])
            elif spin_pol=='dw':
               rbfi = Rbf(x,y,weights_bz_sum_dw,function='linear')
               x = y = np.linspace(-lim,lim,npts)
               weights_bz_sum_dw = np.zeros([npts,npts])
            else:
               rbfi = Rbf(x,y,weights_bz_sum,function='linear')
               x = y = np.linspace(-lim,lim,npts)
               weights_bz_sum = np.zeros([npts,npts])

            for col in range(npts):
                if spin_pol=='up':
                   weights_bz_sum_up[:,col] = rbfi(x,np.ones_like(x)*y[col])
                elif spin_pol=='dw':
                   weights_bz_sum_dw[:,col] = rbfi(x,np.ones_like(x)*y[col])
                else:
                   weights_bz_sum[:,col] = rbfi(x,np.ones_like(x)*y[col])
            # NB we have to take the transpose of the imshow data to get the correct plot
            if spin_pol=='up':
               s=ax.imshow(weights_bz_sum_up.T,interpolation=interp_method,extent=[-lim,lim,-lim,lim])
            elif spin_pol=='dw':
               s=ax.imshow(weights_bz_sum_dw.T,interpolation=interp_method,extent=[-lim,lim,-lim,lim])
            else:
               s=ax.imshow(weights_bz_sum.T,interpolation=interp_method,extent=[-lim,lim,-lim,lim])
        title = kwargs.pop('title',str(excitons))

        ax.set_title(title)
        ax.set_aspect('equal')
        ax.set_xticks([])
        ax.set_yticks([])
       
        return ax,s

    def get_exciton_3D(self,excitons,f=None):
        """get data of the exciton in 2D"""
        weights = self.get_exciton_weights(excitons)
        #sum all the bands
        weights_bz_sum = np.sum(weights,axis=1)
        if f: weights_bz_sum = f(weights_bz_sum)

        kmesh_full, kmesh_idx = replicate_red_kmesh(self.lattice.red_kpoints,repx=range(-1,2),repy=range(-1,2))
        x,y,z = red_car(kmesh_full,self.lattice.rlat)[:,:3].T
        weights_bz_sum = weights_bz_sum[kmesh_idx]
        return x,y,z,weights_bz_sum
 
    def plot_exciton_3D_ax(self,ax,excitons,f=None,mode='hexagon',limfactor=0.8,**kwargs):
        """
        Plot the exciton weights in a 3D Brillouin zone
       
           Arguments:
            excitons -> list of exciton indexes to plot
            f -> function to apply to the exciton weights. Ex. f=log will compute the 
                 log of th weight to enhance the small contributions
            mode -> possible values are 'hexagon' to use hexagons as markers for the 
                    weights plot and 'rbf' to interpolate the weights using radial basis functions.
            limfactor -> factor of the lattice parameter to choose the limits of the plot 
            scale -> size of the markers
        """
        x,y,z,weights_bz_sum = self.get_exciton_3D(excitons,f=f)
        print(x,y,z,weights_bz_sum)
        exit()

        #filter points outside of area
        lim = np.max(self.lattice.rlat)*limfactor
        dlim = lim*1.1
        filtered_weights = [[xi,yi,di] for xi,yi,di in zip(x,y,weights_bz_sum) if -dlim<xi and xi<dlim and -dlim<yi and yi<dlim]
        x,y,z,weights_bz_sum = np.array(filtered_weights).T

        #plotting
        if mode == 'hexagon': 
            scale = kwargs.pop('scale',1)
            ax.scatter(x,y,s=scale,marker='H',c=weights_bz_sum,rasterized=True,**kwargs)
            ax.set_xlim(-lim,lim)
            ax.set_ylim(-lim,lim)
        elif mode == 'square': 
            scale = kwargs.pop('scale',1)
            ax.scatter(x,y,s=scale,marker='s',c=weights_bz_sum,rasterized=True,**kwargs)
            ax.set_xlim(-lim,lim)
            ax.set_ylim(-lim,lim)
        elif mode == 'rbf':
            from scipy.interpolate import Rbf
            npts = kwargs.pop('npts',100)
            interp_method = kwargs.pop('interp_method','bicubic')
            rbfi = Rbf(x,y,weights_bz_sum,function='linear')
            x = y = np.linspace(-lim,lim,npts)
            weights_bz_sum = np.zeros([npts,npts])
            for col in range(npts):
                weights_bz_sum[:,col] = rbfi(x,np.ones_like(x)*y[col])
            # NB we have to take the transpose of the imshow data to get the correct plot
            ax.imshow(weights_bz_sum.T,interpolation=interp_method,extent=[-lim,lim,-lim,lim])
        title = kwargs.pop('title',str(excitons))
        
        ax.set_title(title)
        ax.set_aspect('equal')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.add_patch(BZ_Wigner_Seitz(self.lattice))
        return ax

    def plot_nbrightest_2D(self,emin=0,emax=10,estep=0.001,broad=0.1,
                           mode='rbf',scale=3,nrows=2,ncols=2,eps=1e-5):
        """
        Create a plot with chi and vertical bars for the brightest excitons
        Also plot the 2D wavefunctions of the brightest excitons.

          Arguments:
            emin,emax -> minimum and maximum energy range to plot chi
            estep -> energy step to plot chi
            broad -> broadening of the exciton peaks
            mode -> possible values are 'hexagon' to use hexagons as markers for the 
                    weights plot and 'rbf' to interpolate the weights using radial basis functions.
            scale -> size of the markers
            nrows,ncols -> number of rows and colums for the 2D plots (default: 2x2)
            eps -> threshold to find degenerate states
        """
        import matplotlib.pyplot as plt
        figexc = plt.figure()
        n_brightest = nrows*ncols
        figchi = self.plot_chi(emin=emin,emax=emax,estep=estep,broad=broad,n_brightest=n_brightest,show=False)
        #plot vertical bar on the brightest excitons
        exc_e,exc_i = self.get_sorted()
        sorted_exc = sorted(exc_i[:n_brightest],key = lambda x: x[1])
        for n,(i,idx) in enumerate(sorted_exc):
            ax = figexc.add_subplot(nrows,ncols,n+1)
            excitons = self.get_degenerate(idx,eps)
            self.plot_exciton_2D_ax(ax,excitons,scale=scale,mode=mode)
        return figchi,figexc

    def get_exciton_bs(self,energies_db,path,excitons,size=1,space='bands',f=None,debug=False):
        """
        Get a YambopyBandstructure object with the exciton band-structure
        
            Arguments:
            ax          -> axis extance of matplotlib to add the plot to
            lattice     -> Lattice database
            energies_db -> Energies database, can be either a SaveDB or QPDB
            path        -> Path in the brillouin zone
        """
        from qepy.lattice import Path
        if not isinstance(path,Path): 
            raise ValueError('Path argument must be a instance of Path. Got %s instead'%type(path))
    
        if space == 'bands':
            if self.spin_pol=='no':
               bands_kpoints, energies, weights = self.exciton_bs(energies_db, path.kpoints, excitons, debug)
               nkpoints = len(bands_kpoints)
               plot_energies = energies[:,self.start_band:self.mband]
               plot_weights  = weights[:,self.start_band:self.mband]
        #    elif spin_pol=='pol':
               
        else:
            raise NotImplementedError('TODO')
            eh_size = len(self.unique_vbands)*len(self.unique_cbands)
            nkpoints = len(bands_kpoints)
            plot_energies = np.zeros([nkpoints,eh_size])
            plot_weights = np.zeros([nkpoints,eh_size])
            for eh,(v,c) in enumerate(product(self.unique_vbands,self.unique_cbands)):
                plot_energies[:,eh] = energies[:,c]-energies[:,v]
                plot_weights[:,eh] = weights[:,c] 

        if f: plot_weights = f(plot_weights)
        size *= 1.0/np.max(plot_weights)
        ybs = YambopyBandStructure(plot_energies, bands_kpoints, weights=plot_weights, kpath=path, size=size)
        return ybs

    def get_magnon_bs(self,energies_db,path,excitons,size=1,space='bands',f=None,debug=False):
        """
        Get a YambopyBandstructure object with the exciton band-structure
        
            Arguments:
            ax          -> axis extance of matplotlib to add the plot to
            lattice     -> Lattice database
            energies_db -> Energies database, can be either a SaveDB or QPDB
            path        -> Path in the brillouin zone
        """
        from qepy.lattice import Path
        if not isinstance(path,Path): 
            raise ValueError('Path argument must be a instance of Path. Got %s instead'%type(path))
    
        if space == 'bands':
            bands_kpoints, energies, weights = self.magnon_bs(energies_db, path.kpoints, excitons, debug)
            nkpoints = len(bands_kpoints)
            plot_energies = energies[:,self.start_band:self.mband]
            plot_weights  = weights[:,self.start_band:self.mband]
        else:
            raise NotImplementedError('TODO')
            eh_size = len(self.unique_vbands)*len(self.unique_cbands)
            nkpoints = len(bands_kpoints)
            plot_energies = np.zeros([nkpoints,eh_size])
            plot_weights = np.zeros([nkpoints,eh_size])
            for eh,(v,c) in enumerate(product(self.unique_vbands,self.unique_cbands)):
                plot_energies[:,eh] = energies[:,c]-energies[:,v]
                plot_weights[:,eh] = weights[:,c] 

        if f: plot_weights = f(plot_weights)
        size *= 1.0/np.max(plot_weights)
        ybs = YambopyBandStructure(plot_energies, bands_kpoints, weights=plot_weights, kpath=path, size=size)
        return ybs

    def plot_exciton_bs_ax(self,ax,energies_db,path,excitons,size=1,space='bands',f=None,debug=None):
        ybs = self.get_exciton_bs(energies_db,path,excitons,size=size,space=space,f=f,debug=debug)
        return ybs.plot_ax(ax) 

    @add_fig_kwargs
    def plot_exciton_bs(self,energies_db,path,excitons,size=1,space='bands',f=None,debug=False,**kwargs):
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        self.plot_exciton_bs_ax(ax,energies_db,path,excitons,size=size,space=space,f=f,debug=debug)
        return fig

    def interpolate(self,energies,path,excitons,lpratio=5,f=None,size=1,verbose=True,**kwargs):
        """ Interpolate exciton bandstructure using SKW interpolation from Abipy
            This function is still with some bugs...
            for instance, kpoints_indexes should be read from savedb

        """

        if verbose:
            print("This interpolation is provided by the SKW interpolator implemented in Abipy")

        #print(energies)
        #print(k2)
        #print(len(k2))

        lattice = self.lattice
        cell = (lattice.lat, lattice.red_atomic_positions, lattice.atomic_numbers)
        nelect = 0
        # Here there is something strange...

        fermie = kwargs.pop('fermie',0)
        ##
        symrel = [sym for sym,trev in zip(lattice.sym_rec_red,lattice.time_rev_list) if trev==False ]
        time_rev = True
 
        weights = self.get_exciton_weights(excitons)
        #print('weights.shape')
        #print(weights.shape)
        weights = weights[:,self.start_band:self.mband]
        #print('self.start_band')
        #print(self.start_band)
        #print('self.mband')
        #print(self.mband)
        if f: weights = f(weights)
        size *= 1.0/np.max(weights)
        ibz_nkpoints = max(lattice.kpoints_indexes)+1
        #k1,k2,k3=self.lattice.expand_kpts()
        #print(len(k2))
        #print(len(lattice.kpoints_indexes))
        #print(ibz_nkpoints)
        #print('error')
        kpoints = lattice.red_kpoints

        #map from bz -> ibz:
        # bug here? it is self.mband, but why?
        ibz_weights = np.zeros([ibz_nkpoints,self.mband-self.start_band]) 
        #print('ibz_weights.shape')
        #print(ibz_weights.shape)
        #print(self.nbands)
        #print(self.mband)
        #print(self.start_band)
        #exit()
        ibz_kpoints = np.zeros([ibz_nkpoints,3])
        #print(ibz_weights.shape)
        #print(weights.shape)
        #print(lattice.kpoints_indexes)
        #print('just before error')
        #print(len(lattice.kpoints_indexes))

        # Kpoints indexes must be read from a SAVEDB Class
        k1,k2,k3 = energies.expand_kpts()
        kpoints_indexes = k2

        # Fijar este error
        for idx_bz,idx_ibz in enumerate(kpoints_indexes):
            ibz_weights[idx_ibz,:] = weights[idx_bz,:] 
            ibz_kpoints[idx_ibz] = lattice.red_kpoints[idx_bz]
        #get eigenvalues along the path
        if isinstance(energies,(YamboSaveDB,YamboElectronsDB)):
            ibz_energies = energies.eigenvalues[0,:,self.start_band:self.mband]
        elif isinstance(energies,YamboQPDB):
            ibz_energies = energies.eigenvalues_qp # to be done for spin-UP channel
        else:
            raise ValueError("Energies argument must be an instance of YamboSaveDB,"
                             "YamboElectronsDB or YamboQPDB. Got %s"%(type(energies)))

        #interpolate energies
        na = np.newaxis
        skw = SkwInterpolator(lpratio,ibz_kpoints,ibz_energies[na,:,:],fermie,nelect,cell,symrel,time_rev,verbose=verbose)
        kpoints_path = path.get_klist()[:,:3]
        energies = skw.interp_kpts(kpoints_path).eigens

        #interpolate weights
        na = np.newaxis
        skw = SkwInterpolator(lpratio,ibz_kpoints,ibz_weights[na,:,:],fermie,nelect,cell,symrel,time_rev,verbose=verbose)
        kpoints_path = path.get_klist()[:,:3]
        exc_weights = skw.interp_kpts(kpoints_path).eigens

        #create band-structure object
        exc_bands = YambopyBandStructure(energies[0],kpoints_path,kpath=path,weights=exc_weights[0],size=size,**kwargs)
        #shift top v_band to zero
        exc_bands.set_fermi(self.nvbands)

        return exc_bands

    def interpolate_transitions(self,energies,path,excitons,lpratio=5,f=None,size=1,verbose=True,**kwargs):
        """ Interpolate exciton bandstructure using SKW interpolation from Abipy
        """

        if verbose:
            print("This interpolation is provided by the SKW interpolator implemented in Abipy")

        lattice = self.lattice
        cell = (lattice.lat, lattice.red_atomic_positions, lattice.atomic_numbers)
        nelect = 0
        # Here there is something strange...
        fermie = kwargs.pop('fermie',0)
        ##
        symrel = [sym for sym,trev in zip(lattice.sym_rec_red,lattice.time_rev_list) if trev==False ]
        time_rev = True

        #vmin, vmax = self.unique_vbands[0], self.unique_vbands[1]
        #cmin, cmax = self.unique_cbands[0], self.unique_cbands[1]

        transitions = self.get_exciton_transitions(excitons)
        transitions = transitions[:,:,:]

        ibz_nkpoints = max(lattice.kpoints_indexes)+1
        kpoints = lattice.red_kpoints

        #map from bz -> ibz:
        ibz_transitions = np.zeros([ibz_nkpoints,self.nvbands,self.ncbands])
        ibz_kpoints = np.zeros([ibz_nkpoints,3])
        for idx_bz,idx_ibz in enumerate(lattice.kpoints_indexes):
            ibz_transitions[idx_ibz,:,:] = transitions[idx_bz,:,:] 
            ibz_kpoints[idx_ibz] = lattice.red_kpoints[idx_bz]

        #get eigenvalues along the path
        if isinstance(energies,(YamboSaveDB,YamboElectronsDB)):
            ibz_energies = energies.eigenvalues[:,self.start_band:self.mband]
        elif isinstance(energies,YamboQPDB):
            ibz_energies = energies.eigenvalues_qp
        else:
            raise ValueError("Energies argument must be an instance of YamboSaveDB,"
                             "YamboElectronsDB or YamboQPDB. Got %s"%(type(energies)))

        #interpolate energies
        na = np.newaxis
        skw = SkwInterpolator(lpratio,ibz_kpoints,ibz_energies[na,:,:],fermie,nelect,cell,symrel,time_rev,verbose=verbose)
        kpoints_path = path.get_klist()[:,:3]
        energies = skw.interp_kpts(kpoints_path).eigens
     
        #interpolate transitions
        na = np.newaxis
        skw = SkwInterpolator(lpratio,ibz_kpoints,ibz_transitions[na,:,:],fermie,nelect,cell,symrel,time_rev,verbose=verbose)
        kpoints_path = path.get_klist()[:,:3]
        exc_transitions = skw.interp_kpts(kpoints_path).eigens

        print(exc_transitions.shape)
        exit()

        #create band-structure object
        exc_bands = YambopyBandStructure(energies[0],kpoints_path,kpath=path,weights=exc_weights[0],size=size,**kwargs)
        exc_bands.set_fermi(self.nvbands)

        return exc_transitions

    def interpolate_spin(self,energies,spin_proj,path,excitons,lpratio=5,f=None,size=1,verbose=True,**kwargs):
        """ Interpolate exciton bandstructure using SKW interpolation from Abipy
        """

        if verbose:
            print("This interpolation is provided by the SKW interpolator implemented in Abipy")

        lattice = self.lattice
        cell = (lattice.lat, lattice.red_atomic_positions, lattice.atomic_numbers)
        nelect = 0
        # Here there is something strange...

        fermie = kwargs.pop('fermie',0)
        ##
        symrel = [sym for sym,trev in zip(lattice.sym_rec_red,lattice.time_rev_list) if trev==False ]
        time_rev = True
 
        weights = self.get_exciton_weights(excitons)
        weights = weights[:,self.start_band:self.mband]
        if f: weights = f(weights)
        size *= 1.0/np.max(weights)
        ibz_nkpoints = max(lattice.kpoints_indexes)+1
        kpoints = lattice.red_kpoints

        #map from bz -> ibz:
        print("ibz_nkpoints")
        print(ibz_nkpoints)
        print("weights.shape")
        print(weights.shape)
        print(self.unique_vbands)
        print(self.unique_cbands)
        v_1 = self.unique_vbands[ 0]
        v_2 = self.unique_cbands[-1] + 1
        #exit()
        ibz_weights = np.zeros([ibz_nkpoints,self.nbands])
        ibz_kpoints = np.zeros([ibz_nkpoints,3])
        ibz_spin    = np.zeros([ibz_nkpoints,self.nbands])
        for idx_bz,idx_ibz in enumerate(lattice.kpoints_indexes):
            ibz_weights[idx_ibz,:] = weights[idx_bz,:] 
            ibz_kpoints[idx_ibz]   = lattice.red_kpoints[idx_bz]
            ibz_spin[idx_ibz,:]    = spin_proj[idx_bz,v_1:v_2]
        #get eigenvalues along the path
        if isinstance(energies,(YamboSaveDB,YamboElectronsDB)):
            ibz_energies = energies.eigenvalues[:,self.start_band:self.mband]
        elif isinstance(energies,YamboQPDB):
            ibz_energies = energies.eigenvalues_qp
        else:
            raise ValueError("Energies argument must be an instance of YamboSaveDB,"
                             "YamboElectronsDB or YamboQPDB. Got %s"%(type(energies)))

        #interpolate energies
        na = np.newaxis
        print("na")
        print(na)
        skw = SkwInterpolator(lpratio,ibz_kpoints,ibz_energies[na,:,:],fermie,nelect,cell,symrel,time_rev,verbose=verbose)
        kpoints_path = path.get_klist()[:,:3]
        energies = skw.interp_kpts(kpoints_path).eigens
     
        #interpolate weights
        na = np.newaxis
        skw = SkwInterpolator(lpratio,ibz_kpoints,ibz_weights[na,:,:],fermie,nelect,cell,symrel,time_rev,verbose=verbose)
        kpoints_path = path.get_klist()[:,:3]
        exc_weights = skw.interp_kpts(kpoints_path).eigens

        #interpolate spin projection
        na = np.newaxis
        print("na")
        print(na)
        skw = SkwInterpolator(lpratio,ibz_kpoints,ibz_spin[na,:,:],fermie,nelect,cell,symrel,time_rev,verbose=verbose)
        kpoints_path = path.get_klist()[:,:3]
        spin_inter   = skw.interp_kpts(kpoints_path).eigens
        print("spin_inter")
        print(spin_inter)

        #create band-structure object
        exc_bands = YambopyBandStructure(energies[0],kpoints_path,kpath=path,weights=exc_weights[0],spin_proj=spin_inter[0],size=size,**kwargs)
        exc_bands.set_fermi(self.nvbands)

        return exc_bands
 
    def get_amplitudes_phases(self,excitons=(0,),repx=list(range(1)),repy=list(range(1)),repz=list(range(1))):
        """ get the excitonic amplitudes and phases
        """
        if self.eigenvectors is None:
            raise ValueError('This database does not contain Excitonic states,'
                             'please re-run the yambo BSE calculation with the WRbsWF option in the input file.')
        if isinstance(excitons, int):
            excitons = (excitons,)
       
        car_kpoints = self.lattice.car_kpoints
        nkpoints = len(car_kpoints)
        print(nkpoints)
        amplitudes = np.zeros([nkpoints])
        phases     = np.zeros([nkpoints],dtype=np.complex64)
        for exciton in excitons:
            #the the eigenstate
            eivec = self.eigenvectors[exciton-1]
           
            total = 0
            for eh,kvc in enumerate(self.table):
                ikbz, v, c = kvc-1
                Acvk = eivec[eh]
                phases[ikbz]     += Acvk
                amplitudes[ikbz] += np.abs(Acvk)

        #replicate kmesh
        red_kmesh,kindx = replicate_red_kmesh(self.lattice.red_kpoints,repx=repx,repy=repy,repz=repz)
        car_kpoints = red_car(red_kmesh,self.lattice.rlat)

        return car_kpoints, amplitudes[kindx], np.angle(phases)[kindx]

    def get_chi(self,dipoles=False,dir=[1,1,0],emin=0,emax=10,estep=0.01,broad=0.1,q0norm=1e-5, nexcitons='all',spin_degen=2,verbose=0,**kwargs):
        """
        Calculate the dielectric response function using excitonic states
        """
        if nexcitons == 'all': nexcitons = self.nexcitons

        if 'gauge' in kwargs:
            gauge = kwargs['gauge']
        else:
            gauge = 'length'

        #energy range
        w = np.arange(emin,emax,estep,dtype=np.float32)
        nenergies = len(w)
        
        if verbose:
            print("energy range: %lf -> +%lf -> %lf "%(emin,estep,emax))
            print("energy steps: %lf"%nenergies)

        #initialize the susceptibility intensity
        chi = np.zeros([len(w)],dtype=np.complex64)

        if not dipoles:
            #get dipole
            EL1 = self.l_residual
            EL2 = self.r_residual
        else:
            #calculate exciton-light coupling
            if verbose: print("calculate exciton-light coupling")
            if 'folder' in kwargs:
                folder = kwargs['folder']
            if 'filename' in kwargs:
                filename = kwargs['filename']
            if 'statelist' in kwargs:
                statelist = kwargs['statelist']  
                nexcitons = len(statelist)  
            try:    
                EL2_vec = self.get_exc_dipoles_from_file(folder,filename, statelist, gauge)
            except:
                EL2_vec = self.get_exc_dipoles(folder,filename, statelist, gauge)
            EL2 = np.dot(EL2_vec,dir)
            EL1 = np.conj(EL2) 

        if isinstance(broad,float): broad = [broad]*nexcitons

        if isinstance(broad,tuple): 
            broad_slope = broad[1]-broad[0]
            min_exciton = np.min(self.eigenvalues.real)
            broad = [ broad[0]+(es-min_exciton)*broad_slope for es in self.eigenvalues[:nexcitons].real]

        if "gaussian" in broad or "lorentzian" in broad:
            i = broad.find(":")
            if i != -1:
                value, eunit = broad[i+1:].split()
                if eunit == "eV": sigma = float(value)
                else: raise ValueError('Unknown unit %s'%eunit)

            f = gaussian if "gaussian" in broad else lorentzian
            broad = np.zeros([nexcitons])
            for s in range(nexcitons):
                es = self.eigenvalues[s].real
                broad += f(self.eigenvalues.real,es,sigma)
            broad = 0.1*broad/nexcitons

        #iterate over the excitonic states
        for s in range(nexcitons):
            #get exciton energy
            es = self.eigenvalues[s]
 
            #calculate the green's functions
            G1 = -1/(   w - es + broad[s]*I)
            G2 = -1/( - w - es - broad[s]*I)

            r = EL1[s]*EL2[s]
            chi += r*G1 + r*G2

        #dimensional factors
        try:
            if not self.Qpt=='1': q0norm = 2*np.pi*np.linalg.norm(self.car_qpoint)
        except:
            print("[WARNING] 1/q^2 set to 1 in eps2")
            q0norm=1
        try:
            if self.q_cutoff is not None: q0norm = self.q_cutoff
        except:
            print("[WARNING] 1/q^2 set to 1 in eps2")
            q0norm=1

        d3k_factor = self.lattice.rlat_vol/self.lattice.nkpoints
        cofactor = ha2ev*spin_degen/(2*np.pi)**3 * d3k_factor * (4*np.pi)  / q0norm**2
        if (gauge == 'velocity') : cofactor = cofactor*q0norm**2
        chi = 1. + chi*cofactor #We are actually computing the epsilon, not the chi.

        return w,chi

    def get_abs_lifetimes(self,dipoles=None,dir=0,emin=0,emax=10,estep=0.01,broad=0.1,q0norm=1e-5, nexcitons='all',spin_degen=2,verbose=0,**kwargs):
        """
        Calculate the dielectric response function using excitonic states
        """
        if 'gauge' in kwargs:
            gauge = kwargs['gauge']
        else:
            gauge = 'length'
        if 'pl_res' in kwargs: 
            pl_res = kwargs['pl_res']
        else:
            pl_res = False
        if 'degen_step' in kwargs: 
            degen_step = kwargs['degen_step']
        else:
            degen_step = 0.001
        if nexcitons == 'all': nexcitons = self.nexcitons
        
        #energy range
        w = np.arange(emin,emax,estep,dtype=np.float32)
        nenergies = len(w)
        
        if verbose:
            print("energy range: %lf -> +%lf -> %lf "%(emin,estep,emax))
            print("energy steps: %lf"%nenergies)

        #initialize the susceptibility intensity
        abs = np.zeros([len(w)],dtype=np.complex64)

        if dipoles is None:
            #get dipole
            EL1 = self.l_residual
            EL2 = self.r_residual
        else:
            #calculate exciton-light coupling
            if verbose: print("calculate exciton-light coupling")
            EL1,EL2 = self.project1(dipoles.dipoles[:,dir],nexcitons) 

        a1 = self.lattice.lat[0]
        a2 = self.lattice.lat[1]
        A = np.linalg.norm(np.cross(a1,a2))
        # q0 norm factor 
        q0_norm = 1e-5
        tmp_factor = 4*np.pi/ha2ev/(q0_norm**2*self.nkpoints*A*speed_of_light)
        #iterate over the excitonic states
        for s in range(nexcitons):
            #get exciton energy
            es = self.eigenvalues[s]
 
            # calculate absorption
            r = EL1[s]*EL2[s]
#            if (gauge=='length'):
            gg = tmp_factor*es*r
            broad_s = hbar*gg.real/AU2S
            # elif (gauge =='velocity'):
            #     gg = 4.*np.pi*es/ha2ev*(r/(es/ha2ev)**2*self.nkpoints)/(A*speed_of_light)
            #     broad_s = hbar*gg.real/AU2S                
            abs += 1j/(np.pi*self.nexcitons*self.nkpoints)*(broad_s*r)/((w-es)**2-broad_s**2)
        #dimensional factors
        # try:
        #     if not self.Qpt=='1': q0norm = 2*np.pi*np.linalg.norm(self.car_qpoint)
        # except:
        #     print("[WARNING] 1/q^2 set to 1 in eps2")
        #     q0norm=1
        # try:
        #     if self.q_cutoff is not None: q0norm = self.q_cutoff
        # except:
        #     print("[WARNING] 1/q^2 set to 1 in eps2")
        #     q0norm=1

        # d3k_factor = self.lattice.rlat_vol/self.lattice.nkpoints
        # cofactor = ha2ev*spin_degen/(2*np.pi)**3 * d3k_factor * (4*np.pi)  / q0norm**2
        
        # abs = 1. + abs*cofactor #We are actually computing the epsilon, not the chi.

        return w,abs    

    def get_absorbance(self,deltat,dipoles=None,dir=0,emin=0,emax=10,estep=0.01,broad=0.1,q0norm=1e-5, nexcitons='all',spin_degen=2,verbose=0,**kwargs):
        """
        Calculate absorbance from polarizability 2D
            Arguments:
            deltat ->  thickness of 2D material in Bohr
        """        
        if 'gauge' in kwargs:
            gauge = kwargs['gauge']
        else:
            gauge = 'length'

        w, chi = self.get_chi(dipoles=None,dir=0,emin=0,emax=10,estep=0.01,broad=0.1,q0norm=1e-5, nexcitons='all',spin_degen=2,verbose=0,gauge='length')        
        absorbance = np.zeros_like(chi)

        absorbance = w/speed_of_light*deltat + w/speed_of_light*4*np.pi*chi
    
        return absorbance

    def plot_chi_ax(self,ax,reim='im',n_brightest=-1,**kwargs):
        """Plot chi on a matplotlib axes"""
        w,chi = self.get_chi(**kwargs)
        #cleanup kwargs variables
        cleanup_vars = ['dipoles','dir','emin','emax','estep','broad',
                        'q0norm','nexcitons','spin_degen','verbose']
        for var in cleanup_vars: kwargs.pop(var,None)
        if 're' in reim: ax.plot(w,chi.real,**kwargs)
        if 'im' in reim: ax.plot(w,chi.imag,**kwargs)
        ax.set_ylabel('$Im(\chi(\omega))$')
        ax.set_xlabel('Energy (eV)')
        #plot vertical bar on the brightest excitons
        if n_brightest>-1:
            exc_e,exc_i = self.get_sorted()
            for i,idx in exc_i[:n_brightest]:
                exciton_energy,idx = exc_e[idx]
                ax.axvline(exciton_energy,c='k')
                ax.text(exciton_energy,0.1,idx,rotation=90)
        return w,chi

    @add_fig_kwargs
    def plot_chi(self,n_brightest=-1,**kwargs):
        """Produce a figure with chi"""
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        self.plot_chi_ax(ax,n_brightest=n_brightest,**kwargs)
        return fig

    def save_chi(self,filename,**kwargs):
        """Compute chi and dump it to file"""
        w,chi = self.get_chi(**kwargs)
        np.savetxt(filename,np.array([w,chi.imag,chi.real]).T)

    ##########################################
    #  SPIN DEPENDENT PART UNDER DEVELOPMENT #
    ##########################################

    def exciton_bs_spin_pol(self,energies,path,excitons=(0,),debug=False):
        """
        Calculate exciton band-structure
        This function should contains the case of non-polarized calculations.
        Now is a first version
            
            Arguments:
            energies -> can be an instance of YamboSaveDB or YamboQBDB
            path     -> path in reduced coordinates in which to plot the band structure
            exciton  -> exciton index to plot
            spin     -> ??
        """
        if self.eigenvectors is None:
            raise ValueError('This database does not contain Excitonic states,'
                              'please re-run the yambo BSE calculation with the WRbsWF option in the input file.')
        if isinstance(excitons, int):
            excitons = (excitons,)
        #get full kmesh
        kpoints = self.lattice.red_kpoints
        path = np.array(path)

        rep = list(range(-1,2))
        kpoints_rep, kpoints_idx_rep = replicate_red_kmesh(kpoints,repx=rep,repy=rep,repz=rep)
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
        if isinstance(energies,(YamboSaveDB,YamboElectronsDB)):
            #expand eigenvalues to the full brillouin zone
            energies_up = energies.eigenvalues[0,self.lattice.kpoints_indexes]
            energies_dw = energies.eigenvalues[1,self.lattice.kpoints_indexes]
            
        elif isinstance(energies,YamboQPDB):
            #expand the quasiparticle energies to the bull brillouin zone
            # Check this is OK
            # To-do
            print('todo')
            print(energies.eigenvalues_qp[self.lattice.kpoints_indexes].shape)
            pad_energies_up = energies.eigenvalues_qp[self.lattice.kpoints_indexes,:,0]
            pad_energies_dw = energies.eigenvalues_qp[self.lattice.kpoints_indexes,:,1]
            print(energies.max_band)
            min_band = energies.min_band       # check this
            print(pad_energies_up.shape)
            print(pad_energies_dw.shape)
            nkpoints, nbands = pad_energies_up.shape
            #exit()
            #pad_energies = energies.eigenvalues_qp[self.lattice.kpoints_indexes]
            #min_band = energies.min_band       # check this
            #nkpoints, nbands = pad_energies.shape
            energies_up = np.zeros([nkpoints,energies.max_band])
            energies_dw = np.zeros([nkpoints,energies.max_band])
            energies_up[:,min_band-1:] = pad_energies_up 
            energies_dw[:,min_band-1:] = pad_energies_dw
        else:
            raise ValueError("Energies argument must be an instance of YamboSaveDB,"
                             "YamboElectronsDB or YamboQPDB. Got %s"%(type(energies)))

        energies_up, energies_dw = energies_up[band_indexes],energies_dw[band_indexes]

        weights_up, weights_dw = self.get_exciton_weights_spin_pol(excitons)      
        weights_up, weights_dw = weights_up[band_indexes], weights_dw[band_indexes]

        #make top valence band to be zero
        fermi_level = max([ max(energies_up[:,max(self.unique_vbands)]), max(energies_dw[:,max(self.unique_vbands)] ) ])
        energies_up -= fermi_level  
        energies_dw -= fermi_level  
        
        return np.array(band_kpoints), energies_up, energies_dw, weights_up, weights_dw

    def get_exciton_bs_spin_pol(self,energies_db,path,excitons,size_up=1,size_dw=1,space='bands',f=None,debug=False):
        """
        Get a YambopyBandstructure object with the exciton band-structure SPIN POLARIZED
        
            Arguments:
            ax          -> axis extance of matplotlib to add the plot to
            lattice     -> Lattice database
            energies_db -> Energies database, can be either a SaveDB or QPDB
            path        -> Path in the brillouin zone
        """
        from qepy.lattice import Path
        if not isinstance(path,Path): 
            raise ValueError('Path argument must be a instance of Path. Got %s instead'%type(path))
        if space == 'bands':
            if self.spin_pol=='pol':
               bands_kpoints, energies_up, energies_dw, weights_up, weights_dw = self.exciton_bs_spin_pol(energies_db, path.kpoints, excitons, debug)
               nkpoints = len(bands_kpoints)
               plot_energies_up = energies_up[:,self.start_band:self.mband]
               plot_energies_dw = energies_dw[:,self.start_band:self.mband]
               plot_weights_up  = weights_up[:,self.start_band:self.mband]
               plot_weights_dw  = weights_dw[:,self.start_band:self.mband]
        #    elif spin_pol=='pol':
               
        else:
            raise NotImplementedError('TODO')
            eh_size = len(self.unique_vbands)*len(self.unique_cbands)
            nkpoints = len(bands_kpoints)
            plot_energies = np.zeros([nkpoints,eh_size])
            plot_weights = np.zeros([nkpoints,eh_size])
            for eh,(v,c) in enumerate(product(self.unique_vbands,self.unique_cbands)):
                plot_energies[:,eh] = energies[:,c]-energies[:,v]
                plot_weights[:,eh] = weights[:,c] 

        
        if f: plot_weights_up, plot_weights_dw = f(plot_weights_up), f(plot_weights_dw)
        size_plot_up = 100.0 # 1.0/np.max(plot_weights_up)*100.0
        size_plot_dw = 100.0 # 1.0/np.max(plot_weights_dw)*100.0
        ybs_up = YambopyBandStructure(plot_energies_up, bands_kpoints, weights=plot_weights_up, kpath=path, size=size_plot_up)
        ybs_dw = YambopyBandStructure(plot_energies_dw, bands_kpoints, weights=plot_weights_dw, kpath=path, size=size_plot_dw)
        
        #from numpy import arange
        #x = arange(nkpoints)
        #import matplotlib.pyplot as plt
        #for ib1 in range(17):
        #    plt.plot(x,plot_energies_up[:,ib1],'r--')
        #    plt.scatter(x,plot_energies_up[:,ib1],s=plot_weights_up[:,ib1]*size_up*1000,c='red')
        #for ib1 in range(17):
        #    plt.plot(x,plot_energies_dw[:,ib1],'b--')
        #    plt.scatter(x,plot_energies_dw[:,ib1],s=plot_weights_dw[:,ib1]*size_dw*1000,c='blue')
        #plt.show()
        #exit()

        return ybs_up, ybs_dw

    def get_exciton_weights_spin_pol(self,excitons):
    
        """get weight of state in each band for spin-polarized case"""
        table_up, table_dw, table_updw = [] , [], []
        for t,kcv in enumerate(self.table):
            k,c,v,c_s,v_s = kcv-1   # We substract 1 to be consistent with python numbering of arrays
            if c_s == 0 and v_s == 0:
               table_up.append(np.array(kcv[0:3]))
            if c_s == 1 and v_s == 1:
               table_dw.append(np.array(kcv[0:3]))
            if c_s == 1 and v_s == 0:
               table_updw.append(np.array(kcv[0:3]))
        table_up=np.array(table_up)
        table_dw=np.array(table_dw)
        table_updw=np.array(table_updw)

        self.unique_vbands_up = np.unique(table_up[:,1]-1)
        self.unique_cbands_up = np.unique(table_up[:,2]-1)
        self.unique_vbands_dw = np.unique(table_dw[:,1]-1)
        self.unique_cbands_dw = np.unique(table_dw[:,2]-1)
        self.mband_up = max(self.unique_cbands_up) + 1
        self.mband_dw = max(self.unique_cbands_dw) + 1
        self.start_band_up = min(self.unique_vbands_up)
        self.start_band_dw = min(self.unique_vbands_dw)

        weights_up = np.zeros([self.nkpoints,self.mband_up])
        weights_dw = np.zeros([self.nkpoints,self.mband_dw])
        
        for exciton in excitons:
            #get the eigenstate
            eivec = self.eigenvectors[exciton-1]
            #add weights
            sum_weights = 0
            for t,kcv in enumerate(self.table):
                k,c,v,c_s,v_s = kcv-1   # We substract 1 to be consistent with python numbering of arrays
                this_weight = abs2(eivec[t])
                if c_s == 0 and v_s == 0:
                   weights_up[k,c] += this_weight
                   weights_up[k,v] += this_weight
                elif c_s == 1 and v_s == 1: 
                   weights_dw[k,c] += this_weight
                   weights_dw[k,v] += this_weight
                sum_weights += this_weight
            if abs(sum_weights - 1) > 1e-3: raise ValueError('Excitonic weights does not sum to 1 but to %lf.'%sum_weights)
 
        return weights_up, weights_dw

    def interpolate_spin_pol(self,energies,path,excitons,lpratio=5,f=None,size_up=1.0,size_dw=1.0,verbose=True,**kwargs):
        """ Interpolate exciton bandstructure using SKW interpolation from
        Abipy and SPIN-POLARIZED CALCULATIONS
        """

        if verbose:
            print("This interpolation is provided by the SKW interpolator implemented in Abipy")

        lattice = self.lattice
        cell = (lattice.lat, lattice.red_atomic_positions, lattice.atomic_numbers)
        nelect = 0
        # Here there is something strange...

        fermie = kwargs.pop('fermie',0)
        ##
        symrel = [sym for sym,trev in zip(lattice.sym_rec_red,lattice.time_rev_list) if trev==False ]
        time_rev = True
 
        weights_up, weights_dw = self.get_exciton_weights_spin_pol(excitons)

        weights_up = weights_up[:,self.start_band_up:self.mband_up]
        weights_dw = weights_dw[:,self.start_band_dw:self.mband_dw]

        if f: weights_up, weights_dw = f(weights_up), f(weights_dw)

        size_up *= 1.0/np.max(weights_up)
        size_dw *= 1.0/np.max(weights_dw)

        ibz_nkpoints = max(lattice.kpoints_indexes)+1
        kpoints = lattice.red_kpoints

        #map from bz -> ibz:
        # bug here? it is self.mband, but why?
        ibz_weights_up = np.zeros([ibz_nkpoints,self.mband_up-self.start_band_up]) 
        ibz_weights_dw = np.zeros([ibz_nkpoints,self.mband_dw-self.start_band_dw]) 
        
        ibz_kpoints = np.zeros([ibz_nkpoints,3])
        print(self.mband_up,self.start_band_up)
        print(ibz_weights_up.shape)
        print(weights_up.shape)
        print(lattice.kpoints_indexes)
        print('just before error')
        for idx_bz,idx_ibz in enumerate(lattice.kpoints_indexes):
            print(weights_up[idx_bz,:])
            ibz_weights_up[idx_ibz,:], ibz_weights_dw[idx_ibz,:]= weights_up[idx_bz,:], weights_dw[idx_bz,:] 
            ibz_kpoints[idx_ibz] = lattice.red_kpoints[idx_bz]

        #get eigenvalues along the path
        # DFT values from SAVE
        if isinstance(energies,(YamboSaveDB,YamboElectronsDB)):
            ibz_energies_up = energies.eigenvalues[0,:,self.start_band:self.mband] # spin-up channel
            ibz_energies_dw = energies.eigenvalues[1,:,self.start_band:self.mband] # spin-dw channel
            ibz_kpoints_qp  = ibz_kpoints
        # GW values from ndb.QP
        elif isinstance(energies,YamboQPDB):
            ibz_nkpoints_gw=len(energies.kpoints_iku)
            if not ibz_nkpoints == ibz_nkpoints_gw :
               print('GW and BSE k-grid are differents!')
               kpoints_gw_iku = energies.kpoints_iku
               kpoints_gw_car = np.array([ k/lattice.alat for k in kpoints_gw_iku ])
               kpoints_gw_red = car_red( kpoints_gw_car,lattice.rlat)
               ibz_kpoints_qp = kpoints_gw_red
            else:
               ibz_kpoints_qp = ibz_kpoints
            pad_energies_up = energies.eigenvalues_qp[:,:,0]
            pad_energies_dw = energies.eigenvalues_qp[:,:,1]
            min_band = energies.min_band
            nkpoints, nbands = pad_energies_up.shape

            ibz_energies_up = pad_energies_up 
            ibz_energies_dw = pad_energies_dw
            #print('ibz',ibz_energies_up.shape)
        else:
            raise ValueError("Energies argument must be an instance of YamboSaveDB,"
                             "YamboElectronsDB or YamboQPDB. Got %s"%(type(energies)))

        #interpolate energies
        na = np.newaxis

        skw_up = SkwInterpolator(lpratio,ibz_kpoints_qp,ibz_energies_up[na,:,:],fermie,nelect,cell,symrel,time_rev,verbose=verbose)
        skw_dw = SkwInterpolator(lpratio,ibz_kpoints_qp,ibz_energies_dw[na,:,:],fermie,nelect,cell,symrel,time_rev,verbose=verbose)
        kpoints_path = path.get_klist()[:,:3]
        energies_up = skw_up.interp_kpts(kpoints_path).eigens
        energies_dw = skw_dw.interp_kpts(kpoints_path).eigens
     
        #interpolate weights
        na = np.newaxis
        skw_up = SkwInterpolator(lpratio,ibz_kpoints,ibz_weights_up[na,:,:],fermie,nelect,cell,symrel,time_rev,verbose=verbose)
        skw_dw = SkwInterpolator(lpratio,ibz_kpoints,ibz_weights_dw[na,:,:],fermie,nelect,cell,symrel,time_rev,verbose=verbose)
        kpoints_path = path.get_klist()[:,:3]
        exc_weights_up = skw_up.interp_kpts(kpoints_path).eigens
        exc_weights_dw = skw_dw.interp_kpts(kpoints_path).eigens

        # Find and set the up-dw Fermi energy to zero
        self.nvbands_up = len(self.unique_vbands_up)
        self.nvbands_dw = len(self.unique_vbands_dw)
        fermi_up_dw = max([max(energies_up[0][:,self.nvbands_up-1]), max(energies_dw[0][:,self.nvbands_dw-1])])

        #create band-structure object
        exc_bands_up = YambopyBandStructure(energies_up[0],kpoints_path,kpath=path,fermie=fermi_up_dw,weights=exc_weights_up[0],size=size_up,**kwargs)
        exc_bands_dw = YambopyBandStructure(energies_dw[0],kpoints_path,kpath=path,fermie=fermi_up_dw,weights=exc_weights_dw[0],size=size_dw,**kwargs)

        return exc_bands_up, exc_bands_dw

    # Riccardo Reho methods
    def eigenvectors_transitions(self):
        '''
        Get the BSE eigenvectors in lambda,v,c,k space as [t, nk, nv, nc, i_s_v, i_s_v]
        For now, this utility is not used and not tested.
        '''
        if (self.spin_vars == "no"):
            Alambdavck = np.zeros((self.ntransitions, self.nkpoints, self.nbandsv, self.nbandsc,1,1), dtype = complex)
        elif(self.spin_vars == "yes"):
            # Not sure this is the case for spin
            Alambdavck = np.zeros((self.ntransitions, self.nbandsv, self.nbandsc,2,2), dtype = complex)
        # Loop to create the array over indices
        for t in range(0, self.nexcitons):
            for tp in range(0,self.nexcitons):
                ik = self.table[tp][0]
                iv = self.table[tp][1]
                ic = self.table[tp][2]
                isc = self.table[tp][3]
                isv = self.table[tp][4]

                Alambdavck[t, ik,iv, ic, isv, isc] = self.eigenvectors[t][tp]
        return Alambdavck
    
    def get_exc_dipoles(self, folder, filename, trange, gauge = 'length'):
        import netCDF4 as ns
        '''
        get exciton dipole according to \sum_{cvk} A^\lambda_cvk *d_cvk
        The dipole matrix in ndb.dipoles has the following indexes (to be checked with Yambo version):
        [nspin, nkpoints, nbands valence, nbands conduction, cartesian direction, 2]
        2 is for complex.        
        '''
        exc_dipoles = np.zeros((self.ntransitions, 3),dtype = complex)
        A_ttp = np.ma.asarray(self.eigenvectors) # A^lambda_cvk in transiton space
        # Get dipoles from ndb.dipoles
        
        path_filename = os.path.join(folder,filename)
        if not os.path.isfile(path_filename):
            raise FileNotFoundError("File %s not found in YamboExcitonDB"%path_filename)        
        dip_db = ns.Dataset(f'{path_filename}')
        if (gauge == 'length'):
            dip_db = YamboDipolesDB(self.lattice, folder, filename=filename, dip_type='iR',expand=True )
        if (gauge == 'velocity'):
            dip_db = YamboDipolesDB(self.lattice, folder, filename=filename, dip_type='P', expand=True )
        dkcv = dip_db.dipoles
        nbandsv = dkcv.shape[3]
        #dkvc = dipoles[:,:,:,:,0] + 1j*dipoles[:,:,:,:,1]
        for it, t in enumerate(trange):
            for tp in np.arange(0,self.nexcitons):
                # I need to to -1 for Python counting
                ik = self.table[tp][0]-1
                iv = self.table[tp][1]-1
                ic = self.table[tp][2]-1
                isc = self.table[tp][3]-1
                isv = self.table[tp][4]-1
                #ikibz = self.lattice.kpoints_indexes[ik]
                exc_dipoles[t,:] += A_ttp[t,tp]*dkcv[ik,:, ic,iv]
        return exc_dipoles  

    def get_exc_dipoles_from_file(self, folder, filename, trange, gauge = 'length'):
        import netCDF4 as ns
        '''
        get exciton dipole according to \sum_{cvk} A^\lambda_cvk *d_cvk
        The dipole matrix in ndb.dipoles has the following indexes (to be checked with Yambo version):
        [nspin, nkpoints, nbands valence, nbands conduction, cartesian direction, 2]
        2 is for complex.        
        '''
        exc_dipoles = np.zeros((self.ntransitions, 3),dtype = complex)
        A_ttp = np.ma.asarray(self.eigenvectors) # A^lambda_cvk in transiton space
        # Get dipoles from ndb.dipoles
        
        path_filename = os.path.join(folder,filename)
        if not os.path.isfile(path_filename):
            raise FileNotFoundError("File %s not found in YamboExcitonDB"%path_filename)        
        if (gauge == 'length'):
            dkcv = np.load(f'{path_filename}')
        if (gauge == 'velocity'):
            dkcv = np.load(f'{path_filename}')
        #dkvc = dipoles[:,:,:,:,0] + 1j*dipoles[:,:,:,:,1]
        for it, t in enumerate(trange):
            for tp in np.arange(0,self.nexcitons):
                # I need to to -1 for Python counting
                ik = self.table[tp][0]-1
                iv = self.table[tp][1]-1
                ic = self.table[tp][2]-1
                isc = self.table[tp][3]-1
                isv = self.table[tp][4]-1
                #ikibz = self.lattice.kpoints_indexes[ik]
                exc_dipoles[t,:] += A_ttp[t,tp]*dkcv[ik,:, ic,iv]
        return exc_dipoles       
    
    def write_exc_dipoles(self, folder, filename = 'ndb.dipoles', prefix = 'exc_dipole', trange=None):
        """
        Write exc dipole to a file
        """
        if (trange is None ):
            trange = np.arange(0,self.nexcitons)
        #get exc_dipole
        eig = self.eigenvalues.real
        exc_dipoles = self.get_exc_dipoles(folder,filename, trange)
        print(exc_dipoles.shape)
        #write excitons sorted by energy
        with open('%s.dat'%prefix, 'w') as f:
            f.write("# transition index exc dipole[Bohr]  Re(x) Im(x) Re(y) Im(y) Re(z) Im(z)\n")
            for t,dip in enumerate(exc_dipoles):
                f.write(f"{t+1:d} {dip[0]:.12f} {dip[1]:.12f} {dip[2]:.12f}\n") 

    def get_gamma0(self, statelist = None, degen_step = 0.001, gauge = 'length', 
                   pl_res = False):
        """get gamma0  
        """
        if (statelist is None ):
            statelist = np.arange(1,self.nexcitons+1)
        Omega = self.lattice.lat_vol
        a1 = self.lattice.lat[0]
        a2 = self.lattice.lat[1]
        A = np.linalg.norm(np.cross(a1,a2))
        muS2 = 0 
        gamma0 = 0
        tmpE = np.ma.asarray(self.eigenvalues.real)
        if (pl_res):
            gauge = 'velocity' # force velocity gauge for PL
            tmpI = np.ma.asarray((self.pl_l_residual*self.pl_r_residual))
        else:
            tmpI = np.ma.asarray(self.l_residual * self.r_residual)
        excE = sorted(tmpE)
        excI = tmpI

        gamma0_merged = np.empty(len(statelist), dtype='object')
        # q0 norm factor 
        q0_norm = 1e-5

        for l, state in enumerate(statelist):
            state_internal = state 

            # find states within degen_step window from state
            mesk = np.logical_and(excE>=excE[state_internal]-degen_step,excE<=excE[state_internal]+degen_step)
            states_idx = np.where(mesk==True)[0]
            #compute gamma
            for i, st in enumerate(states_idx):
                ES = excE[st]/ha2ev
                if (gauge=='length'):
                    muS2 = excI[st]/(q0_norm**2) # if you inspect the Yambo code you might expect another 1/((2*np.pi)**3) but I think that d3k_factor/((2np.pi)**3) is actually 1/Omega
                elif (gauge == 'velocity'):
                    muS2 = excI[st]/(ES**2)
                gg = 4.*np.pi*ES*(muS2/self.nkpoints)/(A*speed_of_light)
                gamma0 += gg
            #compute tau
            gamma0_merged[l]  = gamma0.real
        return gamma0_merged/(AU2S*1e12) # convert in 1/ps

    # Function to calculate γS^IP(θ, φ) for a range of φ values
    def gamma_S_IP(self,phi, gamma_S_0,p_S, p_Sx, p_Sy):
        tmp_array = -(p_Sx/p_S) * np.sin(phi) + (p_Sy/p_S) * np.cos(phi)
        return gamma_S_0 * np.conj(tmp_array)*tmp_array
    
    # Function to calculate γS^OOP(θ, φ) for a range of φ values
    def gamma_S_OOP(self,phi, gamma_S_0,p_S, p_Sx, p_Sy):
        tmp_array = (p_Sx/p_S) * np.cos(phi) + (p_Sy/p_S) * np.sin(phi)
        return gamma_S_0 * np.conj(tmp_array)*tmp_array
     
    def gamma_S_OOP_fixed_phi(self,theta, phi, gamma_S_0,p_S, p_Sx, p_Sy):
        tmp_array = (p_Sx/p_S) * np.cos(phi) + (p_Sy/p_S) * np.sin(phi)
        return gamma_S_0 * np.cos(theta)**2*np.conj(tmp_array)*tmp_array  

    def plot_polar_pl_IP(self, theta_deg, phi_range_deg=np.linspace(0, 360, 1000), 
                               statelist = None, degen_step = 0.001, gauge = 'length', 
                               pl_res = False, folder='./', filename = './ndb.dipoles'):
        import matplotlib.pyplot as plt

        # Convert degrees to radians
        if (statelist is None ):
            statelist = np.arange(1,self.nexcitons+1)
        theta = np.radians(theta_deg)
        phi = np.radians(phi_range_deg)
        
        # Assuming some function of theta and phi for demonstration
        # Here we simply use r = 1 for demonstration, as the primary interest is in phi
        gamma0 = self.get_gamma0(statelist, degen_step, gauge, pl_res)
        p_S_array = self.get_exc_dipoles(folder, filename, trange=statelist)
        for i_s, s in enumerate(statelist):
            print(f'Polar plot for exciton-state {s+1}')
            p_S_x = p_S_array[i_s][0]
            p_S_y = p_S_array[i_s][1]
            p_S = np.linalg.norm(p_S_array[i_s])
            gamma_values = self.gamma_S_IP(phi, gamma0[i_s], p_S, p_S_x, p_S)
            # Plotting
            # Plotting using fig, ax
            fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(8, 6))
            ax.plot(phi, gamma_values)  # Plotting on polar axes
            ax.set_title(f'Polar plot IP exciton-state {s+1}')


    def plot_polar_pl_OOP(self, theta_deg, phi_range_deg=np.linspace(0, 360, 1000), 
                               statelist = None, degen_step = 0.001, gauge = 'length', 
                               pl_res = False, folder='./', filename = './ndb.dipoles'):
        import matplotlib.pyplot as plt
        
        # Convert degrees to radians
        if (statelist is None ):
            statelist = np.arange(1,self.nexcitons+1)
        theta = np.radians(theta_deg)
        phi = np.radians(phi_range_deg)
        
        # Assuming some function of theta and phi for demonstration
        # Here we simply use r = 1 for demonstration, as the primary interest is in phi
        gamma0 = self.get_gamma0(statelist, degen_step, gauge, pl_res)*np.cos(theta_deg)**2
        p_S_array = self.get_exc_dipoles(folder, filename, trange=statelist)
        for i_s, s in enumerate(statelist):
            print(f'Polar plot  OOP for exciton-state {s+1}')
            p_S_x = p_S_array[i_s][0]
            p_S_y = p_S_array[i_s][1]
            p_S = np.linalg.norm(p_S_array[i_s])
            gamma_values = self.gamma_S_OOP(phi, gamma0[i_s], p_S, p_S_x, p_S)
            # Plotting
            # Plotting using fig, ax
            fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(8, 6))
            ax.plot(phi, gamma_values)  # Plotting on polar axes
            ax.set_title(f'Polar plot OOP exciton-state {s+1}')

    def plot_polar_fixed_phi(self, phi_deg=0.0, theta_range_deg=np.linspace(0, 360, 1000), 
                               statelist = None, degen_step = 0.001, gauge = 'length', 
                               pl_res = False, folder='./', filename = './ndb.dipoles'):
        import matplotlib.pyplot as plt
        # Convert degrees to radians
        if (statelist is None ):
            statelist = np.arange(1,self.nexcitons+1)
        theta = np.radians(theta_range_deg)
        phi = np.radians(phi_deg)
        
        # Assuming some function of theta and phi for demonstration
        # Here we simply use r = 1 for demonstration, as the primary interest is in phi
        gamma0 = self.get_gamma0(statelist, degen_step, gauge, pl_res)
        p_S_array = self.get_exc_dipoles(folder, filename, trange=statelist)
        for i_s, s in enumerate(statelist):
            print(f'Polar plot  OOP fixed phi for exciton-state {s+1}')
            p_S_x = p_S_array[i_s][0]
            p_S_y = p_S_array[i_s][1]
            p_S = np.linalg.norm(p_S_array[i_s])
            gamma_values = self.gamma_S_OOP_fixed_phi(theta, phi, gamma0[i_s], p_S, p_S_x, p_S)
            # Plotting
            # Plotting using fig, ax
            fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(8, 6))
            ax.plot(theta, gamma_values)  # Plotting on polar axes
            ax.set_title(f'Polar plot OOP fixed phi exciton-state {s+1}')

    ##############################################
    #  END SPIN DEPENDENT PART UNDER DEVELOPMENT #
    ##############################################

    def get_string(self,mark="="):
        lines = []; app = lines.append
        app( marquee(self.__class__.__name__,mark=mark) )
        app( "BSE solved at Q:            %s"%self.Qpt )
        app( "number of excitons:         %d"%self.nexcitons )
        if self.table is not None: 
            app( "number of transitions:      %d"%self.ntransitions )
            app( "number of kpoints:          %d"%self.nkpoints  )
            app( "number of valence bands:    %d"%self.nvbands )
            app( "number of conduction bands: %d"%self.ncbands )
        return '\n'.join(lines)
    
    def __str__(self):
        return self.get_string()
