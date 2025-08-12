from yambopy import *
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import os
from yambopy.plot.plotting import BZ_Wigner_Seitz
from yambopy.kpoints import expand_kpoints
from yambopy.units import ha2ev, ev2cm1

def phonon_overlap(mode1,mode2):
    """
    Compute the overlap between phonon eigenvectors e_ka,n (q).
    This can be useful to disentangle phonon modes along a q-direction.

    - See Zacharias and Giustino PRR 2020, Eq. (55)

    k -> atomic index
    a -> Cartesian index
    n -> branch index
    q -> momentum index

    (the eigenvectors are complex)

    :: mode1 is e (q1)
    :: mode2 is e (q2)
    
    Eigenvectors are read from prefix.dynq files (QE)

    The overlap matrix is:

    M_nm(q1,q2) = \\sum_ka e_ka,n (q1) e^*_ka,m (q2)
    """
    
    # Overlap matrix
    M = np.einsum('nka,mka->nm',mode1,np.conj(mode2))

    # Debug (explicit equivalent of einsum)
    #M2 = np.zeros([nmodes,nmodes],dtype=np.complex)
    #for n in range(nmodes):
    #    for m in range(nmodes):
    #        for k in range(natoms):
    #            for a in range(3):
    #                M2[n,m] += mode1[n,k,a]*np.conj(mode2[m,k,a])

    return M    

class YamboRefineElphDB():
    """
    This class:
    
    (i) replaces the phonon energies stored in the ndb.elph_gkkp_expanded
    databases with those resulting from matdyn.x, which include ASR and LOTO
    proper treatment

    (ii) for 2D systems, replaces the values of g_LO(q=0), which is zero in a QE
    calculation, with the average absolute value of the closest q/=0 elements, 
    approximating the correct behaviour (cusp) for dense enough q-sampling

    - Inputs: 
    :: YamboLatticeDB object
    :: folder of ndb.elph_gkkp_expanded databases
    :: path of prefix.bn file (from matdyn.x) including filename
    :: [optional] index of LO branch to correct at q=0 
    :: [optional] path of prefix.dyn files (from matdyn.x) including filename w/out q-label
    :: REPLACE: if True, will prompt the user about editing the databases

    - Usage and main variables:
        :: yelph = YamboRefineElphDB(ylat, ... args ... )

        :: Has all relevant attributes for easy testing
        :: Contains testing functions to plot q-BZ map and phonon dispersion to be manually edited

        :: Uses branch disentangling to check for nearest-neighbour LO modes. This can be done
           by reading eigenvectors either from yambo or from prefix.dyn files.
           BEWARE!! Orderings may change in the two cases!!
    """ 
    def __init__(self,lattice,filename='ndb.elph_gkkp_expanded',folder_gkkp='SAVE',out_matdyn='prefix.freq',LO_ind=None,matdyn_filename='prefix.dyn',REPLACE=False):
        # Find correct database names
        filename = "%s/%s"%(folder_gkkp,filename)
        self.frag_filename = filename + "_fragment_"
        self.matdyn = out_matdyn        
        self.dyn_files = matdyn_filename
        if LO_ind is not None: self.LO = LO_ind

        # Check if databases exist. Exit if header is missing.
        try: database = Dataset(filename)
        except: raise FileNotFoundError("error opening header %s"%filename)
        try: 
            database_frag = Dataset("%s1"%self.frag_filename)
            database_frag.close()
        except: 
            raise FileNotFoundError("Database fragment at q=0 not detected")
    
        #read qpoints
        self.alat    = lattice.alat
        self.rlat    = lattice.rlat
        self.sym_car = lattice.sym_car
        self.qpoints = database.variables['PH_Q'][:].T
        self.car_qpoints = np.array([ q/self.alat for q in self.qpoints ])
        #read dimensions of electron phonon parameters
        self.nmodes, self.nqpoints, self.nkpoints, self.nbands = database.variables['PARS'][:4].astype(int)
        self.nmodes, self.nqpoints, self.nkpoints, b_1, b_2 = database.variables['PARS'][:5].astype(int)
        if b_1>b_2: # Old database (no GkkpBands in PARS)
            self.nbands = b_1
            self.b_in, self.b_out = [0,self.nbands-1]
        else: # New database (PARS with GkkpBands)
            self.b_in, self.b_out = [b_1-1,b_2-1]
            self.nbands = b_2-b_1+1
        self.natoms = int(self.nmodes/3)
        try: # Check if K-point list is provided (upon expansion)
            self.kpoints_elph = database.variables['PH_K'][:].T
            self.car_kpoints = np.array([ k/self.alat for k in self.kpoints_elph ])
            database.close()
        except KeyError:
            database.close()

        #Check how many databases are present
        self.nfrags = self.nqpoints
        for iq in range(self.nqpoints):
            if not os.path.isfile("%s%d"%(self.frag_filename,iq+1)):
                self.nfrags = iq
                break
        if self.nfrags!=self.nqpoints: raise FileNotFoundError("Found: %d fragments. Needed: %d"%(self.nfrags,self.nqpoints))

        #Check presence of matdyn output
        if not os.path.isfile(self.matdyn): raise FileNotFoundError("Matdyn 'prefix.freq' file not found")
        # Now the checks are done, we can continue
        # Frequencies part 
        self.read_frequencies_yambo()
        self.read_matdyn()
        self.expand_kpoints()
        self.expand_frequencies()
        # LO mode part
        if LO_ind is not None:
            self.find_nearest_qneighbor()
            self.find_nearest_modes(mode='yambo')
            self.get_LO_modes()

        # Now that we have everything, we can edit the DBs
        if REPLACE:
            self.replace_frequencies()
            self.replace_LO_modes()

    def read_frequencies_yambo(self):
        """
        Read phonon frequencies in eV
        """
        self.ph_energies  = np.zeros([self.nfrags,self.nmodes])

        for iq in range(self.nfrags):
            fil = self.frag_filename + "%d"%(iq+1)
            database = Dataset(fil,'r')
            self.ph_energies[iq] = np.sqrt(database.variables['PH_FREQS%d'%(iq+1)][:])*ha2ev
            database.close()

    def read_eigenvectors(self,iq,mode='matdyn'):
        """ Wrapper for either matdyn or yambo reading 
        """
        if mode=='matdyn': return self.read_eigenvectors_matdyn(iq)
        if mode=='yambo':  return self.read_eigenvectors_yambo(iq)
        

    def read_eigenvectors_yambo(self,iq):
        """
        Read phonon eigenvectors
        """
        ph_eigenvectors = np.zeros([self.nmodes,self.natoms,3],dtype=np.complex64)

        fil = self.frag_filename + "%d"%(iq+1)
        database= Dataset(fil)
        eigs_q = database.variables['POLARIZATION_VECTORS'][:].T
        ph_eigenvectors = eigs_q[0,:,:,:] + eigs_q[1,:,:,:]*1j
        database.close()
        
        return ph_eigenvectors

    def read_matdyn(self):
        """ Read prefix.freq file
            - Read qpoints (ibz)
            - Read frequencies (ibz)
            - Transform into bz
        """
        invcm2eV = 1.23981e-4 
        f = open(self.matdyn,'r')
        lines = f.readlines()
        Nq_ibz = int(lines[0].split()[-2])
        lines = lines[1:]
        Nlines = len(lines)
        if Nlines % Nq_ibz != 0: raise ValueError("Something wrong in matdyn file!")
        nstep = int(Nlines/Nq_ibz)
        qpts_matdyn = []
        freqs_matdyn = []
        for iq in range(Nq_ibz):
            qpts_matdyn.append([float(i) for i in lines[iq*nstep].split()])
            aux = []
            for i1 in range(1,nstep): aux += [ float(i) for i in lines[iq*nstep+i1].split() ]
            freqs_matdyn.append(aux)
        self.qpts_matdyn  = np.array(qpts_matdyn)/self.alat[0]
        self.freqs_matdyn = np.array(freqs_matdyn)*invcm2eV
        if np.any(self.freqs_matdyn<0.): print("[WARNING] NEGATIVE frequencies found in the matdyn file!")

    def expand_kpoints(self,atol=1e-6,verbose=1):
        """
        Expand matdyn qpoints
        """

        weights, kpoints_indexes, symmetry_indexes, kpoints_full = expand_kpoints(self.qpts_matdyn,self.sym_car,self.rlat,atol=atol)

        if verbose: print("%d kpoints expanded to %d"%(len(self.car_kpoints),len(kpoints_full)))

        #set the variables
        self.weights_ibz        = weights
        self.qpoints_indices    = kpoints_indexes
        self.symmetry_indices   = symmetry_indexes
        self.iku_matdyn_qpoints = np.array([k*self.alat for k in kpoints_full])
        self.car_matdyn_qpoints = kpoints_full

    def expand_frequencies(self):
        self.matdyn_ph_energies = np.zeros([self.nqpoints,self.nmodes])
        for iq in range(self.nqpoints): self.matdyn_ph_energies[iq] = self.freqs_matdyn[self.qpoints_indices[iq]]
    
    def read_eigenvectors_matdyn(self,iq):
        """
        Read eigenvectors from dyn files
        """
        ph_raw_energies = np.zeros(self.nmodes)
        ph_eigen        = np.zeros([self.nmodes,self.natoms,3],dtype=np.complex)
        f = open(self.dyn_files+str(iq+1), "r")
        start_line = -(self.natoms+1)*self.nmodes-1
        end_line = -1
        data = f.readlines()[start_line:end_line]
        # First read frequencies (convert to eV, so far unused)
        im=0
        for line in data:
            if 'freq' in line:
                ph_raw_energies[im] = float(line.split()[-2])/ev2cm1
                im+=1
        # Next read eigenmodes
        im=0
        for im in range(self.nmodes):
            for iat in range(self.natoms):
                ind = im*(self.natoms+1)+1+iat
                tmp = data[ind].split()
                ph_eigen[im,iat,0] = float(tmp[1])+1j*float(tmp[2])
                ph_eigen[im,iat,1] = float(tmp[3])+1j*float(tmp[4])
                ph_eigen[im,iat,2] = float(tmp[5])+1j*float(tmp[6])

        return ph_eigen 

    def find_nearest_qneighbor(self, max_NN=12, atol=1e-4):
        """ Finds the first nearest-neighbour points
            to q=(0,0,0) in the q-grid. This is intended to be
            the first shell of NN in a *regular grid*. 

            Function can be improved to find first NN across x,y,z if the
            grid is anisotropic.
            
            :: max_NN must be >= number of first NN
        """

        from sklearn.neighbors import KDTree
        
        qpts = self.car_qpoints

        # Build a KDTree from the grid points
        tree = KDTree(qpts)

        # Find the indices of the max_NN nearest neighbors of G
        distances, indices = tree.query(np.array([qpts[0]]), k=10)
        distances, indices = [distances[0][1:],indices[0][1:]]

        # Select only the first NN up to atol
        min_dist = min(distances)
        q_NN = []
        for iq, qdist in enumerate(distances):  
            if np.isclose(qdist,min_dist,atol=atol,rtol=atol): q_NN.append(indices[iq])

        self.q_NN = sorted(q_NN)
        self.N_NN = len(self.q_NN)

    def find_nearest_modes(self,mode='matdyn',verbose=1):
        """
        Find the indices of the LO-adjacent mode at q/=0
        corresponding to the first NN of q=0.

        We use the overlap matrix.
        """
        # q=0 modes
        eig0 = self.read_eigenvectors(0,mode=mode)
        # NN modes
        self.modes_NN = []
        for iq_bz in self.q_NN:
            iq_ibz = self.qpoints_indices[iq_bz]
            # (i) Read all the modes at iq_NN (in the IBZ)
            eigq = self.read_eigenvectors(iq_ibz,mode=mode)
            # (ii) Get the overlap of all the modes with those at q=0
            Overlap_matrix = phonon_overlap(eig0,eigq)
            # (iii) Find the mode with highest overlap with LO mode
            overlap = 0.
            best_mode = self.LO
            for im in range(self.nmodes):
                ref_value = np.abs(Overlap_matrix[self.LO,im])
                if ref_value > overlap:
                    overlap   = ref_value
                    best_mode = im
            # (iv) Obtain mode index of NN of LO modes
            self.modes_NN.append(best_mode)
            if verbose:
                if best_mode==self.LO:
                    print("==== No crossings detected ===")
                    print("Phonon mode #%d at q 0 remains phonon mode #%d at q_ibz #%d / q_bz #%d"%(self.LO,best_mode,iq_ibz,iq_bz))
                else:
                    print("==== Crossings detected ===")
                    print("Phonon mode #%d at q 0 becomes phonon mode #%d at q_ibz #%d / q_bz $%d"%(self.LO,best_mode,iq_ibz,iq_bz))

    def get_LO_modes(self):
        """ Obtain values for LO modes
        """
        # gkkp[q][k][bnd2][bnd1][mode][cmplx]
        dvscf_NN = np.zeros([self.N_NN,self.nkpoints,self.nbands,self.nbands,self.nmodes],dtype=np.complex64)

        for iq,Q in enumerate(self.q_NN):
            fil = self.frag_filename + "%d"%(Q+1)
            database = Dataset(fil)
            dvscf_aux = database.variables['ELPH_GKKP_Q%d'%(Q+1)][:]
            dvscf_NN[iq] = dvscf_aux[:,:,:,:,0] + 1j*dvscf_aux[:,:,:,:,1]
            database.close()

        # Check integrity of elph values
        if np.isnan(dvscf_NN).any(): print('[WARNING] NaN values detected in elph database.')
        
        # Compute modulus of NN g_LO(q_NN) and average it to g_LO_new
        g_LO_new = np.zeros([self.nkpoints,self.nbands,self.nbands])
        for iq in range(self.N_NN): g_LO_new += np.abs(dvscf_NN[iq,:,:,:,self.modes_NN[iq]])
        self.g_LO_new = g_LO_new/len(self.q_NN)    
        
        # Read dvscf at q=0    
        fil = self.frag_filename + "1"
        database = Dataset(fil)
        self.dvscf_0 = database.variables['ELPH_GKKP_Q1'][:]
        database.close()
        self.g_LO_old = self.dvscf_0[:,:,:,self.LO,0] + 1j*self.dvscf_0[:,:,:,self.LO,1]
        
    def replace_frequencies(self):
        """ 
        """
        usr_inp = input("Smooth frequencies from matdyn.x will replace raw ones in the ndb.elph_gkkp* databases. Proceed ['y','n']? ")
        if usr_inp != 'y':
            print("Frequencies not replaced.")
        else:
        
            aux_ph_energies = (self.matdyn_ph_energies/ha2ev)**2. 
        
            for iq in range(self.nfrags):
                fil = self.frag_filename + "%d"%(iq+1)
                database = Dataset(fil,'r+')
                database.variables['PH_FREQS%d'%(iq+1)][:] = aux_ph_energies[iq]
                database.close()
            
            print("Frequencies replaced.")

    def replace_LO_modes(self):
        """
        """
        usr_inp = input("The LO matrix elements at q=0 will be replaced with the |average| of their q nearest neighbors. Proceed ['y','n']? ")
        if usr_inp != 'y':
            print("LO matrix elements not replaced.")
        else:

            # Replace values for LO mode
            self.dvscf_0[:,:,:,self.LO,0] = self.g_LO_new
            self.dvscf_0[:,:,:,self.LO,1] = 0.
            
            # Write dvscf at q=0
            fil = self.frag_filename + "1"
            database = Dataset(fil,"r+")
            database.variables['ELPH_GKKP_Q1'][:] = self.dvscf_0
            database.close()        
                
            print("LO phonon modes replaced.")
        
    ### MATERIAL-DEPENDENT TEST FUNCTIONS ###
    def plot_qBZ(self,title="",save=False,mode='yambo'):
        """Plot annotated q-points in qBZ
        """
        if mode=='yambo':  points = self.car_qpoints
        if mode=='matdyn': points = self.car_matdyn_qpoints
        fig = plt.figure(figsize=(9,9))
        ax = plt.gca()
        ax.add_patch(BZ_Wigner_Seitz(self.rlat,color='black',linewidth=1.))
        ax.set_aspect('equal')
        ax.set_title(title)
        ax.scatter(points[:,0],points[:,1],marker='H',s=80,color='teal',linewidth=0.5,edgecolors='black',label='expanded')
        for i_k,kpt in enumerate(points):
            kx,ky = kpt[0],kpt[1]
            ax.annotate(i_k, (kx,ky), color='red', xytext=(kx+0.0005,ky+0.001), fontsize=5)
        if save: plt.savefig('%s.pdf'%title)
        plt.show()

    def find_path(self,indx_0,indx_1,band_indx,path='GMKG',mode='yambo'):
        """
        path: only GMKG

        :: indx_0 = list of high-symmetry point indices (e.g. G,M,K)
        :: indx_1 = list of high-symm. indx. in band path including end point
        :: band_indices = indices of band path 

        The above variables must be taken manually from plot_qBZ maps.

        Example for 12x12x1 grid:

        indx_0 = [0,31,142]
        indx_1 = [0,6,8,12]
        band_indx = [0,1,7,13,19,25,31,118,142,124,88,34,0] 
       
        """
        if path!='GMKG': raise ValueError("Only 'GMKG' path option is possible")

        if mode=='yambo': 
            pts=self.car_qpoints
            energies=self.ph_energies
        if mode=='matdyn':
            pts=self.car_matdyn_qpoints
            energies=self.matdyn_ph_energies

        iG, iM, iK   = indx_0
        Sindx        = indx_1
        band_indices = band_indx
        pts[iM] = np.abs(pts[iM])

        Np = len(band_indices)
        Nb = energies.shape[1]
        band_momenta = np.zeros(Np)
        band_energies= np.zeros([Np,Nb])
        nrm = np.linalg.norm
        for ip in range(1,Np): band_momenta[ip]  = band_momenta[ip-1]+nrm(pts[band_indices[ip]]-pts[band_indices[ip-1]])
        for ip in range(Np):   band_energies[ip] = energies[band_indices[ip],:] 

        return Sindx, band_momenta, band_energies

    def plot_bands(self,bnd_rng,indx_0,indx_1,band_indx,plt_show=True):
        """
        Phonon bands from gkkp (1) and matdyn (2)
        bnd_rng: [b_first, b_last]

        See above function for meaning of indx arguments.
        """
        Sindx, x1, y1 = self.find_path(indx_0,indx_1,band_indx)
        _,     x2, y2 = self.find_path(indx_0,indx_1,band_indx,mode='matdyn') 
        y1 = y1*1000.
        y2 = y2*1000.

        title='Phonon dispersion'
        clr1='teal'
        clr2='orange'
        flnm='Phon_Bands.pdf'
        ms1=10
        ms2=6
        lab1 = 'ndb.elph'
        lab2 = 'matdyn'

        fig = plt.figure(figsize=(8,6))
        ax = plt.gca()
        b_i, b_o = bnd_rng
        ax.set_xlim(x1[0]-0.01,x1[-1]+0.01)

        Spts = [x1[Sindx[0]],x1[Sindx[1]],x1[Sindx[2]],x1[Sindx[-1]]]
        Slabs = ['G','M','K','G']

        for PT in Spts: ax.axvline(PT,c='black')
        ax.set_xticks(Spts)
        ax.set_xticklabels(Slabs)
        ax.set_ylabel('Energy (meV)')
        for ib in range(b_i,b_o+1):
            if ib==b_i:
                ax.plot(x1,y1[:,ib],'.-',markersize=ms1,c=clr1,label=lab1)
                ax.plot(x2,y2[:,ib],'^-',markersize=ms2,c=clr2,label=lab2)
            else:
                ax.plot(x1,y1[:,ib],'.-',markersize=ms1,c=clr1)
                ax.plot(x2,y2[:,ib],'^-',markersize=ms2,c=clr2)

        plt.legend()
        plt.savefig(flnm)
        if plt_show: plt.show()

    ### END MATERIAL-DEPENDENT TEST FUNCTIONS ###
