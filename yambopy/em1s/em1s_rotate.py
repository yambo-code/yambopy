# First version by Matteo Zanfrognini
# Revised and expanded by FP

from yambopy import *
from netCDF4 import Dataset

def expand_kpoints(kpoints,syms,rlat,atol=1e-6,verbose=0):
    """
    Expand reciprocal-space BZ vectors using lattice symmetries

    The expansion is consistent with the yambo expansion

    == Inputs ==
    :: kpoints: points in the IBZ to be expanded [CC]
    :: syms   : list of symmetry operations [CC]
    :: rlat   : reciprocal lattice vectors

    == Outputs ==
    :: List of expanded points [CC]
    :: Index table to go from unexpanded to expanded points
    :: Index table for symmetries
    :: List of weights

    """
    #check if the kpoints were already exapnded
    kpoints_indexes  = []
    kpoints_full     = []
    symmetry_indexes = []

    #kpoints in the full brillouin zone organized per index
    kpoints_full_i = {}

    sym_car_to_apply = syms

    #expand using symmetries
    for nk,k in enumerate(kpoints):

        #if the index in not in the dictionary add a list
        if nk not in kpoints_full_i:
            kpoints_full_i[nk] = []

        for ns,sym in enumerate(sym_car_to_apply):

            new_k = np.dot(sym,k)

            #check if the point is inside the bounds
            k_red = car_red([new_k],rlat)[0]
            k_bz = (k_red+atol)%1

            #if the vector is not in the list of this index add it
            if not vec_in_list(k_bz,kpoints_full_i[nk]):
                kpoints_full_i[nk].append(k_bz)
                kpoints_full.append(new_k)
                kpoints_indexes.append(nk)
                symmetry_indexes.append(ns)
                continue

    #calculate the weights of each of the kpoints in the irreducible brillouin zone
    nkpoints_full = len(kpoints_full)
    weights = np.zeros([nkpoints_full])
    for nk in kpoints_full_i:
        weights[nk] = float(len(kpoints_full_i[nk]))/nkpoints_full

    if verbose: print("%d kpoints expanded to %d"%(len(kpoints),len(kpoints_full)))

    #set the variables
    expanded_car_kpoints  = np.array(kpoints_full)
    kpoints_indices       = np.array(kpoints_indexes)
    symmetry_indices      = np.array(symmetry_indexes)
    weights_ibz           = np.array(weights)

    return expanded_car_kpoints,kpoints_indices,symmetry_indices,weights_ibz

def find_inversion_type(n_atoms,atom_pos,syms):
    """
    Check if the crystal has spatial inversion symmetry.

    Input:
    :: n_atoms is an array with the number of atoms per species
    :: atom_pos is an array of atomic positions per atom per species, already written in the COM frame

    Output:
    :: inv_type : 'spatial' or 'trev'
    :: inv_index: index of -1 symmetry
    """

    # Search for spatial inversion
    for i_species in range(len(n_atoms)): # N. species
        for i_atom in range(n_atoms[i_species]): # N. atoms/species
            atom1 = atom_pos[i_species][i_atom]
            for j_atom in range(n_atoms[i_species]):
                inv_type='spatial' # Assume we have it
                atom2 = atom_pos[i_species][j_atom]
                if np.allclose(atom1,-atom2): break            
                inv_type='trev' # If we get here, we don't have it

    # Search for inversion symmetry index
    for i_s,sym in enumerate(syms):
        if np.allclose(sym, -1.*np.identity(3)):
            inv_index = i_s
            break
    
    return inv_type, inv_index

class YamboEm1sRotate():
    """
    This class expands the em1s computed by yambo in the IBZ to the
    full BZ. 

    NOTE: the full BZ q-grid must be the one used by yambo, i.e., in the full Wigner-Seitz          cell of the crystal. 
          This grid can be generated in many ways (e.g. ypp -k in yambo) and has to be fed
          to the DFT codes for the no-symmetry calculation replacing the automatic grid

    :: Input 
    - YamboStaticScreeningDB object (IBZ calculation)
    - Location of ns.db1 database (default is inside SAVE)
    - [OPTIONAL] Location of output databases

    :: Output
    - numpy array with expanded static screening 
    - [OPTIONAL] ndb.em1s and ndb.em1s_fragment_* databases corresponding to the 
    full BZ.

    TODO:
    - table of S^-1 G [Not working in 3D as G-shell is not symmetric!?]
    - write new netCDF4 database  
    """

    def __init__(self,yem1s,save_path="SAVE",db1='ns.db1',path_output_DBs=None):

        supported_cutoffs = ['none','slab z']
        self.cutoff       = yem1s.cutoff

        if yem1s.cutoff not in supported_cutoffs: raise NotImplementedError("[ERROR] The em1s rotation is not currently implemented for cutoff %s."%yem1s.cutoff)

        # Attributes imported from StaticScreeningDB    
        alat              = yem1s.alat
        rlat              = yem1s.rlat
        self.qpoints_ibz  = yem1s.car_qpoints
        self.nqpoints_ibz = yem1s.nqpoints
        self.gvectors     = yem1s.gvectors
        self.red_gvectors = car_red(yem1s.gvectors,rlat)
        self.ngvectors    = yem1s.ngvectors
        self.X_ibz        = yem1s.X

        # Get symmetries in CC and real-space atomic positions
        if not os.path.isfile('%s/%s'%(save_path,db1)): raise FileNotFoundError("File %s not found."%db1)
        database = Dataset("%s/%s"%(save_path,db1), 'r')
        sym_car = np.array(database.variables['SYMMETRY'][:])
        self.syms = np.transpose(sym_car, (0, 2, 1))
        n_atoms =  database.variables['N_ATOMS'][:].astype(int)
        atom_pos = database.variables['ATOM_POS'][:]
        database.close()
        self.nsyms = len(self.syms)

        # Obtain transformed qpoints q'=Sq in the full BZ
        self.qpoints, self.qpoints_indices, self.syms_indices, _ = \
        expand_kpoints(self.qpoints_ibz,self.syms,rlat)
        self.nqpoints = len(self.qpoints)

        # Obtain transformed gvectors G'=S^{-1}G
        self.Sm1G_table = self.inverse_Gvector_table(rlat)

        # Spatial inversion or T-rev?
        # [WARNING] We assume one of the two is used!
        self.inv_type, self.inv_index = \
        find_inversion_type(n_atoms,atom_pos,self.syms)

        # Rotate em1s from IBZ to BZ
        print("=== Rotating em1s... ===")
        self.rotate_em1s()
        print("===      Done.       ===")

        if path_output_DBs is not None:
            print("=== Saving database... ==="
            self.outpath = path_output_DBs
            self.saveDBS()
            print("===        Done.       ==="
        else: print("(Enter value for path_output_DBs to print expanded ndb.em1s)")

    def inverse_Gvector_table(self,rlat):
        """
        Build table Sm1G_table such as:

        if ig_S = Sm1G_table[ig,iS], then S^{-1}G[ig]=G[ig_S]
        """ 
        inv_syms = np.linalg.inv(self.syms)
        Sm1G_table = np.zeros((self.ngvectors,len(inv_syms)),dtype=np.int)
        
        for iG,G in enumerate(self.gvectors):
            for i_S,sym in enumerate(inv_syms):
                G_S = np.dot(sym,G) 
                iG_S = self.get_g_index(G_S)
                if iG_S is None: raise ValueError("[ERROR] Problem in mapping inverse G-vectors. Try reducing isclose() tolerance in get_g_index")
                else: Sm1G_table[iG,i_S]=iG_S

        return Sm1G_table

    def get_g_index(self,g):
        """
        get the index of the gvectors.
        If the gvector is not present return None
        """

        for ng,gvec in enumerate(self.gvectors):
            if np.isclose(g,gvec).all():
                return ng
        return None

    def rotate_em1s(self):
        """
        Rotation of static screening

        The quantity rotated (content of ndb.em1s) is
 
        :math: D_{g1,g2}(q) = √v_g1(q) X_{g1g2}(q) √v_g2(q)   

        The rotations are performed according to

        - no time reversal (S)

            :math: D_{g1,g2}(Sq) = D_{S^-1g1,S^-1g2}(q)

        - time reversal is used (S->IS)

            :math:  D_{g1,g2}(ISq) = [ D_{(IS)^-1g1,(IS)^-1g2}(q) ]^*

        """
        X = np.zeros([self.nqpoints,self.ngvectors,self.ngvectors],dtype=np.complex64)
        
        for iq in range(self.nqpoints):
            iq_ibz = self.qpoints_indices[iq] # Index of untransformed q_ibz
            iS     = self.syms_indices[iq]    # Index of symmetry Sq=q_ibz
            for ig1 in range(self.ngvectors):
                Sm1_ig1 = self.Sm1G_table[ig1,iS] # index of G' such as G'=S-1G    
                for ig2 in range(self.ngvectors):
                    Sm1_ig2 = self.Sm1G_table[ig2,iS] 
                    # No TR
                    if self.inv_type=='spatial' or iS<self.inv_index:
                        X[iq,ig1,ig2]=self.X_ibz[iq_ibz,Sm1_ig1,Sm1_ig2]
                    # TR
                    if self.inv_type=='trev' and iS>=self.inv_index:
                        X[iq,ig1,ig2]=np.conj(self.X_ibz[iq_ibz,Sm1_ig1,Sm1_ig2])

        self.X = X                   

    def saveDBS(self,path):
        """
        Save the database
        """
        if os.path.isdir(path): shutil.rmtree(path)
        os.mkdir(path)

        #copy all the files
        oldpath = self.save
        filename = self.filename
        shutil.copyfile("%s/%s"%(oldpath,filename),"%s/%s"%(path,filename))
        for nq in range(self.nqpoints):
            fname = "%s_fragment_%d"%(filename,nq+1)
            shutil.copyfile("%s/%s"%(oldpath,fname),"%s/%s"%(path,fname))

        #edit with the new wfs
        X = self.X
        for nq in range(self.nqpoints):
            fname = "%s_fragment_%d"%(filename,nq+1)
            database = Dataset("%s/%s"%(path,fname),'r+')
            database.variables['X_Q_%d'%(nq+1)][0,0,:] = X[nq].real
            database.variables['X_Q_%d'%(nq+1)][0,1,:] = X[nq].imag
            database.close()

    def __str__(self):

        lines = []; app=lines.append
        app(marquee(self.__class__.__name__))

        app('nqpoints (ibz):   %d'%self.nqpoints_ibz)
        app('nqpoints (bz):   %d'%self.nqpoints)
        app('X size (G-space): %d'%self.ngvectors)
        app('cutoff:           %s'%self.cutoff)
        app('inversion type:   %s'%self.inv_type)
        app('nsymmetries:      %s'%self.nsyms)
        app('inversion sym is: %s'%self.inv_index)

        return "\n".join(lines)

