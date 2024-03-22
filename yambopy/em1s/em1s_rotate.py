#
# License-Identifier: GPL
#
# Copyright (C) 2024 The Yambo Team
#
# Authors: MZ, FP
# First version by MZ (2022), revised and expanded by FP (2023)
#
# This file is part of the yambopy project
#
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
            k_red[np.abs(k_red) < atol] = 0. # Set to zero values < atol to avoid mistakes
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
    full BZ (i.e., Wigner-Seitz cell). 

    === Usage and variables ===

    >> yem1s   = YamboStaticScreeningDB(save=db_path,em1s=db_path)
    >> yexpand = YamboEm1sRotate(yem1s,save_path=db_path,path_output_DBs=db_path.split('/')[0]+'/Expanded',verbose=1)

    :: Input 
    - yem1s           -> YamboStaticScreeningDB object (IBZ calculation)
    - save_path       -> [OPTIONAL] Location of ns.db1 database (default is inside SAVE)
    - path_output_DBs -> [OPTIONAL] Print expanded databases at location (default: do not print)
    - verbose         -> [OPTIONAL] If True, prints a list of kpoints for QE calculation of nosym system

    :: Output
    - numpy array with expanded static screening 
    - [OPTIONAL] ndb.em1s and ndb.em1s_fragment_* databases corresponding to the 
    full BZ.
    - [OPTIONAL] data file 'kpoints_bz.dat' containing expanded rlu kpoint coordinates in PW format

    """

    def __init__(self,yem1s,save_path="SAVE",db1='ns.db1',path_output_DBs=None,verbose=0):

        # Attributes imported from StaticScreeningDB    
        supported_cutoffs = ['none','slab z']
        self.cutoff       = yem1s.cutoff

        if yem1s.filename != 'ndb.em1s': raise NotImplementedError("[ERROR] The screening rotation is only implemented for ndb.em1s.")
        if yem1s.cutoff not in supported_cutoffs: raise NotImplementedError("[ERROR] The em1s rotation is not currently implemented for cutoff %s."%yem1s.cutoff)

        self.rlat         = yem1s.rlat
        self.alat         = yem1s.alat
        self.qpoints_ibz  = yem1s.car_qpoints
        self.nqpoints_ibz = yem1s.nqpoints
        self.gvectors     = yem1s.gvectors
        self.red_gvectors = yem1s.red_gvectors
        self.ngvectors    = yem1s.ngvectors
        self.X_ibz        = yem1s.X
        self.em1s_path    = yem1s.em1s

        self.k_output     = 'kpoints_bz.dat'

        # Get symmetries in CC and real-space atomic positions
        if not os.path.isfile('%s/%s'%(save_path,db1)): raise FileNotFoundError("File %s not found."%db1)
        database = Dataset("%s/%s"%(save_path,db1), 'r')
        sym_car = np.array(database.variables['SYMMETRY'][:])
        self.syms = np.transpose(sym_car, (0, 2, 1))
        n_atoms =  database.variables['N_ATOMS'][:].astype(int)
        atom_pos = database.variables['ATOM_POS'][:]
        if verbose: iku_kpoints_ibz = database.variables['K-POINTS'][:].T
        database.close()
        self.nsyms = len(self.syms)

        print("=== Rotating em1s... ===")
        print(" * Getting q-map...  ")

        # Obtain transformed qpoints q'=Sq in the full BZ
        self.qpoints, self.qpoints_indices, self.syms_indices, _ = \
        expand_kpoints(self.qpoints_ibz,self.syms,self.rlat)
        self.nqpoints = len(self.qpoints)

        print(" * Getting G-map ...  ")
        # Obtain transformed gvectors G'=S^{-1}G
        self.Sm1G_table = self.inverse_Gvector_table()

        print(" * Getting new em1s... ")

        # Spatial inversion or T-rev?
        # [WARNING] We assume one of the two is used!
        self.inv_type, self.inv_index = \
        find_inversion_type(n_atoms,atom_pos,self.syms)

        # Rotate em1s from IBZ to BZ
        self.rotate_em1s()

        if path_output_DBs is not None:
            print(" * Saving databases... [Remember that you might have to change the serial number!]")
            self.outpath = path_output_DBs
            self.saveDBS()
        else: print(" [!] Enter value for path_output_DBs to print expanded ndb.em1s")

        # Print kpts for PW nosym run
        if verbose: self.print_kpts_PW_format(iku_kpoints_ibz,units='rlu')

        print("===      Done.       ===")

    def inverse_Gvector_table(self,tol=1e-5):
        """
        Build table Sm1G_table such as:

        if ig_S = Sm1G_table[ig,iS], then S^{-1}G[ig]=G[ig_S]

        - kwargs are atol and rtol for np.isclose
        """
        inv_syms = np.linalg.inv(self.syms)
        Sm1G_table = np.zeros((self.ngvectors,len(inv_syms)),dtype=int)

        self.rotated_gvectors = np.zeros([len(self.syms),len(self.gvectors),3])
        for iG,G in enumerate(self.gvectors):
            for i_S,sym in enumerate(inv_syms):
                check = np.sum( (self.gvectors - np.dot(sym,G))**2., axis=1) < tol
                if np.sum(check) == 1.0: # One G-vector G' has been found to correspond to sym^{-1}G
                    Sm1G_table[iG,i_S]=int(np.where(check==1.0)[0])
                elif np.sum(check)!= 1.0: #None or multiple G-vectors have been found
                    raise ValueError("\n[ERROR] Problem in mapping inverse G-vectors. Try:\n - (i) reducing isclose() tolerance in get_g_index (easy case) \n - (ii) check that yambo packs G-shells correctly for your lattice type and G-cutoff (difficult case)")

        return Sm1G_table

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

    def saveDBS(self):
        """
        Write yambo-compatible ndb.em1s and ndb.em1s_fragment* 
        databases containing the expanded static screening.

        TODO: 
        - enable restarts
        - automatically change serial number to new SAVE

        WARNING:
            There are many quirks and special rules and version-dependent ifs in the yambo IO of ndb.em1s 
            and its fragments.

            The implementation here has to overcome the difficulty of indulging these quirks so that
            the yambopy-created database is not rejected by yambo.
            
            The present implementation tries to be as general as reasonably possible, and it has been tested,
            but instances may still arise with the appearance of errors and incompatibilities while trying to 
            get yambo to read these dbs. Also, in the event of an update in the yambo IO of ndb.em1s, these 
            writing functions will almost certainly be broken. 
        """
        path = self.outpath

        # If it exists, always remove outdir for safety
        if os.path.isdir(path): shutil.rmtree(path)
        os.mkdir(path)

        self.write_em1s_header(path)

        self.write_em1s_fragments(path)

    def write_em1s_header(self,path):
        """ Write ndb.em1s

            :: Dimensions ([F]ixed or [V]ariable):
               0 - D_0000000003 = 3   [F]-> HEAD_VERSION, HEAD_QPT, X_RL_vecs
               1 - D_0000000001 = 1   [F]-> HEAD_REVISION, SERIAL_NUMBER, HEAD_WF, CUTOFF, FRAGMENTED,
                                            Xs_energies_xc_KIND, Xs_wavefunctions_xc_KIND, GAUGE, 
                                            X_Time_ordering, X_TDDFT_KERNEL, X_OPTICAL_AVERAGE
               2 - D_0000000002 = 2   [F]-> SPIN_VARS, TEMPERATURES, X_DRUDE
               3 - D_0000000004 = 4   [F]-> HEAD_R_LATT
               4 - D_00000000Nk = Nk  [V]-> HEAD_QPT
               5 - D_0000000100 = 100 [F]-> CUTOFF, Xs_energies_xc_KIND, Xs_wavefunctions_xc_KIND,
                                            GAUGE, X_Time_ordering, X_TDDFT_KERNEL, X_OPTICAL_AVERAGE
               6 - D_0000000005 = 5   [F]-> X_PARS_1
               7 - D_0000000010 = 10  [F]-> X_PARS_3
               8 - D_00000000NG = NG  [V]-> X_RL_vecs

            :: Contents of X_PARS_1:
                - [V]X matrix size, X%ng 
                - [F]X band range, tmp_ib (start, end)
                - [F]X e/h energy range, X%ehe (start, end)  

            :: Contents of X_PARS_2: [Not usually written, we assume it's not present]

            :: Contents of X_PARS_3:
                - [F]X poles, X%cg_percentual,
                - [V]RL vectors in the sum, X%ngostnts
                - [F][r,Vnl] included, X%Vnl_included
                - [F]Longitudinal Gauge, local_long_gauge [normally absent]
                - [F]Field direction, X%q0 (x,y,z)
                - [F]BZ energy Double Grid, use_X_DbGd
                - [F]BZ energy DbGd points, X_DbGd_nkpts
                - [F]BZ Q point size factor, X_DbGd_percentual

        """
        def netcdftype(var_type):
            """ Distinguish between double and float
            """
            if var.dtype=='float32': return 'f4'
            elif var.dtype=='float64': return 'f8'
            else: raise TypeError('\n[ERROR] Variable type not recognized. It should be either float (float32) or double (float64).\n')
        
        # New database
        dbs = Dataset(path+'/ndb.em1s',mode='w',format='NETCDF4')

        # Old database
        dbs_ibz = Dataset(self.em1s_path+'/ndb.em1s')
        ibz_dims = dbs_ibz.dimensions.values()
        ibz_vars = dbs_ibz.variables.values() 

        # New dimensions

        # check if we are in an unfortunate cases due to how yambo writes the databases
        hard_dims = [3,1,2,4,100,5,10]  # fixed dimension sizes
        ibz_ = self.nqpoints_ibz in hard_dims #q_ibz number coincides with one of the fixed dimensions
        bz_  = self.nqpoints     in hard_dims #q_bz number coincides with one of the fixed dimensions
        gv_  = self.ngvectors    in hard_dims #g-vecs number coincides with one of the fixed dimensions

        sizes     = [dim.size for dim in ibz_dims]
        n_dims    = len(sizes)        
        ind_q_ibz_dim = sizes.index(self.nqpoints_ibz) 
        # Nq-related dimension is in fifth place, Gv-related one in last place       
        error_conditions= np.array([n_dims==9 and ind_q_ibz_dim!=4, n_dims<9 and ibz_==False and gv_==False]) 
        if error_conditions.any(): raise ValueError("[ERROR] something wrong with em1s dbs dimensions.")
        
        # Create dimensions
        iaux=0
        for dim in ibz_dims:
            # manage new q_bz dimension
            if iaux==4: 
                if ibz_==False and bz_==False:
                    # Replace old Nq_ibz dimension with new Nq_bz one
                    dbs.createDimension('D_%010d'%self.nqpoints,self.nqpoints)
                    iaux+=1
                if ibz_==True and bz_==False: 
                    # Add new Nq_bz dimension in between existing ones
                    dbs.createDimension('D_%010d'%self.nqpoints,self.nqpoints) # Create Nq_bz
                    dbs.createDimension(dim.name,dim.size)
                    iaux+=1
                if ibz_==False and bz_==True: 
                    # Do not copy the old unneeded Nq_ibz dimension
                    iaux+=1
                if ibz_==True and bz_==True:
                    # Copy existing dimensions as normal
                    dbs.createDimension(dim.name,dim.size)
                    iaux+=1
            # just copy existing dimensions including G-size
            else:
                dbs.createDimension(dim.name,dim.size)
                iaux+=1

        # Create variables  
        for var in ibz_vars:
            if var.name=='HEAD_QPT': dbs.createVariable('HEAD_QPT',netcdftype(var.dtype), ('D_%.10d'%3, 'D_%.10d'%self.nqpoints))
            elif var.name=='HEAD_R_LATT': dbs.createVariable('HEAD_R_LATT',netcdftype(var.dtype), ('D_%.10d'%4))
            else: dbs.createVariable(var.name, var.dtype, var.dimensions)
            
        # Store values in new DB, including new qpt coords in iku
        for var in ibz_vars: 
            if var.name=='HEAD_QPT':      dbs[var.name][:] = np.array([q*self.alat for q in self.qpoints]).T
            elif var.name=='HEAD_R_LATT': dbs[var.name][:] = [ self.nqpoints for i in range(4) ]
            else:                         dbs[var.name][:] = dbs_ibz[var.name][:]

        dbs.close()
        dbs_ibz.close()

    def write_em1s_fragments(self,path):
        """ Write ndb.em1s_fragment*

            :: Dimensions ([F]ixed or [V]ariable):
               0 - D_0000000006 = 6   [F]-> FREQ_PARS_sec_iqN
               1 - D_0000000001 = 2   [F]-> FREQ_sec_iqN, X_Q_N
               2 - D_0000000002 = 1   [F]-> FREQ_sec_iqN, X_Q_N
               8 - D_00000000NG = NG  [V]-> X_Q_N

            :: Contents of FREQ_PARS_sec_iqN:
                - [V]Current Q-pt index, iq
                - [F]X energy range, Xw%er (start, end)
                - [F]X damping range, Xw%dr (start, end)
                - [F]Number of frequencies, Xw%n_freqs
        """

        # Old database (reference for q=1)
        dbs_ibz = Dataset(self.em1s_path+'/ndb.em1s_fragment_1')
        ibz_dims = dbs_ibz.dimensions.values()
        ibz_vars = dbs_ibz.variables.values()
        pars = dbs_ibz.variables['FREQ_PARS_sec_iq1'][:]
        freq = dbs_ibz.variables['FREQ_sec_iq1'][:]

        # Create each fragment
        for iq_bz in range(self.nqpoints):

            iq_aux = iq_bz+1
            pars[0] = iq_aux 

            dbs_qbz = Dataset(path+'/ndb.em1s_fragment_%d'%iq_aux,mode='w',format='NETCDF4') 

            # Dimensions
            for dim in ibz_dims: dbs_qbz.createDimension(dim.name,dim.size)

            # Variables
            for var in ibz_vars: 
                if 'PARS' in var.name: dbs_qbz.createVariable('FREQ_PARS_sec_iq%d'%iq_aux, var.dtype, var.dimensions)
                if 'Q_se' in var.name: dbs_qbz.createVariable('FREQ_sec_iq%d'%iq_aux, var.dtype, var.dimensions)
                if 'X' in var.name:    dbs_qbz.createVariable('X_Q_%d'%iq_aux, var.dtype, var.dimensions)

            # Values
            dbs_qbz['FREQ_PARS_sec_iq%d'%iq_aux][:] = pars   
            dbs_qbz['FREQ_sec_iq%d'%iq_aux][:]      = freq

            # Finally, store expanded screening
            re, im = [self.X[iq_bz].real, self.X[iq_bz].imag]
            dbs_qbz['X_Q_%d'%iq_aux][:] = np.stack((re,im),axis=2)

            dbs_qbz.close()

        dbs_ibz.close()

    def print_kpts_PW_format(self,kpoints_iku,units='rlu'):
        """
        Print expanded kpoints in the Wigner-Seitz cell for PW
        nosym=.true., noinv=.true. calculation

        - kpoints_iku: iku kpoints in the IBZ read from ns.db1
        - units: either rlu or cc
        """
        def format_array(value):
            if value >= 0: return ' '+format(value,'.9f')
            else:          return format(value,'.9f')

        print(" * Printing PW-format kpoints file.")

        kpoints_ibz = np.array([ k/self.alat for k in kpoints_iku ])
        kpoints     = expand_kpoints(kpoints_ibz,self.syms,self.rlat)[0]
        if units=='rlu': points = car_red(kpoints,self.rlat)
        if units=='cc':  points = kpoints*self.alat[0]

        kpts2prnt = np.empty((len(points),4),dtype=object)
        for i in range(len(points)):
            for j in range(4):
                if j<3: 
                    val = points[i,j]
                    if val>=0: kpts2prnt[i,j]=' '+"{:0.9f}".format(val)
                    else:      kpts2prnt[i,j]="{:0.9f}".format(val)   
                if j==3:       kpts2prnt[i,j]='1'

        np.savetxt(self.k_output, kpts2prnt, fmt='%s %s %s %s')

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

