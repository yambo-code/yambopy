# By FP (2023)

from yambopy import *
from yambopy.lattice import point_matching
from netCDF4 import Dataset

class YamboQPRotate():
    """
    This class expands the ndb.QP computed by yambo in the IBZ to the
    full k-BZ (i.e., Wigner-Seitz cell).

    === Usage and variables ===

    >> yqp  = YamboQPDB.from_db(filename='ndb.QP',folder='db_path')
    >> yexpand = YamboQPRotate(yqp,save_path=db_path,path_output_DBs=db_path.split('/')[0]+'/Expanded',verbose=1)

    :: Input
    - yqp             -> YamboQPDB object (IBZ calculation)
    - save_path       -> [OPTIONAL] Location of ns.db1 database (default is inside SAVE)
    - path_output_DBs -> [OPTIONAL] Print expanded databases at location (default: do not print)
    - verbose         -> [OPTIONAL] If True, prints a list of kpoints for QE calculation of nosym system

    :: Output
    - numpy array with expanded quasiparticle energies
    - [OPTIONAL] ndb.QP database corresponding to the full BZ.
    - [OPTIONAL] data file 'kpoints_bz.dat' containing expanded rlu kpoint coordinates in PW format

    """
    def __init__(self,yqp,save_path="SAVE",db1='ns.db1',path_output_DBs=None,verbose=0):

        # Quantities to store from QPDB
        self.nkibz, self.nbands = yqp.eigenvalues_qp.shape
        self.qp_path            = yqp.qp_path

        self.k_output     = 'kpoints_bz.dat'

        # Read geometry and BZ expanded k-grid
        lattice     = YamboLatticeDB.from_db_file(save_path+"/ns.db1")
        self.map_bz2ibz  = lattice.kpoints_indexes
        self.nkbz        = lattice.nkpoints
        self.alat        = lattice.alat
        self.rlat        = lattice.rlat
        self.iku_kpoints = np.array(lattice.iku_kpoints)

        self.nstates     = self.nkbz*self.nbands

        print("=== Rotating ndb.QP... ===")

        print(" * QP kpts and table... ")
        # Check that IBZ k-grid is complete and matches ns.db1
        if len(yqp.qps['Kpoint'])!=self.nkibz: 
            raise ValueError('\n[ERROR] Expansion currently works if ndb.QP is computed for all IBZ kpoints.')
        k_matches = point_matching(yqp.qps['Kpoint'],lattice.ibz_kpoints)
        if not k_matches.tolist() == [ ik for ik in range(self.nkibz) ]:
             raise TypeError('\n[ERROR] Something wrong with QP kpts. It seems k-indices are not matching with SAVE data.')

        self.qp_table = self.expand_QPtable(yqp.qps['qp_table'])
        print(" * KS energies... ")
        self.Eo = self.expand_QP(yqp.qps['Eo'].real)
        print(" * QP energies and lifetimes... ")
        self.E = self.expand_QP(yqp.qps['E'])
        print(" * Z renorm factors... ")
        self.Z = self.expand_QP(yqp.qps['Z'])

        if path_output_DBs is not None:
            print(" * Saving databases... ")
            self.outpath = path_output_DBs
            self.saveDBS()
        else: print(" [!] Enter value for path_output_DBs to print expanded ndb.QP")

        # Print kpts for PW nosym run
        if verbose: self.print_kpts_PW_format(lattice.iku_kpoints,units='rlu')

        print("===      Done.       ===")

    def expand_QP(self,values_ibz):
        """
        Expand a QP quantity
        """

        values_bz = np.zeros((self.nkbz*self.nbands),dtype=values_ibz.dtype)
        for ib in range(self.nbands):
            for ikbz in range(self.nkbz): 
                bz_index  = ib+ikbz*self.nbands
                ibz_index = ib+self.map_bz2ibz[ikbz]*self.nbands
                values_bz[bz_index]= values_ibz[ibz_index]

        return values_bz

    def expand_QPtable(self,table):
        """
        QP table works like this:
        table[0][i_QP]=i_band1
        table[1][i_QP]=i_band2 (= i_band1 since SE is diagonal)
        table[2][i_QP]=i_kpt
        """
        table_b1  = self.expand_QP(table[0])
        table_b2  = self.expand_QP(table[1])
        table_kpt = np.zeros(self.nstates)

        for ib in range(self.nbands):
            for ikbz in range(self.nkbz):
                bz_index  = ib+ikbz*self.nbands
                table_kpt[bz_index]=ikbz

        table_bz = np.array([table_b1,table_b2,table_kpt])

        return table_bz

    def saveDBS(self):
        """
        Write yambo-compatible ndb.QP databases containing the expanded quasiparticles.

        WARNING:
            There are many quirks and special rules and version-dependent ifs in the yambo IO 
            of ndb.QP.

            The implementation here has to overcome the difficulty of indulging these quirks so that
            the yambopy-created database is not rejected by yambo.

            The present implementation tries to be as general as reasonably possible, and it has been tested,
            but instances may still arise with the appearance of errors and incompatibilities while trying to
            get yambo to read these dbs. Also, in the event of an update in the yambo IO of ndb.QP, these
            writing functions will almost certainly be broken.
        
        ==========================================================================================
        Write ndb.QP

            :: Dimensions ([F]ixed or [V]ariable):
               0 - D_0000000001 = 1   [F]-> QP_DB_kind, HEAD_REVISION, SERIAL_NUMBER, CUTOFF,
                                            G_energies_xc_KIND,G_wavefunctions_xc_KIND,Xp_energies_xc_KIND,
                                            Xp_wavefunctions_xc_KIND,QP_N_DESCRIPTORS,QP_GW_solver,
                                            QP_GW_approximation,QP_PPA_imaginary_Energy,QP_GW_SC_iterations,
                                            QP_dS_dw_steps,QP_dS_dw_step,QP_X_Gs,QP_X_poles,QP_X_xc-Kernel,
                                            QP_X_BZ_energy_Double_Grid,QP_Sc_G_damping,QP_Sc_bands_terminator,
                                            QP_Sx_RL_components,QP_EMPTY_STR_Nr18,QP_GF_energies_kind,
                                            QP_GF_WFs_kind,QP_Xs_energies_kind,QP_Xs_WFs_kind
               1 - D_0000000001 = 3   [F]-> HEAD_VERSION, HEAD_D_LATT, QP_table, QP_kpts
               2 - D_0000000002 = 2   [F]-> SPIN_VARS, TEMPERATURES, QP_X_bands, QP_X_e_h_E_range,
                                            QP_Sc_G_bands,QP_QP_@_state_1_K_range,QP_QP_@_state_1_b_range,
                                            QP_E,QP_Z
               3 - D_0000000100 = 100 [F]-> G_energies_xc_KIND,G_wavefunctions_xc_KIND,Xp_energies_xc_KIND,
                                            Xp_wavefunctions_xc_KIND,QP_GW_solver,QP_GW_approximation,
                                            QP_X_xc-Kernel,QP_EMPTY_STR_Nr18,QP_GF_energies_kind,QP_GF_WFs_kind,
                                            QP_Xs_energies_kind,QP_Xs_WFs_kind
               4 - D_0000000006 = 6   [F]-> PARS
               5 - QP_desc_size = 24  [F]-> QP_DESCRIPTORS_SIZES,QP_DESCRIPTORS_NAMES,QP_DESCRIPTORS_KINDS,
                                            QP_DESCRIPTORS_TERMS
               6 - QP_string_len = 100 [F]-> QP_DESCRIPTORS_NAMES,QP_DESCRIPTORS_KINDS,QP_DESCRIPTORS_TERMS
               7 - D_0000000Nqp = Nqp  [V]-> QP_table,QP_E,QP_Eo,QP_Z
               8 - D_00000000Nk = Nk   [V]-> QP_kpts,QP_E,QP_Eo,QP_Z

            :: Contents of PARS:
                - [F]Nbands, qp%nb
                - [V]Nkpts in IBZ, qp%nk
                - [V]Nqp=Nkpts*Nbands, qp%n_states
                - [F]GWo SC iterations, GWo_iterations
                - [F]GW SC iterations, GW_iterations
                - [F]QP_desc_size, qp%desc%n

        ===========================================================================================
        """
        def netcdftype(var_type):
            """ Distinguish between double and float
            """
            if var.dtype=='float32': return 'f4'
            elif var.dtype=='float64': return 'f8'
            else: raise TypeError('\n[ERROR] Variable type not recognized. It should be either float (float32) or double (float64).\n')

        path = self.outpath
        # If it exists, always remove outdir for safety
        if os.path.isdir(path): shutil.rmtree(path)
        os.mkdir(path)

        # New database
        dbs = Dataset(path+'/ndb.QP',mode='w',format='NETCDF4')

        # Old database
        dbs_ibz = Dataset(self.qp_path+'/ndb.QP')
        ibz_dims = dbs_ibz.dimensions.values()
        ibz_vars = dbs_ibz.variables.values()

        # check if we are in an unfortunate cases due to how yambo writes the databases
        sizes     = [dim.size for dim in ibz_dims]
        hard_dims = sizes[:-2] # fixed dimension sizes
        ibz_ = self.nkibz   in hard_dims #k_ibz number coincides with one of the fixed dimensions
        iqp_ = self.nkibz*self.nbands in hard_dims #qp ibz number coincides with one of the fixed dimensions
        bz_  = self.nkbz    in hard_dims #k_bz number coincides with one of the fixed dimensions
        qp_  = self.nstates in hard_dims #qp number coincides with one of the fixed dimensions

        sizes     = [dim.size for dim in ibz_dims]
        n_dims    = len(sizes)
        ind_k_ibz_dim = sizes.index(self.nkibz)
        # Nk-related dimension is in last place, Nqp-related one in 7th place
        error_conditions=np.array([n_dims==9 and ind_k_ibz_dim!=8, n_dims<9 and ibz_==False,n_dims<9 and iqp_==False])
        if error_conditions.any(): raise ValueError("[ERROR] something wrong with QP dbs dimensions.")        

        # Create dimensions
        iaux=0
        for dim in ibz_dims:
            # new k_bz dimension
            if iaux==6 and qp_==False: dbs.createDimension('D_%010d'%self.nstates,self.nstates)
            elif iaux==6 and qp_==True: continue
            if iaux==7 and bz_==False: dbs.createDimension('D_%010d'%self.nkbz,self.nkbz)
            elif iaux==7 and bz_==True: continue
            # fixed dimensions
            else: dbs.createDimension(dim.name,dim.size)
            iaux+=1

        # New variables
        pars = dbs_ibz.variables['PARS'][:]
        pars[1] = self.nkbz
        pars[2] = self.nstates 

        QP_E = np.stack((self.E.real, self.E.imag),axis=1)
        QP_Z = np.stack((self.Z.real, self.Z.imag),axis=1)

        # Create variables
        for var in ibz_vars:
            if var.name=='QP_table': dbs.createVariable('QP_table',netcdftype(var.dtype), ('D_%.10d'%3, 'D_%.10d'%self.nstates))
            elif var.name=='QP_kpts': dbs.createVariable('QP_kpts',netcdftype(var.dtype), ('D_%.10d'%3, 'D_%.10d'%self.nkbz))
            elif var.name=='QP_E': dbs.createVariable('QP_E',netcdftype(var.dtype), ('D_%.10d'%self.nstates, 'D_%.10d'%2))
            elif var.name=='QP_Eo': dbs.createVariable('QP_Eo',netcdftype(var.dtype), ('D_%.10d'%self.nstates))
            elif var.name=='QP_Z': dbs.createVariable('QP_Z',netcdftype(var.dtype), ('D_%.10d'%self.nstates, 'D_%.10d'%2))
            else: dbs.createVariable(var.name, var.dtype, var.dimensions)

        # Store values in new DB, including new qpt coords in iku
        for var in ibz_vars:
            if   var.name=='PARS':     dbs[var.name][:] = pars
            elif var.name=='QP_table': dbs[var.name][:] = self.qp_table
            elif var.name=='QP_kpts':  dbs[var.name][:] = self.iku_kpoints.T
            elif var.name=='QP_E':     dbs[var.name][:] = QP_E
            elif var.name=='QP_Eo':    dbs[var.name][:] = self.Eo
            elif var.name=='QP_Z':     dbs[var.name][:] = QP_Z
            elif var.name=='QP_QP_@_state_1_K_range': dbs[var.name][:] = [ 1, self.nkbz ]
            else:                      dbs[var.name][:] = dbs_ibz[var.name][:]

        dbs.close()
        dbs_ibz.close()

    def print_kpts_PW_format(self,kpoints_iku,units='rlu'):
        """
        Print expanded kpoints in the Wigner-Seitz cell for PW
        nosym=.true., noinv=.true. calculation

        - kpoints_iku: iku kpoints in the BZ read from ns.db1
        - units: either rlu or cc
        """
        def format_array(value):
            if value >= 0: return ' '+format(value,'.9f')
            else:          return format(value,'.9f')

        print(" * Printing PW-format kpoints file.")

        kpoints = np.array([ k/self.alat for k in kpoints_iku ])
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

        app('nkpoints (ibz):  %d'%self.nkibz)
        app('nkpoints (bz):   %d'%self.nkbz)
        app('nbands:          %d'%self.nbands)
        app('nqpstates (ibz): %d'%(self.nkibz*self.nbands))
        app('nqpstates (bz):  %d'%self.nstates)

        return "\n".join(lines)
