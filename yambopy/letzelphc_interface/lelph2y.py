import numpy as np
from netCDF4 import Dataset
import os
from yambopy.tools.string import marquee
from yambopy.lattice import rec_lat, red_car
from yambopy.units import ha2ev

def netcdftype(var_type,prec='float32'):
    """ netcdf4 types
    """
    if prec=='float32': ncprec=4
    if prec=='float64': ncprec=8

    if var_type=='f': return 'f%d'%ncprec
    if var_type=='i': return 'i%d'%ncprec
    if var_type=='s': return 'S1'

    else: raise TypeError('\n[ERROR] Variable type not recognized.\n')

class ConvertElectronPhononDB():
    """
    Convert a netCDF4 electron-phonon database to the yambo format.

    - Input: object related to the input database (OBJ)
    - Input: SAVE path. Needs to find ndb.kindx or ndb.gops and ns.db1
    - Input: OUT_path. Directory where to write dbs
    - Output: ndb.elph_expanded and ndb.elph_expanded_fragments databases

    OBJ must include:
    - nqpoints [coords]
    - nkpoints [coords]
    - nmodes
    - natoms
    - nspin
    - ph. energies [units]
    - ph. eigenvectors
    - el-ph mat. elements
    """
    def __init__(self,OBJ,code,SAVE_path,OUT_path=None):
        """
        1. Check SAVE and read needed variables
        2. Process OBJ
        3. Create header
        4. Create fragments
        """
        # By default place ndb.elph_gkkp_expanded* dbs in the SAVE directory
        if OUT_path is None: OUT_path=SAVE_path

        # Get yambo initialization and lattice info
        self.get_yambo_header_variables(SAVE_path)

        # Get el-ph data from external code
        match code:
            case 'lelphc': self.get_elph_variables_LELPHC(OBJ)
            case _: raise NotImplementedError("Code %s not found or implemented"%code)

        if not os.path.isdir(OUT_path): os.mkdir(OUT_path)

        # Write ndb.elph_gkkp_expanded in OUT_path
        self.write_header(OUT_path)

        # Write ndb.elph_gkkp_expanded_fragments_* in OUT_path
        self.write_fragments(OUT_path)

    def get_yambo_header_variables(self,SAVE_path):
            
        try: ns_db1 = Dataset(SAVE_path+'/ns.db1')
        except: raise FileNotFoundError("error opening %s"%(SAVE_path+'/ns.db1'))

        self.head_kpt     = ns_db1.variables["K-POINTS"][:]
        self.nkpoints_ibz = self.head_kpt.shape[1]
        self.nqpoints_ibz = self.nkpoints_ibz
        self.alat         = ns_db1.variables['LATTICE_PARAMETER'][:].T
        lat               = ns_db1.variables['LATTICE_VECTORS'][:].T
        self.rlat         = rec_lat(lat)
        self.noncollinear = ns_db1.variables["DIMENSIONS"][:][11]
        
        ns_db1.close()

        try: ndb_db = Dataset(SAVE_path+'/ndb.kindx')
        except FileNotFoundError: 
            try: ndb_db = Dataset(SAVE_path+'/ndb.gops')
            except: raise FileNotFoundError("ndb.gops or ndb.kindx not found in %s"%SAVE_path)

        self.head_version  = ndb_db.variables["HEAD_VERSION"][:]
        self.head_revision = ndb_db.variables["HEAD_REVISION"][:]
        self.serial_number = ndb_db.variables["SERIAL_NUMBER"][:]

        ndb_db.close()

    def get_elph_variables_LELPHC(self,OBJ):

        ph_k_red = OBJ.kpoints # To be converted in iku
        ph_q_red = OBJ.qpoints # To be converted in iku
        ph_k_car = red_car(ph_k_red,self.rlat) 
        ph_q_car = red_car(ph_q_red,self.rlat)
        self.ph_k = np.array([ k*self.alat for k in ph_k_car ]) 
        self.ph_q = np.array([ q*self.alat for q in ph_q_car ]) 
        self.ph_freqs    = ( OBJ.ph_energies/ha2ev )**2. # from eV to Hartree then squared
        self.max_ph_freq = np.max( OBJ.ph_energies/ha2ev )
        self.nkpoints_bz = OBJ.nk
        self.nqpoints_bz = OBJ.nq
        self.nmodes      = OBJ.nm
        self.natoms      = OBJ.nat
        self.nbands1     = OBJ.nb1
        self.nbands2     = OBJ.nb2
        self.bands       = OBJ.bands
        self.nspin       = OBJ.ns
        self.prec        = OBJ.ncfloat_type

        if self.nspin!=1: raise NotImplementedError("Spin-polarised case is not yet implemented")

        self.polarization_vectors = OBJ.ph_eigenvectors
        self.polarization_vectors = np.swapaxes(self.polarization_vectors,1,3) # swap to [iq][ix][iat][il] 
        self.elph_gkkp = OBJ.gkkp/2.**(1.5) # Convert from [Ry]^3/2 to [Ha]^3/2
        if self.nspin ==1: self.elph_gkkp = np.squeeze(self.elph_gkkp) # Remove spin dimension
        self.elph_gkkp = np.moveaxis(self.elph_gkkp,2,-1)    # swap to [iq][ik][ib1][ib2][im]

        # For now the unneeded E_K_PLUS_Q variable is filled with zeros
        self.e_k_plus_q = np.zeros([self.nqpoints_bz,self.nspin,self.nkpoints_bz,self.nbands2])

    def write_header(self,OUT_path):
           
        # PARS: modes, qpts, kpts, [bnds / bnds_0 bnds_1], using_q_grid=T, hosting_bare_gkkp=F, hosting_DW=F
        if self.bands[0] == 1: self.pars = [self.nmodes,self.nqpoints_bz,self.nkpoints_bz,self.nbands1,True,False,False]
        else:   self.pars = [self.nmodes,self.nqpoints_bz,self.nkpoints_bz,self.bands[0],self.bands[1],True,False,False]
        len_pars = len(self.pars)
        self.spin_vars = [self.nspin,self.noncollinear] #YAMBO: SPIN_vec_disk=(/n_sp_pol,n_spinor/)
        self.head_r_latt = [self.nkpoints_ibz,self.nkpoints_bz,self.nqpoints_ibz,self.nqpoints_bz]
        self.fragmented = True

        # New database
        dbs = Dataset(OUT_path+'/ndb.elph_gkkp_expanded',mode='w',format='NETCDF4')

        # Create dimensions
        dbs.createDimension('D_%010d'%3,3)
        dbs.createDimension('D_%010d'%1,1)
        dbs.createDimension('D_%010d'%2,2)
        dbs.createDimension('D_%010d'%4,4)
        for value in [self.natoms,self.nkpoints_ibz,len_pars,self.nqpoints_bz,self.nkpoints_bz]:
            if value not in [1,2,3,4]: 
                try: dbs.createDimension('D_%010d'%value,value)
                except RuntimeError: pass # This is when one of the dimensions already exists

        # Create variables
        p=self.prec
        dbs.createVariable('HEAD_VERSION',netcdftype('f',p), ('D_%.10d'%3))
        dbs.createVariable('HEAD_REVISION',netcdftype('f',p), ('D_%.10d'%1))
        dbs.createVariable('SERIAL_NUMBER',netcdftype('f',p), ('D_%.10d'%1))
        dbs.createVariable('SPIN_VARS',netcdftype('f',p), ('D_%.10d'%2))
        dbs.createVariable('HEAD_R_LATT',netcdftype('f',p), ('D_%.10d'%4))
        dbs.createVariable('HEAD_KPT',netcdftype('f',p), ('D_%.10d'%3, 'D_%.10d'%self.nkpoints_ibz))
        dbs.createVariable('FRAGMENTED',netcdftype('f',p), ('D_%.10d'%1))
        PARS = dbs.createVariable('PARS',netcdftype('f',p), ('D_%.10d'%len_pars))
        dbs.createVariable('MAX_PH_FREQ',netcdftype('f',p), ('D_%.10d'%1))
        dbs.createVariable('PH_Q',netcdftype('f',p), ('D_%.10d'%3, 'D_%.10d'%self.nqpoints_bz))
        dbs.createVariable('PH_K',netcdftype('f',p), ('D_%.10d'%3, 'D_%.10d'%self.nkpoints_bz))

        # Fill variables
        dbs['HEAD_VERSION'][:]  = self.head_version
        dbs['HEAD_REVISION'][:] = self.head_revision
        dbs['SERIAL_NUMBER'][:] = self.serial_number
        dbs['SPIN_VARS'][:]     = self.spin_vars
        dbs['HEAD_R_LATT'][:]   = self.head_r_latt
        dbs['HEAD_KPT'][:]      = self.head_kpt
        dbs['FRAGMENTED'][:]    = self.fragmented
        dbs['PARS'][:]          = self.pars
        dbs['MAX_PH_FREQ'][:]   = self.max_ph_freq
        dbs['PH_Q'][:]          = self.ph_q.T
        dbs['PH_K'][:]          = self.ph_k.T

        dbs.close()

    def write_fragments(self,OUT_path):

        # New fragments
        for iq in range(self.nqpoints_bz):

            iq_prnt = iq+1

            dbs = Dataset(OUT_path+'/ndb.elph_gkkp_expanded_fragment_%d'%iq_prnt,mode='w',format='NETCDF4')

            # Create dimensions
            dbs.createDimension('D_%010d'%3,3)
            dbs.createDimension('D_%010d'%1,1)
            dbs.createDimension('D_%010d'%2,2)
            for value in [self.nmodes,self.natoms,self.nbands2,self.nkpoints_bz]:
                if value not in [1,2,3]: 
                    try: dbs.createDimension('D_%010d'%value,value)
                    except RuntimeError: pass # This is when one of the dimensions already exists

            # Create variables
            p=self.prec
            dbs.createVariable('PH_FREQS%d'%iq_prnt,netcdftype('f',p), ('D_%.10d'%self.nmodes))
            dbs.createVariable('POLARIZATION_VECTORS',netcdftype('f',p), ('D_%.10d'%3,'D_%.10d'%self.natoms,'D_%.10d'%self.nmodes,'D_%.10d'%2))
            dbs.createVariable('E_K_PLUS_Q%d'%iq_prnt,netcdftype('f',p), ('D_%.10d'%1,'D_%.10d'%self.nkpoints_bz,'D_%.10d'%self.nbands2))
            dbs.createVariable('ELPH_GKKP_Q%d'%iq_prnt,netcdftype('f',p), ('D_%.10d'%self.nkpoints_bz,'D_%.10d'%self.nbands2,'D_%.10d'%self.nbands2,'D_%.10d'%self.nmodes,'D_%.10d'%2))

            # Fill variables
            pol_tmp = np.reshape(self.polarization_vectors[iq],(3,self.natoms,self.nmodes,1))
            polarization_vectors_Q = np.concatenate((pol_tmp.real,pol_tmp.imag),axis=3)
            elph_tmp = np.reshape(self.elph_gkkp[iq],(self.nkpoints_bz,self.nbands2,self.nbands2,self.nmodes,1))
            elph_gkkp_Q = np.concatenate((elph_tmp.real,elph_tmp.imag),axis=4)

            dbs['PH_FREQS%d'%iq_prnt][:]    = self.ph_freqs[iq]
            dbs['POLARIZATION_VECTORS'][:]  = polarization_vectors_Q
            dbs['E_K_PLUS_Q%d'%iq_prnt][:]  = self.e_k_plus_q[iq]
            dbs['ELPH_GKKP_Q%d'%iq_prnt][:] = elph_gkkp_Q

            dbs.close()
