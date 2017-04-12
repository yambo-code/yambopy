# Copyright (c) 2017, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
from yambopy import *
from math import sqrt
from time import time

def abs2(x):
    return x.real**2 + x.imag**2
    
class YamboDipolesDB():
    """
    Open the dipoles databases and store it in a DipolesDB class
    """
    def __init__(self,lattice,save='SAVE',filename='ndb.dip_iR_and_P',expand=True,dip_type='iR'):
        self.lattice = lattice
        self.filename = "%s/%s"%(save,filename)
        
        #read dipoles
        try:
            database = Dataset(self.filename, 'r')
        except:
            raise IOError("Error opening %s in YamboDipolesDB"%self.filename)
            
        self.nq_ibz, self.nq_ibz, self.nk_ibz, self.nk_bz = database.variables['HEAD_R_LATT'][:].astype(int)
        self.spin = database.variables['SPIN_VARS'][1].astype(int)

        # indexv is the maximum partially occupied band
        # indexc is the minimum partially empty band
        self.min_band, self.max_band, self.indexv, self.indexc = database.variables['PARS'][:4].astype(int)
        database.close()

        # determine the number of bands
        self.nbands  = self.max_band-self.min_band+1
        self.nbandsv = self.indexv-self.min_band+1
        self.nbandsc = self.max_band-self.indexc+1
        
        #read the database
        self.dipoles = self.readDB(dip_type)

        #expand the dipoles to the full brillouin zone
        if expand:
            self.expandDipoles(self.dipoles)

    def normalize(self,electrons):
        """ 
        Use the electrons to normalize the dipole matrix elements
        """
        eiv = electrons.eigenvalues
        nkpoints, nbands = eiv.shape
        for nk in xrange(nkpoints):
            
            eivk = eiv[nk]
            
            #create eigenvalues differences arrays
            norm = np.array([ [ec-ev for ev in eivk] for ec in eivk  ])
            
            #normalize
            for i,j in product(xrange(nbands),repeat=2):
                if norm[i,j] == 0: 
                    self.dipoles[nk,:,i,j] = 0
                else:
                    self.dipoles[nk,:,i,j] = self.dipoles[nk,:,i,j]/norm[i,j]
        dipoles = self.dipoles

    def readDB(self,dip_type):
        """
        The dipole matrix has the following indexes:
        [nkpoints, 3, nspin, nbands conduction, nbands valence]
        """
        self.dip_type = dip_type
        dipoles = np.zeros([self.nk_ibz,3,self.nbandsc,self.nbandsv],dtype=np.complex64)
        
        #check dipole db format
        filename = "%s_fragment_1"%(self.filename)
        db = Dataset(filename)
        tag1 = 'DIP_iR_k_0001_spin_0001'
        tag2 = 'DIP_iR_k_0001_xyz_0001_spin_0001'
        if tag1 in db.variables.keys():
            dipoles_format = 1
        elif tag2 in db.variables.keys():
            dipoles_format = 2
        db.close()
        
        for nk in range(self.nk_ibz):

            #open database for each k-point
            filename = "%s_fragment_%d"%(self.filename,nk+1)
            db = Dataset(filename)

            if dipoles_format == 1:
                dip = db.variables['DIP_%s_k_%04d_spin_%04d'%(dip_type,nk+1,1)][:].view(dtype=np.complex64)[:,:,:,0]
                for i in xrange(3):
                    dipoles[nk,i] = dip[:,:,i].T
            elif dipoles_format == 2:
                for i in xrange(3):
                    dip = db.variables['DIP_%s_k_%04d_xyz_%04d_spin_%04d'%(dip_type,nk+1,i+1,1)][:]
                    dipoles[nk,i] = dip[0].T+dip[1].T*1j

            #close database
            db.close()

        return dipoles
        
    def expandDipoles(self,dipoles=None,field_dir=[1,0,0],field_dir3=[0,0,1]):
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
        indexv = self.min_band-1
        indexc = self.indexc-1
        nbands = self.min_band+self.nbands-1
        
        #Note that P is Hermitian and iR anti-hermitian.
        if self.dip_type == 'P':
            factor =  1.0
        else:
            factor = -1.0
            
        self.dipoles = np.zeros([nkpoints,3,nbands,nbands],dtype=np.complex64)
        for nk_fbz,nk_ibz,ns in zip(xrange(nkpoints),nks,nss):
            
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
            
            for c,v in product(xrange(self.nbandsc),xrange(self.nbandsv)):
                #rotate dipoles
                self.dipoles[nk_fbz,:,indexc+c,indexv+v] = np.dot(tra,dip[:,c,v])
        
            #make hermitian
            for c,v in product(xrange(self.nbandsc),xrange(self.nbandsv)):
                self.dipoles[nk_fbz,:,indexv+v,indexc+c] = factor*np.conjugate(self.dipoles[nk_fbz,:,indexc+c,indexv+v])
                        
        self.field_dirx = field_dirx
        self.field_diry = field_diry
        self.field_dirz = field_dirz
        
        return dipoles, kpts
       
    def plot(self,ax,kpoint=0,dir=0,func=abs2):
        return ax.matshow(func(self.dipoles[kpoint,dir]))
        
    def plot_dp_bs(self,path,bandv,bandc,dir=0,spin=0):
        """ Plot the dipole along a path
        """
        dipoles, kpts = self.get_dipoles(expand=True)

        bands_kpoints, bands_indexes, bands_highsym_qpts = self.get_path(path,kpts)

        #calculate distances
        bands_distances = [0]
        distance = 0
        for nk in range(1,len(bands_kpoints)):
            distance += np.linalg.norm(bands_kpoints[nk-1]-bands_kpoints[nk])
            bands_distances.append(distance)

        #check if bandc is an integer
        if type(bandc) is int: bandc = (bandc,)
        #check if bandv in an integer
        if type(bandv) is int: bandv = (bandv,)

        print "tansitions %s -> %s"%(str(bandv),str(bandc))

        #get kpoints in path
        dipole_bands = dipoles[bands_indexes]

        #get spin
        dipole_bands = dipole_bands[:,:,spin]

        #project along q?
        dipole_bands = dipole_bands[:,dir,:,:]

        #absolute value
        dipole_bands = np.absolute(dipole_bands[:,:,:])**2

        #plot highsymetry qpoints
        distance = 0
        for nk in range(1,len(bands_highsym_qpts)):
            plt.axvline(distance,color='k')
            distance+=np.linalg.norm(bands_highsym_qpts[nk]-bands_highsym_qpts[nk-1])

        for c,v in zip(bandc,bandv):
            plt.plot(bands_distances,dipole_bands[:,c-1,v-1],label="%d -> %d"%(c,v))
        plt.legend()
        plt.show()

    def plot_dp_bz(self,bandv,bandc,size=20,dir=0,spin=0,field_dir=[1,0,0],expand=True,repx=range(1),repy=range(1),repz=range(1)):
        """ Plot the weights in a scatter plot of this exciton
        """
        cmap = plt.get_cmap("viridis")

        dipoles, kpts = self.get_dipoles(field_dir=field_dir,expand=expand)

        #check if bandc is an integer
        if type(bandc) is int: bandc = (bandc,)
        #check if bandv in an integer
        if type(bandv) is int: bandv = (bandv,)

        print "tansitions %s -> %s"%(str(bandv),str(bandc))
        weights = np.zeros([len(kpts)])
        for c,v in zip(bandc,bandv):
            weights += np.absolute(dipoles[:,dir,spin,c-1,v-1])**2
        weights = np.array([x for x in weights])
        weights = weights

        fig = plt.figure(figsize=(10,10))
        print len(kpts)
        for x,y,z in product(range(-1,2),range(-1,2),range(1)):
            #shift the brillouin zone
            d = red_car([np.array([x,y,z])],self.lattice.rlat)[0]

            plt.scatter(kpts[:,0]+d[0], kpts[:,1]+d[1], s=size, marker='H', lw=0, cmap=cmap, c=weights)
        plt.axes().set_aspect('equal', 'datalim')
        plt.colorbar()

        plt.show()

    def __str__(self):
        s = ""
        s += "\nkpoints:\n"
        s += "nk_ibz : %d\n"%self.nk_ibz
        s += "nk_bz  : %d\n"%self.nk_bz
        s += "\nnumber of bands:\n"
        s += "nbands : %d\n" % self.nbands
        s += "nbandsv: %d\n" % self.nbandsv
        s += "nbandsc: %d\n" % self.nbandsc
        s += "indexv : %d\n" % (self.min_band-1)
        s += "indexc : %d\n" % (self.indexc-1)
        s += "field_dirx: %10.6lf %10.6lf %10.6lf\n"%tuple(self.field_dirx)
        s += "field_diry: %10.6lf %10.6lf %10.6lf\n"%tuple(self.field_diry)
        s += "field_dirz: %10.6lf %10.6lf %10.6lf\n"%tuple(self.field_dirz)
        return s

if __name__ == "__main__":
    ddb = DipolesDB()
    ddb.get_databases()
    print ddb
