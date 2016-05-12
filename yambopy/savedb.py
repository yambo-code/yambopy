# Copyright (c) 2016, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
from yambopy import *
from yambopy.netcdf import *
from yambopy.plot import *
from itertools import product
from scipy.optimize import bisect

max_exp = 50
ha2ev = 27.211396132
atol = 1e-3

def isbetween(a,b,c):
    #check if c is between a and b
    return np.isclose(np.linalg.norm(a-c)+np.linalg.norm(b-c)-np.linalg.norm(a-b),0,rtol=1e-05, atol=1e-06)

def err_handler(type, flag):
    print "Floating point error (%s), with flag %s" % (type, flag)

saved_handler = np.seterrcall(err_handler)

np.seterr(all='call')

def expand_kpts_val(kpts,syms,val):
    """ Take a list of qpoints and symmetry operations and return the full brillouin zone
    with the corresponding index in the irreducible brillouin zone
    """
    full_kpts = []
    full_val  = []
    print "nkpoints:", len(kpts)
    for nk,k in enumerate(kpts):
        for sym in syms:
            full_kpts.append((nk,np.dot(sym,k)))
            full_val.append(val[nk])

    return full_kpts, full_val

def vec_in_list(veca,vec_list):
    """ check if a vector exists in a list of vectors
    """
    return np.array([ np.allclose(veca,vecb,rtol=atol,atol=atol) for vecb in vec_list ]).any()

def expand_kpts(kpts,syms):
    """ Take a list of qpoints and symmetry operations and return the full brillouin zone
    with the corresponding index in the irreducible brillouin zone
    """
    full_kpts = []
    print "nkpoints:", len(kpts)
    for nk,k in enumerate(kpts):
        for sym in syms:
            full_kpts.append((nk,np.dot(sym,k)))

    return full_kpts

class YamboSaveDB():
    """ Read information from the SAVE database in Yambo
    """
    def __init__(self,save='SAVE'):
        #read database
        try:
            filename = '%s/ns.db1'%save
            database    = Dataset(filename)
        except:
            print "Error reading %s database"%filename
            exit()
        self.atomic_numbers   = database.variables['atomic_numbers'][:]
        self.atomic_positions = database.variables['ATOM_POS'][0,:]
        self.eigenvalues      = database.variables['EIGENVALUES'][0,:]*ha2ev
        self.sym_car          = database.variables['SYMMETRY'][:]
        self.kpts_iku         = database.variables['K-POINTS'][:].T
        self.lat              = database.variables['LATTICE_VECTORS'][:].T
        self.alat             = database.variables['LATTICE_PARAMETER'][:].T
        dimensions = database.variables['DIMENSIONS'][:]
        self.temperature = dimensions[13]
        self.electrons = dimensions[14]
        self.nkpoints  = int(dimensions[6])
        self.spin = int(dimensions[11])
        self.time_rev = dimensions[9]
        database.close()

        self.natoms = len(self.atomic_positions)
        
        #get a list of symmetries with time reversal
        nsym = len(self.sym_car)
        self.time_rev_list = [False]*nsym
        for i in xrange(nsym):
            self.time_rev_list[i] = ( i >= nsym/(self.time_rev+1) )

        #spin degeneracy if 2 components degen 1 else degen 2
        self.spin_degen = [0,2,1][int(self.spin)]

        #get minimum am maximul energies
        eiv = self.eigenvalues.flatten()
        self.min_eival = min(eiv)
        self.max_eival = max(eiv)

        #caclulate the reciprocal lattice
        self.rlat  = rec_lat(self.lat)
        self.nsym  = len(self.sym_car)

        #convert form internal yambo units to cartesian lattice units
        self.kpts_car = np.array([ k/self.alat for k in self.kpts_iku ])

        #convert cartesian transformations to reduced transformations
        inv = np.linalg.inv
        self.sym_rlu = np.zeros([self.nsym,3,3])
        for n,s in enumerate(self.sym_car):
            a = np.dot(s.T,inv(self.rlat))
            self.sym_rlu[n] = np.dot(inv(self.lat.T),a)

        #convert cartesian transformations to reciprocal transformations
        self.sym_rec = np.zeros([self.nsym,3,3])
        for n,s in enumerate(self.sym_car):
            self.sym_rec[n] = inv(s).T

        #status
        self.expanded = False
        self.efermi = None
        self.save = True

    def get_fermi(self,inv_smear=0.001):
        """ Determine the fermi energy
        """
        if self.efermi: return self.efermi

        def fermi(e):
            """ fermi dirac function
            """
            if e > max_exp:
                return 0
            elif e < -max_exp:
                return 1
            return 1/(np.exp(e)+1)

        def fermi_array(e_array,ef):
            """ Fermi dirac function for an array
            """
            e_array = (e_array-ef)/inv_smear
            return [ fermi(e) for e in e_array]

        def occupation_minus_ne(ef):
            """ The total occupation minus the total number of electrons
            """
            return sum([sum(self.spin_degen*fermi_array(self.eigenvalues[nk],ef))*self.weights[nk] for nk in xrange(self.nkpoints)])-self.electrons

        self.efermi = bisect(occupation_minus_ne,self.min_eival,self.max_eival)

        #import matplotlib.pyplot as plt
        #xmin, xmax = self.min_eival,self.max_eival
        #x = np.arange(xmin,xmax,0.1)
        #y = [occupation_minus_ne(i) for i in x ]
        #plt.plot(x,y)
        #plt.plot([xmin,xmax],[0,0])
        #plt.show()

        print "fermi: %lf eV"%self.efermi

        self.eigenvalues -= self.efermi
        self.min_eival -= self.efermi
        self.max_eival -= self.efermi

        self.occupations = np.zeros([self.nkpoints,self.nbands],dtype=np.float32)
        for nk in xrange(self.nkpoints):
            self.occupations[nk] = fermi_array(self.eigenvalues[nk,:self.nbands],0)

        return self.efermi

    def write_kpoints(self,filename_full='kpts_full.dat',filename='kpts.dat'):
        """ Write the kpoints in a file
        """
        kpts, nks, nss = self.expand_kpts()

        f = open(filename_full,'w')
        for k in kpts:
            f.write(("%12.8lf "*3)%tuple(k)+"\n")
        f.close()

        f = open(filename,'w')
        for k in self.kpts_car:
            f.write(("%12.8lf "*3)%tuple(k)+"\n")
        f.close()

    def get_path(self,path,kpts=None):
        """ Obtain a list of indexes and kpoints that belong to the regular mesh
        """
        if kpts is None:
            kpts, nks, nss = self.expand_kpts()
        else:
            nks = range(len(kpts))

        #points in cartesian coordinates
        path_car = red_car(path, self.rlat)

        #find the points along the high symmetry lines
        distance = 0
        bands_kpoints = []
        bands_indexes = []

        #for all the paths
        for k in range(len(path)-1):

            # store here all the points in the path
            # key:   has the coordinates of the kpoint rounded to 4 decimal places
            # value: index of the kpoint
            #        distance to the starting kpoint
            #        the kpoint cordinate
            kpoints_in_path = {}

            start_kpt = path_car[k]   #start point of the path
            end_kpt   = path_car[k+1] #end point of the path

            #generate repetitions of the brillouin zone
            for x,y,z in product(range(-1,2),range(-1,2),range(1)):

                #shift the brillouin zone
                shift = red_car([np.array([x,y,z])],self.rlat)[0]

                #iterate over all the kpoints
                for index, kpt in zip(nks,kpts):

                    kpt_shift = kpt+shift #shift the kpoint

                    #if the point is collinear we add it
                    if isbetween(start_kpt,end_kpt,kpt_shift):
                        key = tuple([round(kpt,4) for kpt in kpt_shift])
                        value = [ index, np.linalg.norm(start_kpt-kpt_shift), kpt_shift ]
                        kpoints_in_path[key] = value

            #sort the points acoording to distance to the start of the path
            kpoints_in_path = sorted(kpoints_in_path.values(),key=lambda i: i[1])

            #for all the kpoints in the path
            for index, disp, kpt in kpoints_in_path:
                bands_kpoints.append( kpt )
                bands_indexes.append( index )
                print ("%12.8lf "*3)%tuple(kpt), index

        self.bands_kpoints = bands_kpoints
        self.bands_indexes = bands_indexes
        self.bands_highsym_qpts = path_car

        return bands_kpoints, bands_indexes, path_car

    def expand_kpts(self):
        """ Take a list of qpoints and symmetry operations and return the full brillouin zone
        with the corresponding index in the irreducible brillouin zone
        """

        #check if the kpoints were already exapnded
        if self.expanded == True: return self.kpoints_full, self.kpoints_indexes, self.symmetry_indexes

        kpoints_indexes  = []
        kpoints_full     = []
        symmetry_indexes = []

        #kpoints in the full brillouin zone organized per index
        kpoints_full_i = {}

        #expand using symmetries
        for nk,k in enumerate(self.kpts_car):
            for ns,sym in enumerate(self.sym_car):
                new_k = np.dot(sym,k)

                #check if the point is inside the bounds
                k_red = car_red([new_k],self.rlat)[0]
                k_bz = (k_red+atol)%1

                #if the index in not in the dicitonary add a list
                if nk not in kpoints_full_i:
                    kpoints_full_i[nk] = []

                #if the vector is not in the list of this index add it
                if not vec_in_list(k_bz,kpoints_full_i[nk]):
                    kpoints_full_i[nk].append(k_bz)
                    kpoints_full.append(new_k)
                    kpoints_indexes.append(nk)
                    symmetry_indexes.append(ns)

        #calculate the weights of each of the kpoints in the irreducible brillouin zone
        self.full_nkpoints = len(kpoints_full)
        weights = np.zeros([self.nkpoints])
        for nk in kpoints_full_i:
            weights[nk] = float(len(kpoints_full_i[nk]))/self.full_nkpoints

        #set the variables
        self.expanded = True
        self.weights = np.array(weights)
        self.kpoints_full     = np.array(kpoints_full)
        self.kpoints_indexes  = np.array(kpoints_indexes)
        self.symmetry_indexes = np.array(symmetry_indexes)

        print "%d kpoints expanded to %d"%(len(self.kpts_car),len(kpoints_full))

        return self.kpoints_full, self.kpoints_indexes, self.symmetry_indexes

    def plot_bs(self,path):
        """ Plot the difference in energies of two bands
        """
        bands_kpoints, bands_indexes, bands_highsym_qpts = self.get_path(path)
        self.get_fermi()

        #calculate distances
        bands_distances = [0]
        distance = 0
        for nk in range(1,len(bands_kpoints)):
            distance += np.linalg.norm(bands_kpoints[nk]-bands_kpoints[nk-1])
            bands_distances.append(distance)

        #plot highsymetry qpoints
        distance = 0
        bands_highsym_qpts_distances = [0]
        for nk in range(1,len(bands_highsym_qpts)):
            plt.axvline(distance,color='k')
            distance += np.linalg.norm(bands_highsym_qpts[nk]-bands_highsym_qpts[nk-1])
            bands_highsym_qpts_distances.append(distance)

        plt.plot(bands_distances,self.eigenvalues[bands_indexes])
        plt.show()

    def plot_bs_bz(self,size=20,bandc=1,bandv=None,expand=True,repx=range(3),repy=range(3),repz=range(3)):
        """ Plot the difference in energies of two bands
        """
        if bandv is None: bandv = self.nbandsv

        cmap = plt.get_cmap("viridis")

        eigenvalues = self.eigenvalues
        print "tansitions %d -> %d"%(bandv,bandc)
        weights = (eigenvalues[:,bandc-1]-eigenvalues[:,bandv-1])
        print "min:", min(weights)
        print "max:", max(weights)
        weights = weights/max(weights)

        if expand:
            kpts, nks = self.expand_kpts(repx=repx,repy=repy,repz=repz)
            weights = weights[nks]
        else:
            kpts = self.kpts_car

        fig = plt.figure(figsize=(10,10))
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',serif="Computer Modern Roman",size=20)
        plt.scatter(kpts[:,0], kpts[:,1], s=size, marker='H', cmap=cmap, lw=0, c=weights)
        plt.axes().set_aspect('equal', 'datalim')
        plt.show()

    def __str__(self):
        s = ""
        s += "reciprocal lattice:\n"
        s += "\n".join([("%12.8lf "*3)%tuple(r) for r in self.rlat])+"\n"
        s += "lattice:\n"
        s += "\n".join([("%12.8lf "*3)%tuple(r) for r in self.lat])+"\n"
        s += "alat:\n"
        s += ("%12.8lf "*3)%tuple(self.alat)+"\n"
        s += "symmetry operations: %d\n"%len(self.sym_car)
        s += "temperature : %lf\n"%self.temperature
        s += "electrons   : %lf\n"%self.electrons
        return s
