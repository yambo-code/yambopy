# Copyright (c) 2018, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
import os
import numpy as np
from itertools import product
from netCDF4 import Dataset
from yambopy.plot import *
from yambopy.lattice import isbetween, car_red, rec_lat, vol_lat
from yambopy.units import ha2ev

max_exp = 50
atol = 1e-3

def expand_kpts_val(kpts,syms,val):
    """
    ake a list of qpoints and symmetry operations and return the full brillouin zone
    with the corresponding index in the irreducible brillouin zone
    """
    full_kpts = []
    full_val  = []
    print("nkpoints:", len(kpts))
    for nk,k in enumerate(kpts):
        for sym in syms:
            full_kpts.append((nk,np.dot(sym,k)))
            full_val.append(val[nk])

    return full_kpts, full_val

def vec_in_list(veca,vec_list):
    """
    Check if a vector exists in a list of vectors
    """
    return np.array([ np.allclose(veca,vecb,rtol=atol,atol=atol) for vecb in vec_list ]).any()

class YamboSaveDB():
    """
    Reads the information from the SAVE database in Yambo

    Arguments:

        ``save``: Path with the save folder (default:SAVE)
        ``filename``: name of the filename of the ns.db1 database created with yambo (default:ns.db1)

    **Properties:**

        ``atomic_numbers`` : atomic number of the species
        ``eigenvalues`` : eigenvalues of the electrons in eV
        ``nkpoints`` : number of kpoints
    """
    def __init__(self,atomic_numbers,atomic_positions,eigenvalues,sym_car,kpts_iku,
                 lat,alat,temperature,electrons,spin,time_rev):

        self.atomic_numbers   = atomic_numbers   
        self.atomic_positions = atomic_positions
        self.eigenvalues      = eigenvalues     
        self.sym_car          = sym_car         
        self.kpts_iku         = kpts_iku        
        self.lat              = lat             
        self.alat             = alat            
        self.temperature      = temperature     
        self.electrons        = electrons       
        self.spin             = spin            
        self.time_rev         = time_rev        

        #TODO: remove this
        self.expanded = False

    @classmethod
    def from_db_file(cls,folder='.',filename='ns.db1'):
        """
        Read the ns.db1 database
        """
        path_filename = os.path.join(folder,filename)
        if not os.path.isfile(path_filename):
            raise FileNotFoundError( "Error reading %s database in YamboSaveDB"%path_filename )

        with Dataset(path_filename) as database:
            
            dimensions            = database.variables['DIMENSIONS'][:]
            
            natoms = database.variables['N_ATOMS'][:].astype(int).T
            tmp_an = database.variables['atomic_numbers'][:].astype(int)

            args = dict( atomic_numbers   = [[tmp_an[n]]*na for n,na in enumerate(natoms)],
                         atomic_positions = database.variables['ATOM_POS'][0,:],
                         eigenvalues      = database.variables['EIGENVALUES'][0,:]*ha2ev,
                         sym_car          = database.variables['SYMMETRY'][:],
                         kpts_iku         = database.variables['K-POINTS'][:].T,
                         lat              = database.variables['LATTICE_VECTORS'][:].T,
                         alat             = database.variables['LATTICE_PARAMETER'][:].T,
                         temperature      = dimensions[13],
                         electrons        = dimensions[14],
                         spin             = int(dimensions[11]),
                         time_rev         = dimensions[9] )

        return cls(**args)

    @property
    def spin_degen(self):
        """spin degeneracy if 2 components degen 1 else degen 2"""
        return [0,2,1][int(self.spin)]

    @property
    def min_eival(self):
        return np.min(self.eigenvalues) 
    
    @property
    def max_eival(self):
        return np.max(self.eigenvalues)

    @property
    def car_kpoints(self):
        """convert form internal yambo units to cartesian lattice units"""
        return np.array([ k/self.alat for k in self.kpts_iku ])

    @property
    def red_kpoints(self):
        """convert from cartesian coordinates to reduced coordinates"""
        if not hasattr(self,"_red_kpoints"):
            self._red_kpoints = car_red(self.car_kpoints,self.rlat)
        return self._red_kpoints

    @property
    def rlat(self):
        """caclulate the reciprocal lattice"""
        return rec_lat(self.lat)

    @property
    def rlat_vol(self):
        return (2*np.pi)**3 * vol_lat(self.rlat)

    @property
    def lat_vol(self):
        return vol_lat(self.lat)

    @property
    def natoms(self):
        return len(self.atomic_positions)

    @property
    def nbands(self):
        _,nbands = self.eigenvalues.shape
        return nbands

    @property
    def nkpoints(self):
        return len(self.kpts_iku)

    @property
    def time_rev_list(self):
        """get a list of symmetries with time reversal"""
        time_rev_list = [False]*self.nsym
        for i in range(self.nsym):
            time_rev_list[i] = ( i >= self.nsym/(self.time_rev+1) )
        return time_rev_list

    @property
    def sym_rlu(self):
        """convert cartesian transformations to reduced transformations """
        sym_rlu = np.zeros([self.nsym,3,3])
        for n,s in enumerate(self.sym_car):
            a = np.dot(s.T,inv(self.rlat))
            sym_rlu[n] = np.dot(inv(self.lat.T),a)
        return sym_rlu

    @property
    def nsym(self):
        return len(self.sym_car)

    @property
    def sym_rec(self):
        """Convert cartesian transformations to reciprocal transformations"""
        sym_rec = np.zeros([self.nsym,3,3])
        for n,s in enumerate(self.sym_car):
            sym_rec[n] = np.linalg.inv(s).T
        return sym_rec

    @property
    def efermi(self):
        if not hasattr(self,"_efermi"):
            self._efermi = self.get_efermi
        return self._efermi

    def get_fermi(self,inv_smear=0.001,verbose=0):
        """ Determine the fermi energy
        """
        from scipy.optimize import bisect

        kpts, nks, nss = self.expand_kpts()

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
            return sum([sum(self.spin_degen*fermi_array(self.eigenvalues[nk],ef))*self.weights[nk] for nk in range(self.nkpoints)])-self.electrons

        efermi = bisect(occupation_minus_ne,self.min_eival,self.max_eival)

        if verbose: print("fermi: %lf eV"%efermi)

        self.eigenvalues -= efermi

        self.occupations = np.zeros([self.nkpoints,self.nbands],dtype=np.float32)
        for nk in range(self.nkpoints):
            self.occupations[nk] = fermi_array(self.eigenvalues[nk,:self.nbands],0)

        return efermi

    def write_kpoints(self,filename_full='kpts_full.dat',filename='kpts.dat'):
        """ Write the kpoints in a file
        """
        kpts, nks, nss = self.expand_kpts()

        f = open(filename_full,'w')
        for k in kpts:
            f.write(("%12.8lf "*3)%tuple(k)+"\n")
        f.close()

        f = open(filename,'w')
        for k in self.car_kpoints:
            f.write(("%12.8lf "*3)%tuple(k)+"\n")
        f.close()

    def get_path(self,path,kpts=None,debug=False):
        """ Obtain a list of indexes and kpoints that belong to the regular mesh
        """
        if kpts is None:
            kpts, nks, nss = self.expand_kpts()
        else:
            nks = list(range(len(kpts)))

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
            for x,y,z in product(list(range(-1,2)),list(range(-1,2)),list(range(1))):

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
            kpoints_in_path = sorted(list(kpoints_in_path.values()),key=lambda i: i[1])

            #for all the kpoints in the path
            for index, disp, kpt in kpoints_in_path:
                bands_kpoints.append( kpt )
                bands_indexes.append( index )
                if debug: print ("%12.8lf "*3)%tuple(kpt), index

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
        for nk,k in enumerate(self.car_kpoints):
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

        print("%d kpoints expanded to %d"%(len(self.car_kpoints),len(kpoints_full)))

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

    def plot_bs_bz(self,size=20,bandc=1,bandv=None,expand=True,repx=list(range(3)),repy=list(range(3)),repz=list(range(3))):
        """ Plot the difference in energies of two bands
        """
        if bandv is None: bandv = self.nbandsv

        cmap = plt.get_cmap("viridis")

        eigenvalues = self.eigenvalues
        print("tansitions %d -> %d"%(bandv,bandc))
        weights = (eigenvalues[:,bandc-1]-eigenvalues[:,bandv-1])
        print("min:", min(weights))
        print("max:", max(weights))
        weights = weights/max(weights)

        if expand:
            kpts, nks = self.expand_kpts(repx=repx,repy=repy,repz=repz)
            weights = weights[nks]
        else:
            kpts = self.car_kpoints

        fig = plt.figure(figsize=(10,10))
        plt.scatter(kpts[:,0], kpts[:,1], s=size, marker='H', cmap=cmap, lw=0, c=weights)
        ax = plt.axes()
        ax.set_aspect('equal')
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
