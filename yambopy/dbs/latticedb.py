from __future__ import print_function
from __future__ import division
# Copyright (c) 2015, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
from builtins import range
from builtins import object
from past.utils import old_div
from yambopy import *
from netCDF4 import Dataset
import itertools
import operator

atol = 1e-6

class YamboLatticeDB(object):
    """
    Class to read the lattice information from the netcdf file
    """
    def __init__(self, save='SAVE',filename='ns.db1',expand=True):
        self.filename = '%s/%s'%(save,filename)
        self.readDB()
        # generate additional structure using the data read from the DBs
        self.process()
        if expand: self.expandKpoints()
            
    def readDB(self):
        try:
            database = Dataset(self.filename)
        except:
            print("error opening %s in YamboLatticeDB"%self.filename)
            exit()

        self.lat         = database.variables['LATTICE_VECTORS'][:].T
        self.alat        = database.variables['LATTICE_PARAMETER'][:].T
        natoms           = database.variables['N_ATOMS'][:].astype(int).T
        self.atomic_numbers   = database.variables['atomic_numbers'][:].astype(int)
        self.atomic_positions = database.variables['ATOM_POS'][0,:]
        self.sym_car     = database.variables['SYMMETRY'][:]
        self.iku_kpoints = database.variables['K-POINTS'][:].T
        dimensions = database.variables['DIMENSIONS'][:]
        self.temperature = dimensions[13]
        self.nelectrons = dimensions[14]
        self.nkpoints  = int(dimensions[6])
        self.spin = int(dimensions[11])
        self.time_rev = dimensions[9]
        #atomic numbers
        atomic_numbers = [[self.atomic_numbers[n]]*na for n,na in enumerate(natoms)]
        self.atomic_numbers = list(itertools.chain.from_iterable(atomic_numbers))
        #atomic masses
        self.atomic_masses = [atomic_mass[a] for a in self.atomic_numbers]

        database.close()
        
    def process(self):
        inv = np.linalg.inv
        #caclulate the reciprocal lattice
        self.rlat  = rec_lat(self.lat)
        self.nsym  = len(self.sym_car)
        
        #convert form internal yambo units to cartesian lattice units
        self.car_kpoints = np.array([ old_div(k,self.alat) for k in self.iku_kpoints ])
        self.red_kpoints = car_red(self.car_kpoints,self.rlat)
        self.nkpoints = len(self.car_kpoints)
        
        #convert cartesian transformations to reciprocal transformations
        self.sym_rec = np.zeros([self.nsym,3,3])
        for n,s in enumerate(self.sym_car):
            self.sym_rec[n] = inv(s).T
            
        #get a list of symmetries with time reversal
        nsym = len(self.sym_car)
        self.time_rev_list = [False]*nsym
        for i in range(nsym):
            self.time_rev_list[i] = ( i >= old_div(nsym,(self.time_rev+1)) )
        
    def expandKpoints(self):
        """
        Take a list of qpoints and symmetry operations and return the full brillouin zone
        with the corresponding index in the irreducible brillouin zone
        """

        #check if the kpoints were already exapnded
        kpoints_indexes  = []
        kpoints_full     = []
        symmetry_indexes = []

        #kpoints in the full brillouin zone organized per index
        kpoints_full_i = {}

        #expand using symmetries
        for nk,k in enumerate(self.car_kpoints):
            #if the index in not in the dicitonary add a list
            if nk not in kpoints_full_i:
                kpoints_full_i[nk] = []
                    
            for ns,sym in enumerate(self.sym_car):
                
                new_k = np.dot(sym,k)

                #check if the point is inside the bounds
                k_red = car_red([new_k],self.rlat)[0]
                k_bz = (k_red+atol)%1
                
                #if the vector is not in the list of this index add it
                if not vec_in_list(k_bz,kpoints_full_i[nk]):
                    kpoints_full_i[nk].append(k_bz)
                    kpoints_full.append(new_k)
                    kpoints_indexes.append(nk)
                    symmetry_indexes.append(ns)
                    continue

        #calculate the weights of each of the kpoints in the irreducible brillouin zone
        self.full_nkpoints = len(kpoints_full)
        weights = np.zeros([self.nkpoints])
        for nk in kpoints_full_i:
            weights[nk] = old_div(float(len(kpoints_full_i[nk])),self.full_nkpoints)

        print("%d kpoints expanded to %d"%(len(self.car_kpoints),len(kpoints_full)))

        #set the variables
        self.weights_ibz      = np.array(weights)
        self.car_kpoints      = np.array(kpoints_full)
        self.red_kpoints      = car_red(self.car_kpoints,self.rlat)
        self.kpoints_indexes  = np.array(kpoints_indexes)
        self.symmetry_indexes = np.array(symmetry_indexes)

    def get_path(self,path,kpts=None,debug=False):
        """
        Obtain a list of indexes and kpoints that belong to the regular mesh
        """
        nks  = range(self.nkpoints)
        kpts = self.car_kpoints
        print nks
        print kpts

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
                if debug: print ("%12.8lf "*3)%tuple(kpt), index

        self.bands_kpoints = bands_kpoints
        self.bands_indexes = bands_indexes
        self.bands_highsym_qpts = path_car

        return bands_kpoints, bands_indexes, path_car


