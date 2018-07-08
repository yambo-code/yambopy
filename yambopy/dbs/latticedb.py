# Copyright (c) 2018, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
from itertools import product, chain
import numpy as np
from netCDF4 import Dataset
from yambopy.tools.jsonencoder import JsonDumper, JsonLoader
from yambopy.lattice import rec_lat, car_red, red_car, vec_in_list, isbetween
from yambopy.units import atomic_mass
from yambopy.tools.string import marquee
from qepy.lattice import Path

class YamboLatticeDB():
    """
    Class to read the lattice information from the netcdf file
    """
    def __init__(self,lat=None,alat=None,sym_car=None,iku_kpoints=None,
                      atomic_positions=None,atomic_numbers=None,time_rev=None):
        self.lat = np.array(lat)
        self.alat = np.array(alat)
        self.sym_car = np.array(sym_car)
        self.iku_kpoints = np.array(iku_kpoints)
        self.atomic_positions = np.array(atomic_positions)
        self.atomic_numbers = np.array(atomic_numbers)
        self.time_rev = time_rev

    @classmethod
    def from_db_file(cls,filename):
        """ Initialize YamboLattice from a local dbfile """
        y = cls()
        y.read_db(filename)
        y._process()
        y.expand_kpoints()
        return y
 
    @classmethod    
    def from_dict(cls,data):
        """ Initialize instance of the class from a dictionary
        """
        lat = data["lat"]
        alat = data["alat"]
        sym_car = data["sym_car"]
        iku_kpoints = data["iku_kpoints"]
        atomic_positions = data["atomic_positions"]
        atomic_numbers = data["atomic_numbers"]
        time_rev = data["time_rev"]
        y = cls(lat,alat,sym_car,iku_kpoints,atomic_positions,atomic_numbers,time_rev)
        y._process()
        y.expand_kpoints()
        return y

    @classmethod
    def from_json_file(cls,filename):
        data = JsonLoader(filename)
        return cls.from_dict(data)
 
    @property
    def nkpoints(self):
        return len(self.car_kpoints)

    def read_db(self,filename):
        try:
            database = Dataset(filename)
        except:
            raise IOError("error opening %s in YamboLatticeDB"%filename)

        #lattice data
        self.lat         = database.variables['LATTICE_VECTORS'][:].T
        self.alat        = database.variables['LATTICE_PARAMETER'][:].T
        self.sym_car     = database.variables['SYMMETRY'][:]
        self.iku_kpoints = database.variables['K-POINTS'][:].T

        #atomic numbers
        natoms         = database.variables['N_ATOMS'][:].astype(int).T
        self.atomic_positions = database.variables['ATOM_POS'][:,0,:]
        atomic_numbers = database.variables['atomic_numbers'][:].astype(int)
        atomic_numbers = [[atomic_numbers[n]]*na for n,na in enumerate(natoms)]
        self.atomic_numbers = list(chain.from_iterable(atomic_numbers))
        self.atomic_masses = [atomic_mass[a] for a in self.atomic_numbers]
        
        dimensions = database.variables['DIMENSIONS'][:]
        self.temperature = dimensions[13]
        self.nelectrons = dimensions[14]
        self.spin = int(dimensions[11])
        self.time_rev = dimensions[9]

        database.close()

    def as_dict(self):
        """ get the information of this class as a dictionary
        """
        data = {"lat" : self.lat,
                "alat" : self.alat,
                "sym_car" : self.sym_car,
                "iku_kpoints" : self.iku_kpoints,
                "atomic_positions" : self.atomic_positions,
                "atomic_numbers" : self.atomic_numbers,
                "time_rev": self.time_rev }
        return data

    def write_json(self,filename):
        """ write a json file with the lattice information """
        JsonDumper(self.as_dict(),filename)

    def _process(self):
        """
        Generate additional information from data in the database
        """
        inv = np.linalg.inv
        #caclulate the reciprocal lattice
        self.rlat  = rec_lat(self.lat)
        self.nsym  = len(self.sym_car)

        #convert form internal yambo units to cartesian lattice units
        self.car_kpoints = np.array([ k/self.alat for k in self.iku_kpoints ])
        self.red_kpoints = car_red(self.car_kpoints,self.rlat)

        #convert cartesian transformations to reciprocal transformations
        self.sym_rec = np.zeros([self.nsym,3,3])
        for n,s in enumerate(self.sym_car):
            self.sym_rec[n] = inv(s).T

        #get a list of symmetries with time reversal
        nsym = len(self.sym_car)
        self.time_rev_list = [False]*nsym
        for i in range(nsym):
            self.time_rev_list[i] = ( i >= nsym/(self.time_rev+1) )

    def expand_kpoints(self,atol=1e-6):
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
            weights[nk] = float(len(kpoints_full_i[nk]))/self.full_nkpoints

        print("%d kpoints expanded to %d"%(len(self.car_kpoints),len(kpoints_full)))

        #set the variables
        self.weights_ibz      = np.array(weights)
        self.car_kpoints      = np.array(kpoints_full)
        self.red_kpoints      = car_red(self.car_kpoints,self.rlat)
        self.kpoints_indexes  = np.array(kpoints_indexes)
        self.symmetry_indexes = np.array(symmetry_indexes)

    def get_path(self,path,debug=False):
        """
        Obtain a list of indexes and kpoints that belong to the regular mesh
        """
        if isinstance(path,Path):
            path = path.get_klist()

        nks  = list(range(self.nkpoints))
        kpts = self.car_kpoints

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
                if debug: print(("%12.8lf "*3)%tuple(kpt), index)

        self.bands_kpoints = bands_kpoints
        self.bands_indexes = bands_indexes
        self.bands_highsym_qpts = path_car

        return bands_kpoints, bands_indexes, path_car

    def __str__(self):
        lines = []; app = lines.append
        app(marquee(self.__class__.__name__))
        app("\nlattice:")
        lines += [("%12.8lf " * 3) % tuple(vec) for vec in self.lat]
        app("\natom positions:")
        for an, pos in zip(self.atomic_numbers, self.atomic_positions):
            app( "%3d " % an + ("%12.8lf " * 3) % tuple(pos) )
        return "\n".join(lines)
