#
# License-Identifier: GPL
#
# Copyright (C) 2024 The Yambo Team
#
# Authors: HPC, FP
#
# This file is part of the yambopy project
#
import os
import numpy as np
from itertools import product
from netCDF4 import Dataset
from yambopy.tools.jsonencoder import JsonDumper, JsonLoader
from yambopy.lattice import vol_lat, rec_lat, car_red, red_car, vec_in_list, isbetween
from yambopy.units import atomic_mass
from yambopy.tools.string import marquee
from qepy.lattice import Path

class YamboLatticeDB(object):
    """
    Class to read the lattice information from the netcdf file
    """
    def __init__(self,lat=None,alat=None,sym_car=None,iku_kpoints=None,
                      car_atomic_positions=None,atomic_numbers=None,time_rev=None):
        self.lat                  = np.array(lat)
        self.alat                 = np.array(alat)
        self.sym_car              = np.array(sym_car)
        self.iku_kpoints          = np.array(iku_kpoints)
        self.car_atomic_positions = np.array(car_atomic_positions)
        self.atomic_numbers       = np.array(atomic_numbers)
        self.time_rev             = time_rev
        self.ibz_nkpoints         = len(iku_kpoints)

    @classmethod
    def from_db(cls,filename='ns.db1',Expand=True,expand_mode=1,atol=1e-6):
        return cls.from_db_file(filename,Expand,expand_mode,atol)
    
    @classmethod
    def from_db_file(cls,filename='ns.db1',Expand=True,expand_mode=1,atol=1e-6):
        """ Initialize YamboLattice from a local dbfile """

        if not os.path.isfile(filename):
            raise FileNotFoundError("error opening %s in YamboLatticeDB"%filename)

        with Dataset(filename) as database:

            dimensions = database.variables['DIMENSIONS'][:]

            natoms_a = database.variables['N_ATOMS'][:].astype(int).T
            tmp_an = database.variables['atomic_numbers'][:].astype(int)
            tmp_apos = database.variables['ATOM_POS'][:,:]
            # Prevent case where n. of atomic species > n. of species in ATOM_POS
            i_at_unused_types = [i_at for i_at, at_type_n in enumerate(natoms_a) if at_type_n == 0]
            if len(i_at_unused_types)>0:
                print('[WARNING] Found unused atomic type(s): ',i_at_unused_types)
                natoms_a = [at_type_n for at_type_n in natoms_a if at_type_n!=0]
                tmp_an = [tmp_an[i] for i in range(len(tmp_an)) if tmp_an[i] not in i_at_unused_types ]

            flatten = lambda l: [item for sublist in l for item in sublist]
            atomic_numbers = flatten([[tmp_an[n]]*na for n,na in enumerate(natoms_a)])
            atomic_positions = np.vstack([[tmp_apos[n,ia] for ia in range(na)] for n,na in enumerate(natoms_a) ])

            args = dict( atomic_numbers       = atomic_numbers,
                         car_atomic_positions = atomic_positions,
                         sym_car              = database.variables['SYMMETRY'][:],
                         iku_kpoints          = database.variables['K-POINTS'][:].T,
                         lat                  = database.variables['LATTICE_VECTORS'][:].T,
                         alat                 = database.variables['LATTICE_PARAMETER'][:].T,
                         time_rev             = dimensions[9] )

        y = cls(**args)
        if Expand: y.expand_kpoints(expand_mode,atol=atol)
        return y
 
    @classmethod    
    def from_dict(cls,data):
        """ Initialize instance of the class from a dictionary
        """
        lat = data["lat"]
        alat = data["alat"]
        sym_car = data["sym_car"]
        iku_kpoints = data["iku_kpoints"]
        atomic_positions = data["car_atomic_positions"]
        atomic_numbers = data["atomic_numbers"]
        time_rev = data["time_rev"]

        y = cls(lat,alat,sym_car,iku_kpoints,atomic_positions,atomic_numbers,time_rev)
        return y

    @classmethod
    def from_json_file(cls,filename):
        data = JsonLoader(filename)
        return cls.from_dict(data)
 
    @property
    def nkpoints(self):
        return len(self.car_kpoints)

    @property
    def red_atomic_positions(self):
        return car_red(self.car_atomic_positions,self.lat)

    def as_dict(self):
        """ get the information of this class as a dictionary
        """
        data = {"lat" : self.lat,
                "alat" : self.alat,
                "sym_car" : self.sym_car,
                "iku_kpoints" : self.iku_kpoints,
                "car_atomic_positions" : self.car_atomic_positions,
                "atomic_numbers" : self.atomic_numbers,
                "time_rev": self.time_rev }
        return data

    def write_json(self,filename):
        """ write a json file with the lattice information """
        JsonDumper(self.as_dict(),filename)

    @property
    def iku_kpoints(self):
        return self._iku_kpoints

    @iku_kpoints.setter
    def iku_kpoints(self,value):
        if hasattr(self,"_red_kpoints"): delattr(self,"_red_kpoints")
        if hasattr(self,"_car_kpoints"): delattr(self,"_car_kpoints")
        self._iku_kpoints = value

    @property
    def nkpoints(self):
        return len(self.iku_kpoints)

    @property
    def nsym(self):
        return len(self.sym_car)

    @property
    def rlat(self):
        """calculate the reciprocal lattice"""
        if not hasattr(self,'_rlat'):
            self._rlat = rec_lat(self.lat)
        return self._rlat 

    @property
    def rlat_vol(self):
        return (2*np.pi)**3 * vol_lat(self.rlat)

    @property
    def lat_vol(self):
        return vol_lat(self.lat)

    @property
    def car_kpoints(self):
        """convert form internal yambo units to cartesian lattice units"""
        if not hasattr(self,"_car_kpoints"):
            self._car_kpoints = np.array([ k/self.alat for k in self.iku_kpoints ])
        return self._car_kpoints

    @property
    def red_kpoints(self):
        """convert from cartesian coordinates to reduced coordinates"""
        if not hasattr(self,"_red_kpoints"):
            self._red_kpoints = car_red(self.car_kpoints,self.rlat)
        return self._red_kpoints
  
    @property
    def sym_red(self):
        """Convert cartesian transformations to reduced transformations"""
        if not hasattr(self,"_sym_red"):
            sym_red = np.zeros([self.nsym,3,3],dtype=int)
            for n,s in enumerate(self.sym_car):
                sym_red[n] = np.round(np.dot(np.dot(self.lat,s.T),np.linalg.inv(self.lat)))
            self._sym_red = sym_red
        return self._sym_red

    @property
    def sym_rec_red(self):
        """Convert reduced transformations to reduced reciprocal transformations"""
        if not hasattr(self,"_sym_rec_red"):
            sym_rec_red = np.zeros([self.nsym,3,3],dtype=int)
            for n,s in enumerate(self.sym_red):
                sym_rec_red[n] = np.linalg.inv(s).T
            self._sym_rec_red = sym_rec_red
        return self._sym_rec_red
         
    @property
    def sym_rec(self):
        """Convert cartesian transformations to reciprocal transformations"""
        if not hasattr(self,"_sym_rec"):
            sym_rec = np.zeros([self.nsym,3,3])
            for n,s in enumerate(self.sym_car):
                sym_rec[n] = np.linalg.inv(s).T
            self._sym_rec = sym_rec
        return self._sym_rec

    @property
    def time_rev_list(self):
        """get a list of symmetries with time reversal"""
        time_rev_list = [False]*self.nsym
        for i in range(self.nsym):
            time_rev_list[i] = ( i >= self.nsym/(self.time_rev+1) )
        return time_rev_list

    def expand_kpoints(self,verbose=1,expand_mode=1,atol=1.e-6):
        """
        Take a list of qpoints and symmetry operations and return the full brillouin zone
        with the corresponding index in the irreducible brillouin zone

        :: expand_mode=0 : k' = S k
        :: expand_mode=1:  k' = S.T k  [default, consistent with yambo]
        """

        #check if the kpoints were already exapnded
        kpoints_indexes  = []
        kpoints_full     = []
        symmetry_indexes = []

        # Store original kpoints in iku coordinates
        self.ibz_kpoints = self.iku_kpoints

        #kpoints in the full brillouin zone organized per index
        kpoints_full_i = {}

        if expand_mode==0: sym_car_to_apply = self.sym_car
        if expand_mode==1:
            sym_car_to_apply = np.copy(self.sym_car)
            sym_car_to_apply = np.transpose(self.sym_car, (0, 2, 1))


        #expand using symmetries
        for nk,k in enumerate(self.car_kpoints):
            #if the index in not in the dicitonary add a list
            if nk not in kpoints_full_i:
                kpoints_full_i[nk] = []

            for ns,sym in enumerate(sym_car_to_apply):

                new_k = np.dot(sym,k)
                #check if the point is inside the bounds
                k_red = car_red([new_k],self.rlat)[0]
                k_red[np.abs(k_red) < atol] = 0. # Set to zero values < atol to avoid mistakes
                k_bz = (k_red+atol)%1

                #if the vector is not in the list of this index add it
                if not vec_in_list(k_bz,kpoints_full_i[nk],atol=atol):
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

        if verbose: print("%d kpoints expanded to %d"%(len(self.car_kpoints),len(kpoints_full)))

        #set the variables
        self.weights_ibz      = np.array(weights)
        self.kpoints_indexes  = np.array(kpoints_indexes)
        self.symmetry_indexes = np.array(symmetry_indexes)
        self.iku_kpoints      = [k*self.alat for k in kpoints_full]

    def get_path(self,path,debug=False):
        """
        Obtain a list of indexes and kpoints that belong to the regular mesh
        """
        if isinstance(path,Path):
            path = path.get_klist()

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
                for index, kpt in enumerate(self.car_kpoints):

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

    def get_units_info(self):

        info_string = \
        "          Yambo cartesian units [cc in yambo]: \n\
                ::   self.car_kpoints*2.*pi\n\
         \n\
          QE cartesian unists [cart. coord. in units 2pi/alat in QE: \n\
                ::   self.car_kpoints*self.alat[0]\n\
         \n\
          Internal yambo units [iku]: \n\
                ::   self.iku_kpoints\n\
         \n\
          Reduced coordinates [rlu in yambo, cryst. coord. in QE]: \n\
                ::   self.red_kpoints\n"
        print(info_string)

    def __str__(self):
        lines = []; app = lines.append
        app(marquee(self.__class__.__name__))
        app("lattice:")
        lines += [("%12.8lf " * 3) % tuple(vec) for vec in self.lat]
        app("atom positions:")
        for an, pos in zip(self.atomic_numbers, self.red_atomic_positions):
            app( "%3d " % an + ("%12.8lf " * 3) % tuple(pos) )
        if self.ibz_nkpoints!=self.nkpoints: app(f"{self.ibz_nkpoints} kpoints expanded to {self.nkpoints}")
        else: app(f"{self.nkpoints} kpoints in the IBZ")
        app(f"Time-reversal symmetry: {bool(self.time_rev)}")
        return "\n".join(lines)
