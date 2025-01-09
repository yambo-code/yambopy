#
# License-Identifier: GPL
#
# Copyright (C) 2024 The Yambo Team
#
# Authors: HPC
#
# This file is part of the yambopy project
#
import numpy as np

def calculate_distances(kpoints):
    """
    take a list of k-points and calculate the distances between all of them
    """
    kpoints = np.array(kpoints)
    distances = [0]
    distance = 0
    for nk in range(1,len(kpoints)):
        distance += np.linalg.norm(kpoints[nk-1]-kpoints[nk])
        distances.append(distance)   
    return distances

class Path(object):
    """ Class that defines a path in the brillouin zone
       
        Input format:
        :: Path(klist,intervals)

        :: klist = [ [[ k1_x, k1_y, k1_z ], 'K1_label'],
                        ...                 ...
                     [[ kN_x, kN_y, kN_z ], 'KN_label'] ]
        
        :: intervals = [ n_steps_1, n_steps_2, ... , n_steps_N-1 ]

        NB: interval steps are the number of kpoints along each high-symmetry direction, must be provided but is USED ONLY FOR INTERPOLATIONS
    """
    def __init__(self,klist,intervals):
        """
        Generation of a path in reciprocal space by specifying a list of k-points
        """

        self.intervals = intervals
        if len(intervals)>len(klist)-1:
            print(f"[WARNING] number of intervals is {len(intervals)}, should be {len(klist)-1} (no. of path lines). Taking the first {len(klist)-1} intervals.")
            self.intervals = intervals[:len(klist)-1]
        if len(intervals)<len(klist)-1:
            Nsteps=20
            print(f"[WARNING] number of intervals (path lines) inconsistent with special points specified, using default of {Nsteps} steps per direction")
            self.intervals = [ Nsteps for ik in range(len(klist)-1) ]

        klabels = []
        kpoints = []
        for kpoint, klabel in klist:
            kpoints.append(kpoint)
            klabels.append(klabel)
        self.kpoints = np.array(kpoints)
        self.klabels = klabels
    
    def as_dict(self):
        d = {'kpoints':self.kpoints.tolist(),
             'klabels':self.klabels,
             'intervals':self.intervals}
        return d

    @classmethod
    def from_dict(cls,d):
        klist = zip(d['kpoints'],d['klabels'])
        return cls(klist,d['intervals'])

    @property
    def distances(self):
        if hasattr(self,'_distances'): return self._distances
        k0 = np.array(self.kpoints[0])
        dist = 0
        distances = [0]
        for kpt in np.array(self.kpoints[1:]):
            dist += np.linalg.norm(kpt-k0)
            k0 = kpt
            distances.append(dist)
        self._distances = distances
        return distances

    def set_xticks(self,ax):
        ax.set_xticks(self.distances)
        ax.set_xticklabels(self.klabels)

    def __iter__(self):
        return iter(zip(self.kpoints,self.klabels,self.distances))

    def get_klist(self):
        """ 
        Output in the format of quantum espresso == [ [kx, ky, kz, 1], ... ]
        """
        kpoints = self.kpoints
        intervals = self.intervals
        kout  = np.zeros([sum(intervals)+1,4])
        kout[:,3] = 1
        io = 0
        for ik,interval in enumerate(intervals):
          for ip in range(interval):
            kout[io,:3] = kpoints[ik] + float(ip)/interval*(kpoints[ik+1] - kpoints[ik])
            io = io + 1
        kout[io,:3] = kpoints[ik] + float(ip+1)/interval*(kpoints[ik+1] - kpoints[ik])

        return kout

    def get_indexes(self):
        """ get the index of each point of the path
        """

        indexes = []
        index = 0
        for n,label in enumerate(self.intervals):
            indexes.append([index,self.klabels[n]])
            index += self.intervals[n] 
        indexes.append([index,self.klabels[-1]])
        return indexes

def vec_in_list(veca,vec_list,atol=1e-6):
    """ check if a vector exists in a list of vectors
    """
    return np.array([ np.allclose(veca,vecb,rtol=atol,atol=atol) for vecb in vec_list ]).any()

def red_car(red,lat):
    """
    Convert reduced coordinates to cartesian
    """
    lat = np.array(lat)
    red = np.array(red)
    if lat.shape != (3,3):
        raise ValueError('Wrong lattice dimensions expected (3,3) got {}'.format(lat.shape))
    return np.array([coord[0]*lat[0]+coord[1]*lat[1]+coord[2]*lat[2] for coord in red])

def car_red(car,lat):
    """
    Convert cartesian coordinates to reduced
    """
    car = np.array(car)
    lat = np.array(lat)
    if lat.shape != (3,3):
        raise ValueError('Wrong lattice dimensions expected (3,3) got {}'.format(lat.shape))
    return np.array([np.linalg.solve(np.array(lat).T,coord) for coord in car])

def rec_lat(lat):
    """
    Calculate the reciprocal lattice vectors
    """
    a1,a2,a3 = np.array(lat)
    v = np.dot(a1,np.cross(a2,a3))
    b1 = np.cross(a2,a3)/v
    b2 = np.cross(a3,a1)/v
    b3 = np.cross(a1,a2)/v
    return np.array([b1,b2,b3])
