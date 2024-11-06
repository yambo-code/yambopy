#
# License-Identifier: GPL
#
# Copyright (C) 2024 The Yambo Team
#
# Authors: FP, HPC
#
# This file is part of the yambopy project
#
import numpy as np
from itertools import product
from yambopy.lattice import red_car, vec_in_list, isbetween, car_red
from qepy.lattice import Path

def expand_kpoints(car_kpoints,sym_car,rlat,atol=1.e-6):
    """
    Take a list of kpoints and symmetry operations and return the full brillouin zone
    with the corresponding index in the irreducible brillouin zone

    Input:
    * IBZ kpoints in Cartesian coodinates
    * Symmetry operations in Cartesian coordinates
    * rlat: reciprocal lattice vectors
    * atol: tolerance for the distance between kpoints

    Output:
    * weights: weights of the kpoints in the irreducible brillouin zone
    * kpoints_indexes: indexes of the kpoints in the irreducible brillouin zone
    * symmetry_indexes: indexes of the symmetries used to generate the kpoints
    * kpoints_full: kpoints in the full brillouin zone
    """

    #check if the kpoints were already exapnded
    kpoints_indexes  = []
    kpoints_full     = []
    symmetry_indexes = []

    #kpoints in the full brillouin zone organized per index
    kpoints_full_i = {}

    #expand using symmetries
    for nk,k in enumerate(car_kpoints):
        #if the index in not in the dicitonary add a list
        if nk not in kpoints_full_i:
            kpoints_full_i[nk] = []

        for ns,sym in enumerate(sym_car):

            new_k = np.dot(sym,k)
            #check if the point is inside the bounds
            k_red = car_red([new_k],rlat)[0]
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

    return np.array(weights), np.array(kpoints_indexes), np.array(symmetry_indexes), np.array(kpoints_full)

def get_path_car(kpts_path_car,path):
    """
    Convert a Path object in reduced coordinates to a path in Cartesian coordinates
    """
    return Path( [[kpts_path_car[i],path.klabels[i]] for i in range(len(kpts_path_car))],path.intervals )

def get_path(car_kpoints,rlat,sym_car,path,debug=False):
    """
    Obtain a list kpoints along a specific high-symmetry path

    Input:
    * car_kpoints: kpoints in Cartesian coordinates [IBZ/Full BZ]
    * rlat: reciprocal lattice vectors
    * sym_car [symmetry ops. if car_kpoints given in IBZ] | None [if car_kpoints in fulll BZ]
    * path: Path object with the high-symmetry path (points in reduced coordinates)

    Output:
    * bands_kpoints: kpoints Cartesian coordinates along the path
    * bands_indexes: indexes of the kpoints in the path
    * path_car: path in Cartesian coordinates
    """

    # expand if symmetries are provided, otherwise the kpoints are considered already expanded
    if sym_car is None: nks = list(range(len(car_kpoints)))
    else:               _, nks, _, car_kpoints = expand_kpoints(car_kpoints,sym_car,rlat)

    # high-symmetry points in cartesian coordinates
    kpts_path_car = red_car(path.kpoints, rlat)
    # Path object in cartesian coordinates (for later plotting)
    path_car = get_path_car(kpts_path_car,path)

    #find the points along the high symmetry lines
    bands_kpoints = []
    bands_indexes = []

    #for all the paths
    for k in range(len(path.kpoints)-1):

        # store here all the points in the path
        # key:   has the coordinates of the kpoint rounded to 4 decimal places
        # value: index of the kpoint
        #        distance to the starting kpoint
        #        the kpoint cordinate
        kpoints_in_path = {}

        start_kpt = kpts_path_car[k]   #start point of the path
        end_kpt   = kpts_path_car[k+1] #end point of the path

        #generate repetitions of the brillouin zone
        for x,y,z in product(list(range(-1,2)),list(range(-1,2)),list(range(1))):

            #shift the brillouin zone
            shift = red_car([np.array([x,y,z])],rlat)[0]

            #iterate over all the kpoints
            for index, kpt in zip(nks,car_kpoints):

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

    return bands_kpoints, bands_indexes, path_car