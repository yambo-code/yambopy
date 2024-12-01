#
# License-Identifier: GPL
#
# Copyright (C) 2024 The Yambo Team
#
# Authors: HPC, FP
#
# This file is part of the yambopy project
#
import numpy as np
from itertools import product

def calculate_distances(kpoints):   # Needs optimization
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
 
def expand_kpts(kpts,syms):
    """ 
    Fast expansion giving only coordinate list in output. See kpoints.py for a more complete version.

    Take a list of qpoints and symmetry operations and return the full brillouin zone
    with the corresponding index in the irreducible brillouin zone
    """
    nkpoints = len(kpts)
    nsym = len(syms)

    # Reshape kpts to (1, N, 3) to enable broadcasting
    kpts_reshaped = kpts[np.newaxis, :, :]
    transformed_kpts = np.einsum('ijk,jk->ij', syms, kpts_reshaped)
    full_kpts_flat = transformed_kpts.reshape(nsym * nkpoints, 3)
    nk_indices = np.repeat(np.arange(nkpoints), nsym)

    # Combine indices and transformed k-points into a structured array or list of tuples if needed
    full_kpts = np.column_stack((nk_indices, full_kpts_flat))

    return full_kpts

def vec_in_list(veca,vec_list,atol=1e-6):
    """
    Check if a vector exists in a list of vectors
    """
    if vec_list:
        vec_list = np.array(vec_list)

        # Use broadcasting to compare veca with all vectors in vec_list
        return np.all(np.isclose(veca, vec_list, atol=atol), axis=1).any()
    else:
        return False    # In case the vector is not in the list, otherwise the np.isclose() fails


def isbetween(a,b,c,eps=1e-5):
    """ Check if c is between a and b
    """
    return np.isclose(np.linalg.norm(a-c)+np.linalg.norm(b-c)-np.linalg.norm(a-b),0,atol=eps)

def red_car(red,lat):
    """
    Convert reduced coordinates to cartesian
    """
    return np.einsum('ij,ji->i', red, lat.T)

def car_red(car,lat):
    """
    Convert cartesian coordinates to reduced
    """

    return np.linalg.solve(np.array(lat).T, np.array(car).T).T

def vol_lat(lat):
    """
    Calculate the volume of a lattice
    """
    a1,a2,a3 = np.array(lat)
    return np.dot(a1,np.cross(a2,a3))

def rec_lat(lat):
    """
    Calculate the reciprocal lattice vectors
    """
    v = vol_lat(lat)
    a1,a2,a3 = np.array(lat)
    b1 = np.cross(a2,a3)/v
    b2 = np.cross(a3,a1)/v
    b3 = np.cross(a1,a2)/v
    return np.array([b1,b2,b3])

def replicate_red_kmesh(kmesh,repx=list(range(1)),repy=list(range(1)),repz=list(range(1))):
    """
    copy a kmesh in the tree directions
    the kmesh has to be in reduced coordinates
    """
    kmesh = np.array(kmesh)
    kmesh_nkpoints = len(kmesh)

    kmesh_full = []
    kmesh_idx  = []
    for x,y,z in product(repx,repy,repz):
        kmesh_shift = kmesh + np.array([x,y,z])
        kmesh_full.append(kmesh_shift)
        kmesh_idx.append(list(range(kmesh_nkpoints)))

    return np.vstack(kmesh_full), np.hstack(kmesh_idx)


def point_matching(a,b,double_check=True,debug=False,eps=1e-8):                 #######
    """
    Matches the points of list a to the points of list b
    using a nearest neighbour finding algorithm

    Arguments:

        double_check: after the nearest neighbours are assigned check further
        if the distance between points is within the precision eps

        eps: precision for the double check (default: 1e-8)

    """
    #karma
    from scipy.spatial import cKDTree
    from time import time
    a = np.array(a)
    b = np.array(b)
    start_time = time()

    #initialize the kdtree
    kdtree = cKDTree(a, leafsize=10)
    map_b_to_a = []
    for xb in b:
        current_dist,index = kdtree.query(xb, k=1, distance_upper_bound=6)
        map_b_to_a.append(index)
    map_b_to_a = np.array(map_b_to_a)
    
    if debug:
        print("took %4.2lfs"%(time()-start_time))

    if double_check:
        for ib,ia in enumerate(map_b_to_a):
            dist = np.linalg.norm(a[ia]-b[ib])
            if dist > eps:
                raise ValueError('point a %d: %s is far away from points b %d: %s  dist: %lf'%(ia,str(a[ia]),ib,str(b[ib]),dist))

    return map_b_to_a

def bravais_types(lats,alat_0):
    """
    Determine Bravais lattice type of unit cell

    :: lats -> lattice vectors from YamboLatticeDB corresponding to lat
    :: alat_0 -> lattice parameter from YamboLatticeDB corresponding to alat[0]

    More lattice types to be implemented
    """
    from math import sqrt,cos,sin

    bravais_types = ['Hexagonal and Trigonal P','Orthorhombic P']

    lats_ = lats/alat_0
 
    if np.array_equal(lats_[0],[1.,0.,0.]):
        
        if np.allclose(lats_[1],[-0.5,sqrt(3.)/2.,0.]):
        
            if np.allclose(lats_[2],[0.,0.,lats_[2,2] ]): return bravais_types[0]

        if np.array_equal(lats_[1],[0.,lats_[1,1],0.]):
         
            if np.allclose(lats_[2],[0.,0.,lats_[2,2] ]): return bravais_types[1]
    else: return 'No type'
