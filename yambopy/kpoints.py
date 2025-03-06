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
from scipy.spatial import KDTree

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


def make_kpositive(klist, tol=1e-6):
    """
    Shifts all k-points into the range [0,1) by applying modulo operation.

    Parameters
    ----------
    klist : np.ndarray
        Array of k-points in crystal coordinates.
    tol : float, optional
        A small tolerance to ensure numerical stability, default is 1e-6.

    Returns
    -------
    np.ndarray
        Array of k-points mapped to the range [0,1).
    """
    kpos = klist - np.floor(klist)  # Ensure k-points are within [0,1)
    return (kpos + tol) % 1  # Apply small tolerance correction


def build_ktree(kpts):
    """
    Builds a k-d tree for efficient k-point searching.

    Parameters
    ----------
    kpts : np.ndarray
        Array of k-points in crystal coordinates.

    Returns
    -------
    KDTree
        A KDTree structure for fast nearest-neighbor lookup of k-points.
    """
    tree = make_kpositive(kpts)  # Normalize k-points to [0,1)
    return KDTree(tree, boxsize=[1, 1, 1])  # Construct KDTree with periodic boundaries


def find_kpt(tree, kpt_search, tol=1e-5):
    """
    Finds the index of a k-point in a KDTree within a given tolerance.

    Parameters
    ----------
    tree : KDTree
        Pre-built KDTree of k-points.
    kpt_search : np.ndarray
        The k-point to search for in crystal coordinates.
    tol : float, optional
        Tolerance for k-point matching, default is 1e-5.

    Returns
    -------
    int
        Index of the matched k-point in the original k-point list.

    Raises
    ------
    SystemExit
        If the k-point is not found within the specified tolerance.
    """
    kpt_search = make_kpositive(kpt_search)  # Normalize k-point
    dist, idx = tree.query(kpt_search, workers=1)  # Perform nearest-neighbor search
    assert np.max(dist) < tol, "Kpoint not found"
    return idx  # Return the index of the found k-point

def regular_grid(nk1,nk2,nk3):
    """
    Generation of positive-coordinate full 'square' grid starting from zero
    Shifted grids not allowed for now
    """
    i, j, k = np.meshgrid(np.arange(nk1), np.arange(nk2), np.arange(nk3), indexing='ij')
    xkg = np.array([
        i.flatten() / nk1,
        j.flatten() / nk2,
        k.flatten() / nk3,
    ])
    return xkg


def kpoint_grid(nk1,nk2,nk3,sym_and_trev,IBZ=True,eps=1.0e-5):
    """
    Generation of gamma-centered Monkhorst-Pack grid.

    This function is the python porting of the Quantum ESPRESSO subroutine `kpoint_grid` found in PW/kpoint_grid.f90

    Input:
        nk1, nk2, nk3 -> grid dimensions
        sym_and_trev  -> (sym_red, time_rev, time_rev_list)

        (Symmetries in reduced coordinates, trev logical, indices of symmetries composed with time reversal)
        These can be obtained from the attributes of the same name in YamboLatticeDB

        IBZ [True]    -> grid is reduced to IBZ by checking equivalent points
            [False]   -> grid is generated with symmetries turned off
                         (i.e. NOT in BZ - the Wigner-Seitz cell - but in the
                          primitive reciprocal unit cell)
        epes [1e-5]   -> numerical error for equivalent points

        In order to obtain a list of kpoints in the full BZ, use IBZ==True and then
        expand_kpoints().

    Output:
        xks -> Number of kpoints in IBZ
        xk  -> kpoints in IBZ in CRYSTAL coordinates (convert to cc with red_car())
        wk  -> kpoints weights
    """
    Nk    = nk1*nk2*nk3
    wkk   = np.ones(Nk)    # Weights
    equiv = np.arange(Nk)  # Initially, each kpoint is equivalent to itself

    # symmetry and time-reversal info
    sym_red, is_trev, is_sym_trev = sym_and_trev
    Nsym = len(sym_red)

    # Generate full regular grid in crystal coordinates
    xkg = regular_grid(nk1,nk2,nk3)

    if IBZ:
        # Now we have to start checking for equivalent points:
        for ik in range(Nk):
            # Check if this k-point has already been found equivalent to another
            # Do not consider points previously found equivalent to another
            # Check both k and -k
            if equiv[ik] == ik:
                for i_s in range(Nsym):
                    # Apply symmetry operations
                    xkr = np.dot(sym_red[i_s,:,:], xkg[:,ik])
                    # Bring back in 1st BZ
                    xkr -= np.round(xkr)
                    # Take opposite if symmetry is composed with TR
                    if is_sym_trev[i_s]: xkr = -xkr

                    # Check if in the list
                    xx, yy, zz = xkr*[nk1,nk2,nk3]
                    in_the_list = all( abs(v-round(v)) <= eps for v in [xx, yy, zz])

                    if in_the_list:
                        i,j,k = [(round(xkr[dim]*nki + 2*nki) % nki) + 1 for dim, nki in enumerate([nk1, nk2, nk3])]
                        n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3

                        if n > ik and equiv[n] == n:
                            equiv[n] = ik
                            wkk[ik] += 1.0
                        elif equiv[n] != ik or n < ik:
                            raise ValueError("Error in the checking algorithm")

                    # Time reversal symmetry check
                    if is_trev:
                        xx, yy, zz = -xkr*[nk1,nk2,nk3]
                        in_the_list = all(abs(v - round(v)) <= eps for v in [xx, yy, zz])

                        if in_the_list:
                            i, j, k = [(round(-xkr[dim]*nki + 2*nki) % nki) + 1 for dim, nki in enumerate([nk1, nk2, nk3])]
                            n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3

                            if n > ik and equiv[n] == n:
                                equiv[n] = ik
                                wkk[ik] += 1.0
                            elif equiv[n] != ik or n < ik:
                                raise ValueError("Error in the checking algorithm")

    # Filter unique k-points
    unique_kpoints = (equiv == np.arange(Nk))
    # Bring back unique points into first BZ
    xk = xkg[:,unique_kpoints] - np.round(xkg[:,unique_kpoints])
    wk = wkk[unique_kpoints] / np.sum(wkk[unique_kpoints])  # Normalize weights
    nks = len(xk[0])

    return nks, xk.T, wk
