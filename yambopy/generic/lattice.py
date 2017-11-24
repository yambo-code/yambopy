#
# Copyright (c) 2017, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
import numpy as np
import itertools
import time
import collections
from itertools import product
from scipy import spatial
from sys import stdout

def load_bar(x,n):
    """
    A Simple loadbar function
    """
    bar_length = 60
    x+=1
    ratio = x/float(n)
    c = int(ratio * bar_length)
    stdout.write("["+"="*c+" "*(bar_length-c)+"] %03.3f%%" % (ratio*100))
    if (x==n): stdout.write("\n")
    stdout.flush()
    stdout.write("\r")

def red_car(red,lat):
    """
    Convert reduced coordinates to cartesian
    """
    lat = np.array(lat)
    return np.array([coord[0]*lat[0]+coord[1]*lat[1]+coord[2]*lat[2] for coord in red])

def car_red(car,lat):
    """
    Convert cartesian coordinates to reduced
    """
    return np.array([np.linalg.solve(np.array(lat).T,coord) for coord in np.array(car)])

def isbetween(a,b,c):
    """check if c is between a and b"""
    return np.isclose(np.linalg.norm(a-c)+np.linalg.norm(b-c)-np.linalg.norm(a-b),0)
    
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

class Lattice():
    """ 
    Define lattice parameters
    """
    def __init__(self, lat, positions=None, numbers=None, red_kpoints=None, car_kpoints=None, weights=None):
        self.lat  = lat
        self.rlat = rec_lat(lat)
        
        if red_kpoints is not None and car_kpoints is not None:
            raise ValueError('Both red_kpoints and car_kpoints are specified at the same time, which does not make sense')
        
        if red_kpoints is None and car_kpoints is not None:
            red_kpoints = car_red(car_kpoints,self.rlat)
            self.nkpoints = len(red_kpoints)
        if car_kpoints is None and red_kpoints is not None: 
            car_kpoints = red_car(red_kpoints,self.rlat)
            self.nkpoints = len(car_kpoints)
            
        self.red_kpoints = red_kpoints #kpoints in reduced coordinates
        self.car_kpoints = car_kpoints #kpoints in cartesian coordinates
        self.positions = positions
        self.numbers = numbers       

    def set_kmesh(self,kmesh):
        """
        set the kpoints to a regular kmesh
        """
        from spglib import get_ir_reciprocal_mesh
        
        cell = (self.lat, self.positions, self.numbers)
        mapping, grid = get_ir_reciprocal_mesh(kmesh, cell, is_shift=[0, 0, 0])

        #number of kpoints
        nkpoints     = np.prod(kmesh)
        nkpoints_ibz = len(np.unique(mapping))

        #counter counts the number of occurences of element in a list
        weights_ibz = np.zeros([nkpoints_ibz])
        counter = collections.Counter(mapping)
        for nk_ibz,inv_weight in enumerate(counter.values()):
            weights_ibz[nk_ibz] = float(inv_weight)/nkpoints

        #irreducible kpoints
        irr_kpoints  = (grid[list(counter.keys())] + [0.0, 0.0, 0.0]) / kmesh
        full_kpoints = (grid + [0.0, 0.0, 0.0]) / kmesh
        mapping_to_ibz = dict(list(zip(list(counter.keys()),list(range(nkpoints_ibz)))))
        ibz_to_full = [mapping_to_ibz[i] for i in mapping]
    
        self.red_kpoints = full_kpoints
        self.car_kpoints = red_car(full_kpoints,self.rlat)
        self.red_kpoints_ibz = irr_kpoints
        self.car_kpoints_ibz = red_car(irr_kpoints,self.rlat)
        self.weights_ibz = weights_ibz
        self.mapping_to_ibz = mapping_to_ibz
        self.ibz_to_full = ibz_to_full
        self.mapping = mapping
    
    def get_wigner(self,mesh,eps=1e-7):
        """
        Get wigner seitz cell
        """
        from math import sqrt
        nx,ny,nz = mesh
        alat = [np.linalg.norm(v) for v in self.lat]
        lat = self.lat/alat
        lat = np.array(lat)
        lat = lat*np.array([1.0/nx,1.0/ny,1.0/nz])
        
        #generate new wigner seitz
        xiter = range(-nx,nx)
        yiter = range(-ny,ny)
        ziter = range(-nz,nz)
        if nx == 1: xiter = range(1)
        if ny == 1: yiter = range(1)
        if nz == 1: ziter = range(1)
                    
        points = np.array([(x,y,z) for x,y,z in product(xiter,yiter,ziter)])
        points = red_car(points,lat)

        #create a list of 26 points corresponding to the neighbouring cells
        neighbours = []
        for x,y,z in product(range(-1,2),repeat=3):
            #do not add zero
            if x==0 and y==0 and x==0:
                continue
            neighbours.append([x*nx,y*ny,z*nz])
        #convert to cartesian cordinates
        neighbours = red_car(np.array(neighbours),lat)

        #for each of the points
        filtered = []
        degeneracies = []
        for point in points:
            #calculate the distance with all the neighbours
            dists = np.array([ np.linalg.norm(neighbour-point) for neighbour in neighbours])
            #get the distance to the center
            dist0 = np.linalg.norm(point)
            # if the distance to the center is smaller then add the point
            if (dist0 <= dists+eps).all():
                filtered.append(point)
                # the degeneracy is obtained checking for how many points are similar
                degeneracies.append( sum(np.isclose(dist0,dists,atol=eps))+1 )

        filtered = np.array(filtered)
        filtered = np.array(np.rint(car_red(filtered,lat)),dtype=int)
        degeneracies = np.array(degeneracies)
        self.filtered = filtered
        self.degeneracies = degeneracies
        return filtered, degeneracies
    
    def _plotBZ(self,ax,kpts,values,dim=1.2,size=50,cmap=None):
        """
        Plot a list of values in the brillouin zone
        ax -> matplotlib axis
        """
        import matplotlib.pyplot as plt
        if cmap is None:
            cmap = plt.get_cmap("viridis")
        
        ax.scatter(kpts[:,0], kpts[:,1], cmap=cmap, c=values, marker="H", s=size, lw=0, rasterized=True)
        
        plt.xlim([-dim,dim])
        plt.ylim([-dim,dim])    
        
        ax.set_aspect('equal', 'datalim')
        ax.xaxis.set_major_locator(plt.NullLocator())
        ax.yaxis.set_major_locator(plt.NullLocator())
    
    def get_path(self,path,kpts=None):
        """ Obtain a list of indexes and kpoints that belong to the regular mesh
        """
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
                for index, kpt in enumerate(kpts):

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
                print(("%12.8lf "*3)%tuple(kpt), index)

        self.bands_kpoints = bands_kpoints
        self.bands_indexes = bands_indexes
        self.bands_highsym_qpts = path_car

        return bands_kpoints, bands_indexes, path_car
                
    def get_double_grid(self,lattice,debug=False,delta_shift=[1e-4,1e-4,1e-4]):
            """
            Generate a double_grid_map array that maps the kpoints from lat to the kpoints of this class
            Arguments:
            lat - lattice instance with the double grid
            """
            print("double grid:", lattice.full_nkpoints)
            print("coarse grid:", self.full_nkpoints)
            
            #repeat the coarse kpoints in three directions
            coarse_kpoints = []
            coarse_indexes = []
            for i,j,k in itertools.product(range(-1,2),repeat=3):
                disp = red_car([[i,j,k]],self.rlat)[0]
                coarse_kpoints.append( self.car_kpoints + disp )
                coarse_indexes.append( range(self.full_nkpoints) )
                
            #shift the coarse kpoints by a infinitesimal
            coarse_kpoints = np.vstack(coarse_kpoints)+delta_shift
            coarse_indexes = np.hstack(coarse_indexes)
            
            if debug:
                import matplotlib.pyplot as plt
                #double grid kpoints
                dgk = lattice.car_kpoints
                plt.scatter(dgk[:,0],dgk[:,1],c='r',s=50,lw=0)
                #coarse kpoints
                #dgk = self.car_kpoints
                dgk = coarse_kpoints
                plt.scatter(dgk[:,0],dgk[:,1],c='b',s=25,lw=0)
                plt.show()
            
            delta_shift = np.array(delta_shift)    
            double_to_coarse_grid = [0]*lattice.full_nkpoints                   #many to one
            coarse_to_double_grid = [ [] for nk in range(self.full_nkpoints) ] #one to many

            print("finding neighbour k-points:")
            start_time = time.time()
        
            #iterate over the points of the double grid
            for dg_nk, dg_k in enumerate(lattice.car_kpoints):
                load_bar(dg_nk,lattice.full_nkpoints)
                
                #the closer k-point from the coarse grid to this k-point gives the index
                dist = spatial.distance.cdist([lattice.car_kpoints[dg_nk]],coarse_kpoints)[0]
                coarse_nk = np.argmin(dist)
                coarse_nk = coarse_indexes[coarse_nk]
                double_to_coarse_grid[dg_nk] = coarse_nk
                coarse_to_double_grid[coarse_nk].append(dg_nk)
                
            print("took %4.2lfs" % (time.time()-start_time))
            return double_to_coarse_grid, coarse_to_double_grid
                    
    def __str__(self):
        s =  ""
        s += "lattice:\n"
        for r in self.lat:
            s += ("%12.8lf "*3+"\n")%tuple(r)
            
        return s
