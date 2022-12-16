# Copyright (c) 2018, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
import os
import numpy as np
from itertools import product
from netCDF4 import Dataset
from yambopy.plot.plotting import add_fig_kwargs
from yambopy.plot import *
from yambopy.tools.string import marquee
from yambopy.lattice import isbetween, car_red, red_car, rec_lat, vol_lat
from yambopy.units import ha2ev

max_exp = 50
atol = 1e-3

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
    def __init__(self,atomic_numbers,car_atomic_positions,eigenvalues,sym_car,kpts_iku,
                 lat,alat,temperature,electrons,spin,time_rev,spinor):

        self.atomic_numbers       = atomic_numbers   
        self.car_atomic_positions = car_atomic_positions
        self.eigenvalues          = eigenvalues     
        self.sym_car              = sym_car         
        self.kpts_iku             = kpts_iku        
        self.lat                  = lat             
        self.alat                 = alat            
        self.temperature          = temperature     
        self.electrons            = electrons       
        self.spin                 = spin            
        self.time_rev             = time_rev
        self.spinor               = spinor

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
            
            natoms_a = database.variables['N_ATOMS'][:].astype(int).T
            tmp_an = database.variables['atomic_numbers'][:].astype(int)
            tmp_apos = database.variables['ATOM_POS'][:,:]

            flatten = lambda l: [item for sublist in l for item in sublist]
            atomic_numbers = flatten([[tmp_an[n]]*na for n,na in enumerate(natoms_a)])
            atomic_positions = np.vstack([[tmp_apos[n,ia] for ia in range(na)] for n,na in enumerate(natoms_a) ])
            # I change the shape of eigenvalues
            # From now on is [spin_index, kpoint_index, band_index]
            # This must be change all along the code
            args = dict( atomic_numbers       = atomic_numbers,
                         car_atomic_positions = atomic_positions,
                         eigenvalues          = database.variables['EIGENVALUES'][:,:]*ha2ev,
                         sym_car              = database.variables['SYMMETRY'][:],
                         kpts_iku             = database.variables['K-POINTS'][:].T,
                         lat                  = database.variables['LATTICE_VECTORS'][:].T,
                         alat                 = database.variables['LATTICE_PARAMETER'][:].T,
                         temperature          = dimensions[13],
                         electrons            = dimensions[14],
                         spin                 = int(dimensions[11]),
                         time_rev             = dimensions[9],
                         spinor               = database.variables['EIGENVALUES'].shape[0],
                         )

        return cls(**args)

    @property
    def red_atomic_positions(self):
        return car_red(self.car_atomic_positions,self.lat)

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
        _,_,nbands = self.eigenvalues.shape
        return nbands

    @property
    def nbandsv(self):
        return int(self.electrons/self.spin_degen)

    @property
    def nbandsc(self):
        return int(self.nbands - self.nbandsv)

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
        sym_rec = np.zeros([self.nsym,3,3])
        for n,s in enumerate(self.sym_car):
            sym_rec[n] = np.linalg.inv(s).T
        return sym_rec

    @property
    def efermi(self):
        if not hasattr(self,"_efermi"):
            # break here??
            # I have changed get_efermi by fermi and now it works. To be check by the author
            self._efermi = self.get_fermi
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
            if self.spinor == 1:
               return sum([sum(self.spin_degen*fermi_array(self.eigenvalues[0,nk],ef))*self.weights[nk] for nk in range(self.nkpoints)])-self.electrons
            elif self.spinor == 2:
               sum_up = sum([sum(self.spin_degen*fermi_array(self.eigenvalues[0,nk],ef))*self.weights[nk] for nk in range(self.nkpoints)]) 
               sum_dw = sum([sum(self.spin_degen*fermi_array(self.eigenvalues[1,nk],ef))*self.weights[nk] for nk in range(self.nkpoints)]) 
               return sum_up + sum_dw -self.electrons

        efermi = bisect(occupation_minus_ne,self.min_eival,self.max_eival)

        if verbose: print("fermi: %lf eV"%efermi)

        self.eigenvalues -= efermi

        self.occupations = np.zeros([self.spinor,self.nkpoints,self.nbands],dtype=np.float32)
        for nspin in range(self.spinor):
            for nk in range(self.nkpoints):
                self.occupations[nspin,nk] = fermi_array(self.eigenvalues[nspin,nk,:self.nbands],0)

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

        # bug
        #
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

    def plot_bs_ax(self,ax,path,bandmin=None,bandmax=None,add_indexes=False,**kwargs):
        """
        Plot this bandstructure on Matpltolib ax
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
            ax.axvline(distance,color='k')
            distance += np.linalg.norm(bands_highsym_qpts[nk]-bands_highsym_qpts[nk-1])
            bands_highsym_qpts_distances.append(distance)
        ax.axvline(distance,color='k')

        #plot bands
        if self.spinor == 2:
           color = kwargs.pop('c','red')
           ax.plot(bands_distances,self.eigenvalues[0,bands_indexes,bandmin:bandmax],c=color,**kwargs)
           color = kwargs.pop('c','blue')
           ax.plot(bands_distances,self.eigenvalues[1,bands_indexes,bandmin:bandmax],c=color,**kwargs)
        ax.set_xlim(0,max(bands_distances))

        if add_indexes:
            ax.set_xticks(bands_distances)
            ax.set_xticklabels(np.array(bands_indexes)+1)
            for d in bands_distances:
                ax.axvline(d,color='k',alpha=0.5)

        return ax

    @add_fig_kwargs
    def plot_bs(self,path,**kwargs):
        """ Plot the difference in energies of two bands
        """
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        self.plot_bs_ax(ax,path,**kwargs)
        return fig

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
        lines = []; app = lines.append
        app(marquee(self.__class__.__name__))
        app("reciprocal lattice:")
        app("\n".join([("%12.8lf "*3)%tuple(r) for r in self.rlat]))
        app("lattice:")
        app("\n".join([("%12.8lf "*3)%tuple(r) for r in self.lat]))
        app("alat:")
        app(("%12.8lf "*3)%tuple(self.alat))
        app("atom positions:")
        for an, pos in zip(self.atomic_numbers, self.red_atomic_positions):
            app( "%3d " % an + ("%12.8lf " * 3) % tuple(pos) )
        app("nkpoints: %d"%self.nkpoints)
        app("symmetry operations: %d\n"%len(self.sym_car))
        app("temperature : %lf"%self.temperature)
        app("electrons   : %lf"%self.electrons)
        return "\n".join(lines)
