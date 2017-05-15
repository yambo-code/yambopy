# Copyright (c) 2016, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
from yambopy import *
from yambopy.plot import *
import os

ha2ev = 27.211396132

def isbetween(a,b,c):
    """
    Check if cartesian point c is between point a and b
    """
    return np.isclose(np.linalg.norm(a-c)+np.linalg.norm(b-c)-np.linalg.norm(a-b),0,rtol=1e-05, atol=1e-06)

class YamboRTDB():
    """
    Open the RT databases and store it in a RTDB class
    """
    def __init__(self,folder='.',calc='.',save=None,referencedb='ndb.RT_reference_components',carriersdb='ndb.RT_carriers'):
        # Find path with RT data
        # Yambopy's realtime scripts folder-structure
        if os.path.exists('%s/%s/pulse/%s'%(folder,calc,referencedb)):
            self.path = '%s/%s/pulse'%(folder,calc)
        # Custom path
        elif os.path.exists('%s/%s/%s'%(folder,calc,referencedb)):
            self.path = '%s/%s'%(folder,calc)
        else:
            raise ValueError('Cannot find file %s in %s/%s'%(referencedb,folder,calc))

        # Set save path
        if save==None:
            if os.path.exists('%s/SAVE'%folder):
                self.save = '%s/SAVE'%folder
            else:
                raise ValueError('Cannot find SAVE in folder %s'%folder)
        else:
            if os.path.exists(save):
                self.save = save
            else:
                raise ValueError('Cannot find SAVE in folder %s'%save)

        self.referencedb = referencedb
        self.carriersdb = carriersdb

        #read save for symmetries
        try:
            filename = '%s/ns.db1'%self.save
            database    = Dataset(filename)
        except:
            raise ValueError( "Error reading %s database"%filename )
        self.alat             = database.variables['LATTICE_PARAMETER'][:].T
        self.lat              = database.variables['LATTICE_VECTORS'][:].T
        self.sym_car          = database.variables['SYMMETRY'][:]
        dimensions = database.variables['DIMENSIONS'][:]
        self.time_rev = dimensions[9]
        database.close()

        #read reference database
        db = Dataset("%s/%s"%(self.path,self.referencedb))
        self.nband_min, self.nband_max, self.nkpoints = db['RT_vars'][:].astype(int)
        self.nbands = self.nband_max - self.nband_min + 1
        db.close()

        #get energies of bands
        db = Dataset("%s/%s"%(self.path,carriersdb))
        self.eigenvalues = db['RT_carriers_E_bare'][:].reshape([self.nkpoints,self.nbands])*ha2ev

        #get kpoints coordinates
        self.kpts_iku = db['RT_kpt'][:].T

        db.close()

        #get a list of symmetries with time reversal
        nsym = len(self.sym_car)

        #caclulate the reciprocal lattice
        self.rlat  = rec_lat(self.lat)
        self.nsym  = len(self.sym_car)

        #convert form internal yambo units to cartesian lattice units
        self.kpts_car = np.array([ k/self.alat for k in self.kpts_iku ])
        #convert cartesian transformations to reduced transformations
        inv = np.linalg.inv
        self.sym_rlu = np.zeros([self.nsym,3,3])
        for n,s in enumerate(self.sym_car):
            a = np.dot(s.T,inv(self.rlat))
            self.sym_rlu[n] = np.dot(inv(self.lat.T),a)

        #convert cartesian transformations to reciprocal transformations
        self.sym_rec = np.zeros([self.nsym,3,3])
        for n,s in enumerate(self.sym_car):
            self.sym_rec[n] = inv(s).T

        #read the databases
        self.readDB()

        #integrate the occupations
        self.integrate()

        #status
        self.expanded = False

    def readDB(self):
        """
        """

        #get how many rt databases exist
        files = [ filename for filename in  os.listdir(self.path) if 'ndb.RT_carriers_Time' in filename]
        print "Number of RT carrier files:", len(files)

        # sorting
        units = {'as':1e-18,'fs':1e-15,'ps':1e-12}
        s = []
        for filename in files:
            for unit in units.keys():
                if unit in filename:
                    factor = units[unit]
            s.append((float(re.findall("\d+\.\d+", filename)[0])*factor,filename))
        ordered_files=sorted(s)
        self.ntimes = len(ordered_files)

        #read all of them
        self.RT_carriers_delta_f        = np.zeros([self.ntimes,self.nkpoints,self.nbands])
        #self.RT_carriers_dE_Self_Energy = np.zeros([self.ntimes,self.nbands,self.nkpoints])
        #self.RT_carriers_dE_V_xc        = np.zeros([self.ntimes,self.nbands,self.nkpoints])
        self.times = [ time for time,filename in ordered_files]

        for n,(time,filename) in enumerate(ordered_files):

            #open database for each k-point
            db = Dataset("%s/%s"%(self.path,filename))

            self.RT_carriers_delta_f[n]          = db['RT_carriers_delta_f'][:].reshape([self.nkpoints,self.nbands])

            #self.RT_carriers_dE_Self_Energy[n]   = db['RT_carriers_dE_Self_Energy'][:].reshape([self.nkpoints,self.nbands])
            #self.RT_carriers_dE_V_xc[n]          = db['RT_carriers_dE_V_xc'][:].reshape([self.nbands,self.nkpoints])

            #close database
            db.close()

    def integrate(self):
        self.occupations = np.zeros([self.ntimes,self.nkpoints,self.nbands])

        for t in xrange(0,self.ntimes):

            #"delta_f" is df(t)-df(t0), so total occupation
            self.occupations[t] = self.RT_carriers_delta_f[t]

    def get_path(self,path,kpts=None):
        """ Obtain a list of indexes and kpoints that belong to the regular mesh
        """
        if kpts is None:
            kpts, nks, nss = self.expand_kpts()
        else:
            nks = range(len(kpts))

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
                #print ("%12.8lf "*3)%tuple(kpt), index

        self.bands_kpoints = bands_kpoints
        self.bands_indexes = bands_indexes
        self.bands_highsym_qpts = path_car

        print 'Path generated using %d kpoints.'%len(bands_kpoints)

        # Calculate distances
        bands_distances = [0]
        distance = 0
        for nk in range(1,len(bands_kpoints)):
            distance += np.linalg.norm(bands_kpoints[nk]-bands_kpoints[nk-1])
            bands_distances.append(distance)

        self.bands_distances = bands_distances

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
        for nk,k in enumerate(self.kpts_car):
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

        print "%d kpoints expanded to %d"%(len(self.kpts_car),len(kpoints_full))

        return self.kpoints_full, self.kpoints_indexes, self.symmetry_indexes

    def __str__(self):
        s = ""
        s += "nkpoints: %d\n"%self.nkpoints
        s += "min_band: %d\n"%self.nband_min
        s += "max_band: %d\n"%self.nband_max
        return s

