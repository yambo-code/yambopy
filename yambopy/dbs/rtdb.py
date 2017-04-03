# Copyright (c) 2016, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
from yambopy import *
from yambopy.netcdf import *
from yambopy.plot import *

class YamboRTDB():
    """
    Open the RT databases and store it in a RTDB class
    """
    def __init__(self,electrons,path='.',referencedb='ndb.RT_reference_components'):
        self.path = path
        self.referencedb = referencedb

        #read reference database
        db = Dataset("%s/%s"%(path,referencedb))
        self.nband_min, self.nband_max, self.nkpoints = db['RT_vars'][:].astype(int)
        self.nbands = self.nband_max - self.nband_min + 1

        #get energies of bands
        self.electrons = electrons
        self.eigenvalues = self.electrons.eigenvalues[:,self.nband_min:self.nband_max+1]

        #read the databases
        self.readDB()

        #integrate the occupations
        self.integrate()

    def readDB(self):
        """
        """

        #get how many rt databases exist
        files = [ filename for filename in  os.listdir(self.path) if 'ndb.RT_carriers_Time' in filename]
        print "number fo files:", len(files)

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
        self.RT_carriers_delta_f        = np.zeros([self.ntimes,self.nbands,self.nkpoints])
        self.RT_carriers_delta_f        = np.zeros([self.ntimes,self.nkpoints,self.nbands])
        #self.RT_carriers_dE_Self_Energy = np.zeros([self.ntimes,self.nbands,self.nkpoints])
        #self.RT_carriers_dE_V_xc        = np.zeros([self.ntimes,self.nbands,self.nkpoints])
        self.times = [ time for time,filename in ordered_files]

        for n,(time,filename) in enumerate(ordered_files):

            #open database for each k-point
            db = Dataset("%s/%s"%(self.path,filename))

            self.RT_carriers_delta_f[n]          = db['RT_carriers_delta_f'][:].reshape([self.nkpoints,self.nbands])
            #self.RT_carriers_delta_f[n]          = db['RT_carriers_delta_f'][:].reshape([self.nbands,self.nkpoints])

            #self.RT_carriers_dE_Self_Energy[n]   = db['RT_carriers_dE_Self_Energy'][:].reshape([self.nkpoints,self.nbands])
            #self.RT_carriers_dE_V_xc[n]          = db['RT_carriers_dE_V_xc'][:].reshape([self.nbands,self.nkpoints])

            #close database
            db.close()

    def integrate(self):
        occupations = np.zeros([self.nkpoints,self.nbands])
        self.occupations_df = np.zeros([self.ntimes,self.nkpoints,self.nbands])
        self.occupations = np.zeros([self.ntimes,self.nkpoints,self.nbands])

        time = self.times
        for t in xrange(1,self.ntimes):
            dt = time[t]-time[t-1]

            #change at time t
            self.occupations_df[t] = self.RT_carriers_delta_f[t]*dt

            #linear integration
            occupations += self.occupations_df[t]

            #store current occupation
            self.occupations[t] = occupations

    def __str__(self):
        s = ""
        s += "nkpoints: %d\n"%self.nkpoints
        s += "min_band: %d\n"%self.nband_min
        s += "max_band: %d\n"%self.nband_max
        return s

