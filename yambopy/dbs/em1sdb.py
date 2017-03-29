# Copyright (c) 2016, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
from yambopy import *
from yambopy.netcdf import *

class YamboStaticScreeningDB():
    """
    Class to handle static screening databases from Yambo
    """
    def __init__(self,save='.',filename='ndb.em1s',db1='ns.db1'):
        self.save = save
        self.filename = "%s/%s"%(save,filename)

        #read the lattice paramaters
        try:
            #posibilities where to find db1
            for filename in ['%s/%s'%(save,db1),'%s/../SAVE/%s'%(save,db1)]:
                if os.path.isfile(filename):
                    break
            database = Dataset(filename, 'r')
            self.alat = database['LATTICE_PARAMETER'][:]
            self.lat  = database['LATTICE_VECTORS'][:].T
            self.volume = np.linalg.det(self.lat)
        except:
            raise IOError("Error opening %s in YamboStaticScreeningDB"%filename)

        #read em1s database
        try:
            database = Dataset(self.filename, 'r')
        except:
            raise IOError("Error opening %s in YamboStaticScreeningDB"%self.filename)

        #read some parameters
        size,nbands,eh = database['X_PARS_1'][:3]
        self.size = int(size)
        self.nbands = int(nbands)
        self.eh = eh

        #read gvectors
        gvectors = np.rint(database['X_RL_vecs'][:].T)
        self.gvectors = np.array([g/self.alat  for g in gvectors])
        self.ngvectors = len(self.gvectors)
        
        #read q-points
        qpoints = database['HEAD_QPT'][:].T
        self.qpoints = np.array([q/self.alat  for q in qpoints])
        self.nqpoints = len(self.qpoints)
        
        #are we usign coulomb cutoff?
        self.cutoff = "".join(database['CUTOFF'][:][0]).strip()
        
        self.readDBs()

    def readDBs(self):
        """
        Read the yambo databases
        """

        #create database to hold all the X data
        self.X = np.zeros([self.nqpoints,self.size,self.size],dtype=np.complex64)
        for nq in range(self.nqpoints):

            #open database for each k-point
            filename = "%s_fragment_%d"%(self.filename,nq+1)
            db = Dataset(filename)

            #static screening means we have only one frequency
            re, im = db['X_Q_%d'%(nq+1)][:][0]
            self.X[nq] = re + 1j*im
 
            #close database
            db.close()

    def saveDBS(self,path):
        """
        Save the database
        """
        if os.path.isdir(path): shutil.rmtree(path)
        os.mkdir(path)

        #copy all the files
        oldpath = self.path
        filename = self.filename
        shutil.copyfile("%s/%s"%(oldpath,filename),"%s/%s"%(path,filename))
        for nk in xrange(self.nkpoints):
            fname = "%s_fragments_%d_1"%(filename,nk+1)
            shutil.copyfile("%s/%s"%(oldpath,fname),"%s/%s"%(path,fname))

        #edit with the new wfs
        wf = self.wf
        for nk in xrange(self.nkpoints):
            fname = "%s_fragments_%d_1"%(filename,nk+1)
            db = Dataset("%s/%s"%(path,fname),'r+')
            db['WF_REAL_COMPONENTS_@_K%d_BAND_GRP_1'%(nk+1)][:] = wf[nk].real
            db['WF_IM_COMPONENTS_@_K%d_BAND_GRP_1'%(nk+1)][:] = wf[nk].imag
            db.close()

    def get_g_index(self,g):
        """
        get the index of the gvectors.
        If the gvector is not present return None
        """
        for ng,gvec in enumerate(self.gvectors):
            if np.isclose(g,gvec).all():
                return ng
        return None
        
    def plot(self,ax,ng1=0,ng2=0,**kwargs):
        """
        Plot the static screening
        
        Arguments
        ax -> Instance of the matplotlib axes or some other object with the plot method
        """
        M1 = np.eye(self.size)
        x = [np.linalg.norm(q) for q in self.qpoints]
        y = [np.abs(xq)[ng2,ng1] for xq in self.X ]
      
        #order according to the distance
        x, y = zip(*sorted(zip(x, y)))        
        y = np.array(y)*self.volume
 
        ax.plot(x,y.real,**kwargs)
        ax.set_xlabel('|q|')

    def __str__(self):
        s = ""
        s += "nqpoints: %d\n"%self.nqpoints
        s += "X size:   %d\n"%self.size
        s += "cutoff: %s\n"%self.cutoff
        return s


if __name__ == "__main__":

    ys = YamboStaticScreeningDB()
    print ys
  
    #plot static screening 
    ax = plt.gca()
    ys.plot(ax)
    plt.show()
