# Copyright (c) 2018, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
from yambopy import *
from netCDF4 import Dataset
from yambopy.lattice import rec_lat, car_red

class YamboStaticScreeningDB(object):
    """
    Class to handle static screening databases from Yambo
    
    This reads the databases ``ndb.em1s*``
    There :math:`√v(q,g1) \chi_{g1,g2} (q,\omega=0) √v(q,g2)` is stored.
    
    To calculate epsilon (static dielectric function) we do:

    .. math::

        \epsilon^{-1}_{g1,g2}(q) = 1-v(q,g1)\chi_{g1,g2}
        
    The symmetric and asymmetric formulations coincide for the head g1=g2=0
    """
    def __init__(self,save='.',em1s='.',filename='ndb.em1s',db1='ns.db1'):
        self.save = save
        self.em1s = em1s
        self.filename = filename

        #read lattice parameters
        if os.path.isfile('%s/%s'%(self.save,db1)):
            try:
                database = Dataset("%s/%s"%(self.save,db1), 'r')
                self.alat = database.variables['LATTICE_PARAMETER'][:]
                self.lat  = database.variables['LATTICE_VECTORS'][:].T
                self.volume = np.linalg.det(self.lat)
                self.rlat = rec_lat(self.lat)
            except:
                raise IOError("Error opening %s."%db1)
        else:
            raise FileNotFoundError("File %s not found."%db1)

        #read em1s database
        if os.path.isfile("%s/%s"%(self.em1s,self.filename)): 
            try:
                database = Dataset("%s/%s"%(self.em1s,self.filename), 'r')
            except:
                raise IOError("Error opening %s/%s in YamboStaticScreeningDB"%(self.save,self.filename))
        else:
            raise FileNotFoundError("File %s not found."%self.filename)

        #read some parameters
        size,nbands,eh = database.variables['X_PARS_1'][:3]
        self.size = int(size)
        self.nbands = int(nbands)
        self.eh = eh

        #read gvectors
        gvectors = np.rint(database.variables['X_RL_vecs'][:].T)
        self.gvectors = np.array([g/self.alat  for g in gvectors])
        self.ngvectors = len(self.gvectors)
        
        #read q-points
        self.iku_qpoints = database.variables['HEAD_QPT'][:].T
        self.car_qpoints = np.array([ q/self.alat for q in self.iku_qpoints ])
        self.red_qpoints = car_red(self.car_qpoints,self.rlat) 
        self.nqpoints = len(self.car_qpoints)

        try:
            database.variables['CUTOFF'][:]
            self.cutoff = str(database.variables['CUTOFF'][:][0],'UTF-8').strip()
        except: IndexError
        
        #read fragments
        read_fragments=True
        for iQ in range(self.nqpoints):
            if not os.path.isfile("%s/%s_fragment_%d"%(self.em1s,self.filename,iQ+1)): read_fragments=False
        if read_fragments: self.readDBs()

    def readDBs(self):
        """
        Read the yambo databases
        """

        #create database to hold all the X data
        self.X = np.zeros([self.nqpoints,self.size,self.size],dtype=np.complex64)
        for nq in range(self.nqpoints):

            #open database for each k-point
            filename = "%s/%s_fragment_%d"%(self.save,self.filename,nq+1)
            try:
                database = Dataset(filename)
            except:
                print("warning: failed to read %s"%filename)


            #static screening means we have only one frequency
            # this try except is because the way this is sotored has changed in yambo
            try:
                re, im = database.variables['X_Q_%d'%(nq+1)][0,:]
            except:
                re, im = database.variables['X_Q_%d'%(nq+1)][0,:].T

            self.X[nq] = re + 1j*im
         
            #close database
            database.close()

    def saveDBS(self,path):
        """
        Save the database
        """
        if os.path.isdir(path): shutil.rmtree(path)
        os.mkdir(path)

        #copy all the files
        oldpath = self.save
        filename = self.filename
        shutil.copyfile("%s/%s"%(oldpath,filename),"%s/%s"%(path,filename))
        for nq in range(self.nqpoints):
            fname = "%s_fragment_%d"%(filename,nq+1)
            shutil.copyfile("%s/%s"%(oldpath,fname),"%s/%s"%(path,fname))

        #edit with the new wfs
        X = self.X
        for nq in range(self.nqpoints):
            fname = "%s_fragment_%d"%(filename,nq+1)
            database = Dataset("%s/%s"%(path,fname),'r+')
            database.variables['X_Q_%d'%(nq+1)][0,0,:] = X[nq].real
            database.variables['X_Q_%d'%(nq+1)][0,1,:] = X[nq].imag
            database.close()

    def writetxt(self,filename='em1s.dat',ng1=0,ng2=0,volume=False):
        """
        Write vVepsilon_{g1=0,g2=0} (q) as a funciton of |q| on a text file
        volume -> multiply by the volume
        """
        x,y = self._geteq(volume=volume)
        np.savetxt(filename,np.array([x,y]).T)
    
    def get_g_index(self,g):
        """
        get the index of the gvectors.
        If the gvector is not present return None
        """
        for ng,gvec in enumerate(self.gvectors):
            if np.isclose(g,gvec).all():
                return ng
        return None

    def _geteq(self,volume=False): 
        """
        Get epsilon_{0,0} = [1/(1+vX)]_{0,0} a function of |q|
        vX is a matrix with size equal to the number of local fields components
 
        In the database we find √vX√v(\omega=0) where:
        v -> coulomb interaction (truncated or not)
        X -> electronic response function

        Arguments:
            ng1, ng2 -> Choose local field components
            volume   -> Normalize with the volume of the cell
        """
        x = [np.linalg.norm(q) for q in self.car_qpoints]
        y = [np.linalg.inv(np.eye(self.ngvectors)+xq)[0,0] for xq in self.X ]
      
        #order according to the distance
        x, y = list(zip(*sorted(zip(x, y))))
        y = np.array(y)

        #scale by volume?
        if volume: y *= self.volume 

        return x,y

    def _getvxq(self,ng1=0,ng2=0,volume=False,symm=True): 
        """
        Get vX_{ng1,ng2} a function of |q|
        vX is a matrix with size equal to the number of local fields components
 
        In the database we find √vX√v(\omega=0) where:
        v -> coulomb interaction (truncated or not)
        X -> electronic response function

        Arguments:
            ng1, ng2 -> Choose local field components
            volume   -> Normalize with the volume of the cell
            symm     -> True:  √v(q,g1) X_{g1,g2}(q) √v(q,g2)
                        False: v(q,g1) X_{g1,g2}(q) TO BE IMPLEMENTED
        """
        x = [np.linalg.norm(q) for q in self.car_qpoints]
        if symm:
            y = [xq[ng2,ng1] for xq in self.X ]
        else: 
            raise NotImplementedError("vXq with symm=False is not presently implemented.")
      
        #order according to the distance
        x, y = list(zip(*sorted(zip(x, y))))
        y = np.array(y)

        #scale by volume?
        if volume: y *= self.volume 

        return x,y
    
    def plot(self,ax,ng1=0,ng2=0,volume=False,symm=True,**kwargs):
        """
        Plot the static screening as a function of |q|
        
        Arguments
        ax   -> Instance of the matplotlib axes or some other object with the plot method
        func -> Function to apply to the dielectric function
        """

        #get vX_{00}
        x,vX = self._getvxq(ng1=ng1,ng2=ng2,volume=volume,symm=symm)
    
        #when plotting we apply a funciton to epsilon to represent it, by default the |x|
        ax.plot(x,(1+vX).real,**kwargs)
        ax.set_xlabel('$|q|$')
        ax.set_ylabel('$\epsilon^{-1}_{%d%d}(\omega=0)$'%(ng1,ng2))

    def __str__(self):
        s = ""
        s += "nqpoints: %d\n"%self.nqpoints
        s += "X size:   %d\n"%self.size
        s += "cutoff: %s\n"%self.cutoff
        return s


if __name__ == "__main__":

    ys = YamboStaticScreeningDB()
    print(ys)
  
    #plot static screening 
    ax = plt.gca()
    ys.plot(ax)
    plt.show()
