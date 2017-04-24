from __future__ import print_function
from __future__ import division
# Copyright (c) 2017, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
from builtins import zip
from builtins import range
from builtins import object
from past.utils import old_div
from yambopy import *
from netCDF4 import Dataset

class YamboStaticScreeningDB(object):
    """
    Class to handle static screening databases from Yambo
    
    This reads the databases ``ndb.em1s*``
    There :math:`v\chi(\omega=0)` is stored.
    
    To calculate epsilon (static dielectric function) we do:

    .. math::

        \epsilon^{-1} = 1-v\chi

    """
    def __init__(self,save='.',filename='ndb.em1s',db1='ns.db1'):
        self.save = save
        self.filename = filename

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
            database = Dataset("%s/%s"%(self.save,self.filename), 'r')
        except:
            raise IOError("Error opening %s/%s in YamboStaticScreeningDB"%(self.save,self.filename))

        #read some parameters
        size,nbands,eh = database['X_PARS_1'][:3]
        self.size = int(size)
        self.nbands = int(nbands)
        self.eh = eh

        #read gvectors
        gvectors = np.rint(database['X_RL_vecs'][:].T)
        self.gvectors = np.array([old_div(g,self.alat)  for g in gvectors])
        self.ngvectors = len(self.gvectors)
        
        #read q-points
        qpoints = database['HEAD_QPT'][:].T
        self.qpoints = np.array([old_div(q,self.alat)  for q in qpoints])
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
            filename = "%s/%s_fragment_%d"%(self.save,self.filename,nq+1)
            try:
                db = Dataset(filename)
            except:
                print("warning: failed to read %s"%filename)


            #static screening means we have only one frequency
            # this try except is because the way this is sotored has changed in yambo
            try:
                re, im = db['X_Q_%d'%(nq+1)][0,:]
            except:
                re, im = db['X_Q_%d'%(nq+1)][0,:].T

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
            db = Dataset("%s/%s"%(path,fname),'r+')
            db['X_Q_%d'%(nq+1)][0,0,:] = X[nq].real
            db['X_Q_%d'%(nq+1)][0,1,:] = X[nq].imag
            db.close()

    def writetxt(self,filename='em1s.dat',ng1=0,ng2=0,volume=False):
        """
        Write vVepsilon_{g1=0,g2=0} (q) as a funciton of |q| on a text file
        volume -> multiply by the volume
        """
        x,y = self._geteq(ng1=ng1,ng2=ng2,volume=volume)
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
 
        In the database we find vX(\omega=0) where:
        v -> coulomb interaction (truncated or not)
        X -> electronic response function

        Arguments:
            ng1, ng2 -> Choose local field components
            volume   -> Normalize with the volume of the cell
        """
        x = [np.linalg.norm(q) for q in self.qpoints]
        y = [np.linalg.inv(1+xq)[0,0] for xq in self.X ]
      
        #order according to the distance
        x, y = list(zip(*sorted(zip(x, y))))
        y = np.array(y)

        #scale by volume?
        if volume: y *= self.volume 

        return x,y

    def _getvxq(self,ng1=0,ng2=0,volume=False): 
        """
        Get vX_{ng1,ng2} a function of |q|
        vX is a matrix with size equal to the number of local fields components
 
        In the database we find vX(\omega=0) where:
        v -> coulomb interaction (truncated or not)
        X -> electronic response function

        Arguments:
            ng1, ng2 -> Choose local field components
            volume   -> Normalize with the volume of the cell
        """
        x = [np.linalg.norm(q) for q in self.qpoints]
        y = [xq[ng2,ng1] for xq in self.X ]
      
        #order according to the distance
        x, y = list(zip(*sorted(zip(x, y))))
        y = np.array(y)

        #scale by volume?
        if volume: y *= self.volume 

        return x,y
    
    def plot(self,ax,volume=False,**kwargs):
        """
        Plot the static screening as a function of |q|
        
        Arguments
        ax   -> Instance of the matplotlib axes or some other object with the plot method
        func -> Function to apply to the dielectric function
        """

        #get vX
        x,vX = self._getvxq(volume=volume)
    
        #when plotting we apply a funciton to epsilon to represent it, by default the |x|
        ax.plot(x,(1+vX).real,**kwargs)
        ax.set_xlabel('$|q|$')
        ax.set_ylabel('$\epsilon^{-1}_{00}(\omega=0)$')

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
