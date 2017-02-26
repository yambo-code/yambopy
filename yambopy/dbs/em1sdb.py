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
    def __init__(self,save='.',filename='ndb.em1s'):
        self.filename = "%s/%s"%(save,filename)

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

        #read number of q-points
        self.qpoints = database['HEAD_QPT'][:].T
        self.nqpoints = len(self.qpoints)
        
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

    def plot(self,ax):
        """
        Plot the static screening
        
        Arguments
        ax -> Instance of the matplotlib axes or some other object with the plot method
        """
        M1 = np.eye(self.size)
        x = [np.linalg.norm(q) for q in self.qpoints]
        y = [np.linalg.inv(M1-xq)[0,0] for xq in self.X ]
      
        #order according to the distance
        x, y = zip(*sorted(zip(x, y)))        
        y = np.array(y)
 
        ax.plot(x,y.real)
        ax.set_xlabel('|q|')

    def __str__(self):
        s = ""
        s += "nqpoints: %d"%self.nqpoints
        return s


if __name__ == "__main__":

    ys = YamboStaticScreeningDB()
    print ys
  
    #plot static screening 
    ax = plt.gca()
    ys.plot(ax)
    plt.show()
