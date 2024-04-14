#
# License-Identifier: GPL
#
# Copyright (C) 2024 The Yambo Team
#
# Authors: HPC, FP
#
# This file is part of the yambopy project
#
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from yambopy.lattice import rec_lat, car_red
from yambopy.tools.string import marquee

class YamboStaticScreeningDB(object):
    """
    Class to handle static screening databases from Yambo
    
    This reads the databases ``ndb.em1s*`` or equivalently the static part of the ``ndb.pp*``.
    There :math:`√v(q,g1) \chi_{g1,g2} (q,\omega=0) √v(q,g2)` is stored.
    
    If a Coulomb truncation is used for :math:`v(q,g)`, then the database ``ndb.cutoff`` is also read.

    We obtain the microscopic inverse dielectric function and the macroscopic dielectric function as:

    .. math::

        \epsilon^{-1}_{g1,g2}(q) = 1+v(q,g1)\chi_{g1,g2}
        \epsilon_{0,0}(q)={\epsilon^{-1}(q)}^{-1}_{0,0}
        
    """
    def __init__(self,save='.',em1s='.',filename='ndb.em1s',db1='ns.db1',do_not_read_cutoff=False):
        self.save = save
        self.em1s = em1s
        self.filename = filename
        self.no_cutoff = do_not_read_cutoff

        #read lattice parameters
        if os.path.isfile('%s/%s'%(self.save,db1)):
            try:
                database = Dataset("%s/%s"%(self.save,db1), 'r')
                self.alat = database.variables['LATTICE_PARAMETER'][:]
                self.lat  = database.variables['LATTICE_VECTORS'][:].T
                gvectors_full = database.variables['G-VECTORS'][:].T
                self.gvectors_full = np.array([ g/self.alat for g in gvectors_full ])
                self.volume = np.linalg.det(self.lat)
                self.rlat = rec_lat(self.lat)
            except:
                raise IOError("Error opening %s."%db1)
        else:
            raise FileNotFoundError("File %s not found."%db1)

        #read em1s database
        if not os.path.isfile("%s/%s"%(self.em1s,self.filename)): 
            if os.path.isfile("%s/%s"%(self.em1s,'ndb.pp')): self.filename = 'ndb.pp' # check for ppa instead of static screening
            else: raise FileNotFoundError("File %s not found."%self.filename)

        try:
            database = Dataset("%s/%s"%(self.em1s,self.filename), 'r')
        except:
            raise IOError("Error opening %s/%s in YamboStaticScreeningDB"%(self.save,self.filename))
        
        #read some parameters
        size,band_0,band_1 = database.variables['X_PARS_1'][:3]
        self.size = int(size)
        self.nbands = int(band_1-band_0+1)
        self.first_band, self.last_band = int(band_0), int(band_1)

        #read gvectors used for em1s
        gvectors          = np.array(database.variables['X_RL_vecs'][:].T)
        self.gvectors     = np.array([g/self.alat for g in gvectors])
        self.red_gvectors = car_red(self.gvectors,self.rlat)
        self.ngvectors    = len(self.gvectors)
        
        #read q-points
        self.iku_qpoints = database.variables['HEAD_QPT'][:].T
        self.car_qpoints = np.array([ q/self.alat for q in self.iku_qpoints ]) #atomic units
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
        if read_fragments: self.readDBs() # get sqrt(v)*X*sqrt(v)

        #get square root of Coulomb potential v(q,G) 
        self.get_Coulomb()

    def readDBs(self):
        """
        Read the yambo databases
        """

        #create database to hold all the X data
        self.X = np.zeros([self.nqpoints,self.size,self.size],dtype=np.complex64)
        for nq in range(self.nqpoints):

            #open database for each k-point
            filename = "%s/%s_fragment_%d"%(self.em1s,self.filename,nq+1)
            try:
                database = Dataset(filename)
            except:
                print("warning: failed to read %s"%filename)


            # static screening means we have only one frequency
            # this try except is because the way this is stored has changed in yambo
            # Reading like this is not the best way to do it, but it works for now 
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

    def writeeps(self,filename='em1s.dat',ng1=0,ng2=0,volume=False):
        """
        Write epsilon_{g1=0,g2=0} (q) as a function of |q| on a text file
        volume -> multiply by the volume
        """
        x,y = self._getepsq(volume=volume)
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

    def get_Coulomb(self):
        """
        By MZ

        If cutoff is present, look for ndb.cutoff and parse it.
        Otherwise, construct bare 3D potential.

        Returns sqrt_V[Nq,Ng]
        """  
  
        if self.cutoff!='none' and not self.no_cutoff:

            if os.path.isfile('%s/ndb.cutoff'%self.em1s):
                try:
                    database = Dataset("%s/ndb.cutoff"%self.em1s, 'r')
                    q_p_G_RE = np.array(database.variables["CUT_BARE_QPG"][:,:,0].T)
                    q_p_G_IM = np.array(database.variables["CUT_BARE_QPG"][:,:,1].T)
                    q_p_G = q_p_G_RE + 1j*q_p_G_IM
                    self.sqrt_V = np.sqrt(4.0*np.pi)/q_p_G
                    database.close()
                except:
                    raise IOError("Error opening ndb.cutoff.")
            else:
                print("[WARNING] Cutoff %s was used but ndb.cutoff not found in %s. Make sure this is fine for what you want!"%(self.cutoff,self.em1s))

        else:

            sqrt_V = np.zeros([self.nqpoints,self.ngvectors])
            nrm = np.linalg.norm
            for iq in range(self.nqpoints):
                for ig in range(self.ngvectors):
                        Q = 2.*np.pi*self.car_qpoints[iq]
                        G = 2.*np.pi*self.gvectors[ig]
                        QPG = nrm(Q+G)
                        if QPG==0.: QPG=1.e-5
                        sqrt_V[iq,ig] = np.sqrt(4.0*np.pi)/QPG        
            self.sqrt_V = sqrt_V

    def _getepsq(self,volume=False,use_trueX=False,indices=None): 
        """
        Get epsilon_{0,0} = [1/(1+vX)]_{0,0} as a function of |q|
        vX is a matrix with size equal to the number of local fields components
 
        In the database we find √vX√v(\omega=0) where:
        v -> coulomb interaction (truncated or not)
        X -> electronic response function

        Arguments:
            ng1, ng2  -> Choose local field components
            volume    -> Normalize with the volume of the cell
            use_trueX -> Use desymmetrised vX [testing]
            indices   -> Use a subset of the total q-points (e.g. along a certain direction)
        """
        if not use_trueX: 
            X = self.X
        if use_trueX:  
            _,_ = self.getem1s()
            X = self.trueX

        if indices is None: indices = [ iq for iq in range(self.nqpoints) ] # Use all q-points

        x = [np.linalg.norm(q) for q in self.car_qpoints[indices]]
        y = [np.linalg.inv(np.eye(self.ngvectors)+xq)[0,0] for xq in X ]
        y = np.array(y)[indices]

        #order according to the distance
        x, y = list(zip(*sorted(zip(x, y))))
        x = np.array(x)
        y = np.array(y)

        #scale by volume?
        if volume: y *= self.volume 

        return x,y
    
    def _getvq(self,ng1=0):
        """
        Get Coulomb potential v_ng1 as a function of |q|

        v -> coulomb interaction (truncated or not)

        The quantity obtained is : v(q,g1)

        Arguments:
            ng1 -> Choose local field component
        """
        x = [np.linalg.norm(q) for q in self.car_qpoints]
        y = [vq[ng1]**2. for vq in self.sqrt_V]

        #order according to the distance
        x, y = list(zip(*sorted(zip(x, y))))
        y = np.array(y)

        return x,y

    def _getvxq(self,ng1=0,ng2=0,volume=False): 
        """
        Get vX_{ng1,ng2} as a function of |q|
        vX is a matrix with size equal to the number of local fields components
 
        In the database we find √vX√v(\omega=0) where:
        v -> coulomb interaction (truncated or not)
        X -> electronic response function

        The quantity obtained is: √v(q,g1) X_{g1,g2}(q) √v(q,g2)

        Arguments:
            ng1, ng2 -> Choose local field components
            volume   -> Normalize with the volume of the cell
        """
        x = [np.linalg.norm(q) for q in self.car_qpoints]
        y = [xq[ng2,ng1] for xq in self.X ]
      
        #order according to the distance
        x, y = list(zip(*sorted(zip(x, y))))
        y = np.array(y)

        #scale by volume?
        if volume: y *= self.volume 

        return x,y
 
    def _getem1s(self,ng1=0,ng2=0,volume=False,indices=None):
        """
        Get eps^-1_{ng1,ng2} a function of |q|

        In the database we find √vX√v(\omega=0) where:
        v -> coulomb interaction (truncated or not)
        X -> electronic response function

        We need to explicitly use √v in order to obtain:

              eps^-1+{g1,g2} = 1+v(q,g1) X_{g1,g2}(q)
                             = 1 + √v_g1 √v_g1 X_g1g2 √v_g2/√v_g2

        This works for 
            - 3D bare √v
            - 2D cutoff √v positive definite (i.e., like slab z)

        Arguments:
            ng1, ng2 -> Choose local field components
            volume   -> Normalize with the volume of the cell
            indices   -> Use a subset of the total q-points (e.g. along a certain direction)
            symm     -> True:  √v(q,g1) X_{g1,g2}(q) √v(q,g2)
                        False: v(q,g1) X_{g1,g2}(q) TO BE IMPLEMENTED
        """
        # First of all store trueX (the desymmetrized form) as attribute
        trueX = np.zeros([self.nqpoints,self.size,self.size],dtype=np.complex64)

        for ig1 in range(self.ngvectors):
            for ig2 in range(self.ngvectors):
                trueX[:,ig1,ig2] = self.sqrt_V[:,ig1]*self.X[:,ig1,ig2]/self.sqrt_V[:,ig2]

        self.trueX = trueX 

        # Now compute quantities for plotting like in _getepsq (but without inversion and for the input g-indices)
        if indices is None: indices = [ iq for iq in range(self.nqpoints) ] # Use all q-points

        x = [np.linalg.norm(q) for q in self.car_qpoints[indices]]
        y = [xq[ng1,ng2] for xq in self.trueX ]
        y = np.array(y)[indices]

        #order according to the distance
        x, y = list(zip(*sorted(zip(x, y))))
        x = np.array(x)
        y = np.array(y)

        #scale by volume?
        if volume: y *= self.volume

        return x,y
   
    def plot_epsm1(self,ax,ng1=0,ng2=0,volume=False,symm=False,**kwargs):
        """
        Plot epsilon^-1_{ng1,ng2} as a function of |q|
        
        Arguments
        ax   -> Instance of the matplotlib axes or some other object with the plot method
        symm -> True:  plot symmetrized version 1 + √vX√v
        symm -> False: plot true em1s as 1+vX [Default]
        """

        #get √vX√v_{ng1,ng2}
        if symm==True:  x,vX = self._getvxq(ng1=ng1,ng2=ng2,volume=volume)
        #get vX_{ng1,ng2}
        if symm==False: x,vX = self._getem1s(ng1=ng1,ng2=ng2,volume=volume)   

        ax.plot(x,(1+vX).real,**kwargs)
        ax.set_xlabel('$|q|$')
        ax.set_ylabel('$\epsilon^{-1}_{%d%d}(\omega=0)$'%(ng1,ng2))

     
    def plot_eps(self,ax,ng1=0,ng2=0,volume=False,use_trueX=False,**kwargs):
        """
        Get epsilon_{0,0} = [1/(1+vX)]_{0,0} as a function of |q|
        """
        x,y = self._getepsq(volume=volume)
        ax.plot(x,y.real,**kwargs)
        ax.set_xlabel('$|q|$')
        ax.set_ylabel('$\epsilon_{%d%d}(\omega=0)$'%(ng1,ng2))

    def plot_v(self,ax,ng1=0,**kwargs):
        """
        Get v_{ng1} (truncated or not) as a function of |q|
        """
        x,y = self._getvq(ng1=ng1)
        ax.plot(x,y.real,**kwargs)
        ax.set_xlabel('$|q|$')
        ax.set_ylabel('$v_{%d}$'%ng1)

    def __str__(self):

        lines = []; app=lines.append
        app(marquee(self.__class__.__name__))

        app('filename:         %s'%self.filename)
        app('nqpoints (ibz):   %d'%self.nqpoints)
        app('X size (G-space): %d'%self.size) 
        app('cutoff:           %s'%self.cutoff) 
        app('bands:            %d (%d-%d)'%(self.nbands,self.first_band,self.last_band))

        return "\n".join(lines)


if __name__ == "__main__":

    ys = YamboStaticScreeningDB()
    print(ys)
  
    #plot static screening 
    ax = plt.gca()
    ys.plot_epsm1(ax)
    plt.show()
