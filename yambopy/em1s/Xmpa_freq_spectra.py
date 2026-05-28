
# License-Identifier: GPL
#
# Copyright (C) 2024 The Yambo Team
#
# Authors: GS
#
# This file is part of the yambopy project
#
"""
This file contains a class to reconstruct and visualize X from the MPA poles
The MPA poles are contained in the mpa_ER database

Here, the poles and residues are summed.
Note that the MPA procedure has changed from yambo 5.4

Before yambo 5.4, the residues reffered to XV.
From yambo 5.4 onward the residues refer to W. As such, also the RIM_qpg or bare_qpg needs to be accessed
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import matplotlib
import os
from yambopy.lattice import vol_lat, rec_lat, car_red, red_car
from yambopy.dbs.latticedb import YamboLatticeDB
from yambopy.units import ha2ev

class XmpaDB(object):
    """
    Class to handle screening databases from Yambo
    
    This reads the databases ``ndb.mpa_ER*`` 
    
    If a Coulomb truncation is used for :math:`v(q,g)`, then the database ``ndb.RIM`` or ``ndb.cutoff`` are also read.

    .. math::

        vX{g1,g2}(q) = sum_n R_n*(1/(w-Omega_n+i delta)-1/(w+Omega_n-i delta))
        
    """
    def __init__(self,save='.',mpa_folder='.',db1='ns.db1',nqs=None,coul_fl=None,do_not_read_cutoff=False):

        self.save = save
        self.mpa = mpa_folder
        self.coul= mpa_folder if coul_fl==None else coul_fl
        self.no_cutoff = do_not_read_cutoff
        self.filename='ndb.mpa_ER'
          
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
        if not os.path.isfile("%s/%s"%(self.mpa,self.filename)): 
           raise FileNotFoundError("File %s not found."%self.filename)

        try:
            database = Dataset("%s/%s"%(self.mpa,self.filename), 'r')
        except:
            raise IOError("Error opening %s/%s in MPA fragments DB"%(self.save,self.filename))
        
        #read some parameters
        size,npoles = database.variables['MPA_PARS_1'][:2]
        self.size = int(size)
        self.npoles = int(npoles)
        head_version=database.variables["HEAD_VERSION"][:]
        self.head_version  = head_version[0]+0.1*head_version[1]+0.01*head_version[2]

        #read gvectors used
        gvectors          = np.array(gvectors_full[:self.size])
        self.gvectors     = np.array([g/self.alat for g in gvectors])
        self.red_gvectors = car_red(self.gvectors,self.rlat)
        self.ngvectors    = len(self.gvectors)
        
        #read q-points
        self.iku_qpoints = database.variables['HEAD_QPT'][:].T
        self.car_qpoints = np.array([ q/self.alat for q in self.iku_qpoints ]) #atomic units
        self.red_qpoints = car_red(self.car_qpoints,self.rlat) 
        self.nqpoints = len(self.car_qpoints)

        # The d3kfactor used here is equal to
        #  
        # d3kfactor= nspin(=2)/(2*pi)**3*RL_vol/Nk
        #
        # It corresponds to the q_weight factor used in yambo 
        # to calculate the gamp(G,G') 
        # 
        ylat = YamboLatticeDB.from_db_file(filename='%s/%s'%(self.save,db1),Expand=True)
        self.d3kfactor = 1/(2*np.pi)*vol_lat(self.rlat)/ylat.nkpoints
        # Multiply weights by Nk to avoid double-counting
        self.weights=ylat.weights_ibz[:self.nqpoints]*ylat.nkpoints 

        # check if rim database present
        self.rim=True
        if not os.path.isfile("%s/ndb.RIM"%(self.coul)):
           self.rim=False 

        # for yambo 5.4 onward cutoff is stored in mpa_ER, prior to it read from Xmpa
        if (self.head_version>= 5.4):
         try:
            database.variables['CUTOFF'][:]
            self.cutoff = str(database.variables['CUTOFF'][:][0],'UTF-8').strip()
         except: IndexError
        else:
         if not os.path.isfile("%s/ndb.Xmpa"%(self.mpa)):
            raise FileNotFoundError("File ndb.Xmpa not found, needed to check if a cutoff is present")
         try:
            database = Dataset("%s/ndb.Xmpa"%(self.mpa), 'r')
            database.variables['CUTOFF'][:]
            self.cutoff = str(database.variables['CUTOFF'][:][0],'UTF-8').strip()
         except: raise IOError("Error in accessing ndb.Xmpa database")

        #set the number of fragments to read
        if isinstance(nqs,int):
           self.range_nqs=[nqs]
        elif isinstance(nqs,list) or isinstance(nqs,range):
           self.range_nqs=nqs
        else:
           self.range_nqs=range(self.nqpoints)                

        #read fragments
        read_fragments=True
        for iQ in range(self.nqpoints):
            if not os.path.isfile("%s/%s_fragment_%d"%(self.mpa,self.filename,iQ+1)): read_fragments=False
        if read_fragments: self.readDBs(self.range_nqs) # get sqrt(v)*X*sqrt(v)

        self.readDBs(self.range_nqs)
        #get qpg
        self.get_Coulomb()

    def readDBs(self,range_nqs):
        """
        Read the yambo databases, we access the fragment 1 to set the number of poles
        """

        #create database to hold all the poles and residues
        self.poles = np.zeros([self.nqpoints,self.npoles,self.size,self.size],dtype=np.complex64)
        self.residues = np.zeros([self.nqpoints,self.npoles,self.size,self.size],dtype=np.complex64)

        for nq in range_nqs:

            #open database for each k-point
            filename = "%s/%s_fragment_%d"%(self.mpa,self.filename,nq+1)
            try:
                database = Dataset(filename)
            except:
                print("warning: failed to read %s"%filename)

            # poles
            re = database.variables['MPA_E_Q_%d'%(nq+1)][...,0]
            im = database.variables['MPA_E_Q_%d'%(nq+1)][...,1]

            self.poles[nq,:] = re[:] + 1j*im[:]

            # residues
            re = database.variables['MPA_R_Q_%d'%(nq+1)][...,0]
            im = database.variables['MPA_R_Q_%d'%(nq+1)][...,1]

            self.residues[nq,:] = re[:] + 1j*im[:]
         
            #close database
            database.close()

    def get_Coulomb(self):
        """
        If rim_qpg database is present, read it

        Otherwise, as second option, look for ndb.cutoff and parse it.
        Otherwise, construct bare 3D potential.

        Returns sqrt_V[Nq,Ng]
        """
        q_p_G=np.zeros((self.nqpoints,self.size),dtype=np.float32)
        self.sqrt_V=np.zeros((self.nqpoints,self.size),dtype=np.float32)

        alloc_dim=0
        if self.rim==True:
           database = Dataset("%s/ndb.RIM"%self.coul, 'r') 
           alloc_dim=database.variables["RIM_qpg"].shape[0]
           q_p_G[:,:alloc_dim] = np.array(np.diagonal(database.variables["RIM_qpg"],axis1=0, axis2=1))
           self.sqrt_V[:,:alloc_dim] = np.sqrt(q_p_G[:,:alloc_dim])
           
        if self.cutoff!='none' and not self.no_cutoff:

            if os.path.isfile('%s/ndb.cutoff'%self.mpa):
              try:
                 database = Dataset("%s/ndb.cutoff"%self.mpa, 'r')
                 q_p_G[:,alloc_dim:self.size] = np.array(database.variables["CUT_BARE_QPG"][alloc_dim:self.size,:,0]).T
                 self.sqrt_V[:,alloc_dim:self.size] = np.sqrt(self.d3kfactor*4.0*np.pi)/q_p_G[:,alloc_dim:self.size]
                 database.close()
              except:
                 raise IOError("Error opening ndb.cutoff.")
            else:
                print("[WARNING] Cutoff %s was used but ndb.cutoff not found in %s. Make sure this is fine for what you want!"%(self.cutoff,self.mpa))

        else:

            nrm = np.linalg.norm
            for iq in range(self.nqpoints):
                for ig in range(alloc_dim,self.size):
                        Q = 2.*np.pi*self.car_qpoints[iq]
                        G = 2.*np.pi*self.gvectors[ig]
                        q_p_G[iq,ig] = nrm(Q+G)
                        if (q_p_G[iq,ig]==0.0): q_p_G[iq,ig]=1.e-8

            self.sqrt_V[:,alloc_dim:self.size] = np.sqrt(self.d3kfactor*4*np.pi)/q_p_G[:,alloc_dim:self.size]

    def get_X(self,ig1=0,ig2=0,w_rnge=None,nw=1000,damp=5e-3):

       # set frequency range of VX
       self.nw =nw
       if w_rnge==None:
        self.w_rnge = np.linspace(0,15,self.nw)
       else:
        self.w_rnge = np.linspace(w_rnge[0],w_rnge[-1],self.nw)
       w_i=self.w_rnge[:,np.newaxis,np.newaxis]/ha2ev
       X_MPA=np.zeros((nw,self.nqpoints),dtype=np.complex64)
     
       if (self.head_version>= 5.4):
          X_MPA[:,self.range_nqs] = np.sum(self.residues[self.range_nqs,:,ig1,ig2]/self.sqrt_V[self.range_nqs,np.newaxis,ig2]/self.sqrt_V[self.range_nqs,np.newaxis,ig1]* \
                                    (1./(w_i-self.poles[self.range_nqs,:,ig1,ig2]+damp)-1./(w_i+self.poles[self.range_nqs,:,ig1,ig2]-damp)),axis=2) 
       else:
          X_MPA[:,self.range_nqs] = np.sum(self.residues[self.range_nqs,:,ig1,ig2]* \
                                    (1./(w_i-self.poles[self.range_nqs,:,ig1,ig2]+damp)-1./(w_i+self.poles[self.range_nqs,:,ig1,ig2]-damp)),axis=2) 

       # return VX(g1,g2)
       return X_MPA
             
    def get_W(self,ig1=0,ig2=0,w_rnge=None,nw=1000,damp=5e-3):

       # set frequency range of W
       self.nw =nw
       if w_rnge==None:
        self.w_rnge = np.linspace(0,15,self.nw)
       else:
        self.w_rnge = np.linspace(w_rnge[0],w_rnge[-1],self.nw)
       w_i=self.w_rnge[:,np.newaxis,np.newaxis]/ha2ev
       W_MPA=np.zeros((nw,self.nqpoints),dtype=np.complex64)

       if (self.head_version>= 5.4):
          W_MPA[:,self.range_nqs] = np.sum(self.residues[self.range_nqs,:,ig1,ig2]* \
                                    (1./(w_i-self.poles[self.range_nqs,:,ig1,ig2]+damp)-1./(w_i+self.poles[self.range_nqs,:,ig1,ig2]-damp)),axis=2) 
       else:
          W_MPA[:,self.range_nqs] = np.sum(self.residues[self.range_nqs,:,ig1,ig2]*self.sqrt_V[iq,np.newaxis,ig2]*self.sqrt_V[iq,np.newaxis,ig1]* \
                                    (1./(w_i-self.poles[self.range_nqs,:,ig1,ig2]+damp)-1./(w_i+self.poles[self.range_nqs,:,ig1,ig2]-damp)),axis=2) 

       #return W(g1,g2) in Ha
       return W_MPA


    def get_W_summed(self,ig1=0,ig2=0,w_rnge=None,nw=1000,damp=5e-3):
       
       # we build WGG' summed on all the q-points of the BZ, the iBZ are multiplied by the proper weights
       # in a 2D WGG' depends on the supercell size, unless also summing on the RL vectors along the truncated direction
       if (len(self.range_nqs)<self.nqpoints):

          raise ValueError("Error: not enough MPA_ER* fragments have been read, all the q-points of the IBZ are needed")

       # set frequency range of W and get weights
       self.nw =nw
       if w_rnge==None:
        self.w_rnge = np.linspace(0,15,self.nw)
       else:
        self.w_rnge = np.linspace(w_rnge[0],w_rnge[-1],self.nw)
       w_i=self.w_rnge[:,np.newaxis,np.newaxis]/ha2ev
       weights=self.weights[:,np.newaxis]
       W_MPA=np.zeros(nw,dtype=np.complex64)

       if (self.head_version>= 5.4):
          W_MPA = np.sum(weights*self.residues[...,ig1,ig2]*(1./(w_i-self.poles[...,ig1,ig2]+damp)-1./(w_i+self.poles[...,ig1,ig2]-damp)),axis=(1,2))

       else:
          W_MPA = np.sum(weights*self.residues[...,ig1,ig2]*(1./(w_i-self.poles[...,ig1,ig2]+damp)-1./(w_i+self.poles[...,ig1,ig2]-damp)) \
                         *self.sqrt_V[...,ig2]*self.sqrt_V[...,ig1],axis=(1,2))

       # return W(g1,g2) q-summed in Ha
       return W_MPA
