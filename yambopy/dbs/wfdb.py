# Copyright (c) 2018, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
from yambopy import *
from netCDF4 import Dataset
import numpy as np
from yambopy.tools.string import marquee
from yambopy.units import I
import shutil
import os

def abs2(x):
    return x.real**2 + x.imag**2

class YamboWFDB():
    """
    Load wavefunctions from yambo (ns.wf)
    
    :: Database: WF[nband,i_sp_pol,ngvect,cmplx] (one db for each k, each sp_pol)

    :: Yambopy: self.wf[nk,nband,nspin,ngvect]

    :: Methods: read(), write(), get_spin_projection() 
    """

    def __init__(self,path=None,save='SAVE',filename='ns.wf'):
        """
        load wavefunction from yambo

        WF[
        """
        if path is None:
            self.path = save
        else:
            self.path = path+f'{save}' # Fix Bug here which made it impossible to read from save with different folder name
        self.filename = filename
        
        #read wf 
        self.read()
        self.nkpoints, self.nbands, self.nspin, self.ng = self.wf.shape

    def read(self):
        path = self.path
        filename = self.filename

        wf = []
        nk = 1
        while True:
            try:
                fname = "%s/%s_fragments_%d_1"%(path,filename,nk)
                database = Dataset(fname)
                aux = database.variables['WF_COMPONENTS_@_SP_POL1_K%d_BAND_GRP_1'%nk][:]
                wf.append( aux[:,:,:,0]+I*aux[:,:,:,1]  )
                nk+=1
            except:
                if nk==1: raise IOError('Could not read %s'%fname)
                break
        self.wf = np.array(wf)

    def write(self,path):
        """
        Write the (new?) wavefunctions in new files
        """
        if os.path.isdir(path): shutil.rmtree(path)
        os.mkdir(path)

        #copy all the files
        oldpath = self.path
        filename = self.filename
        shutil.copyfile("%s/%s"%(oldpath,filename),"%s/%s"%(path,filename))
        for nk in range(self.nkpoints):
            fname = "%s_fragments_%d_1"%(filename,nk+1)
            shutil.copyfile("%s/%s"%(oldpath,fname),"%s/%s"%(path,fname))

        #edit with the new wfs
        wf = self.wf
        for nk in range(self.nkpoints):
            fname = "%s_fragments_%d_1"%(filename,nk+1)
            database = Dataset("%s/%s"%(path,fname),'r+')
            aux = np.array([wf[nk].real,wf[nk].imag])
            aux = np.moveaxis(aux,0,-1)
            database.variables['WF_COMPONENTS_@_SP_POL1_K%d_BAND_GRP_1'%(nk+1)][:] = aux
            database.close()
        print('New wavefunctions written in %s'%path)

    def __str__(self):
        lines = []; app = lines.append
        app(marquee(self.__class__.__name__))
        app("nkpoints: %4d"%self.nkpoints)
        app("nspin:    %4d"%self.nspin)
        app("nbands:   %4d"%self.nbands)
        app("ng:       %4d"%self.ng)
        return "\n".join(lines)

    def get_spin_projections(self,ik,ib):
        """
        By M. Zanfrognini    

        ik : k-point index
        ib : band index

        out: \sum_G|<up|nk>|^2, \sum_G|<down|nk>|^2
        """
        proj_up = np.sum(abs2(self.wf[ik,ib,0]))   
        proj_dn = np.sum(abs2(self.wf[ik,ib,1]))   
        return proj_up, proj_dn

if __name__ == "__main__":
    ywf = YamboWFDB(path='database')
