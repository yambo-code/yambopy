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
from tqdm import tqdm

def abs2(x):
    return x.real**2 + x.imag**2

class YamboWFDB():
    """
    Load wavefunctions from yambo (ns.wf)
    
    :: Database: WF[nband,i_spinor,ngvect,cmplx] (one db for each k, each sp_pol)

    :: Yambopy: self.wf[nk,nspin,nband,nspinor,ngvect]

    :: Methods: read(), write(), get_spin_projection() 
    """

    def __init__(self,path=None,save='SAVE',filename='ns.wf',bands_range = [], kpts_range=[]):
        """
        load wavefunction from yambo

        bands_range = [min_band, max_band], if left empty, all bands loaded. (Both limits are included)
        kpts_range  = [min_kpts, max_kpts], if left empty, all kpoints (in iBZ) loaded. (Both limits are included)
        
        WF[
        """
        if path is None: path = os.getcwd()
        self.path = os.path.join(path,save)
        self.filename = filename
        
        #read wf 
        self.read()

    def read(self):
        path = self.path
        filename = self.filename

        ## open the ns.db1 database to get essential data 
        try :
            ns_db1 = Dataset(os.path.join(path,'ns.db1'), 'r')
            #
            lat_vec = ns_db1['LATTICE_VECTORS'][...].data
            lat_param = ns_db1['LATTICE_PARAMETER'][...].data
            #
            G_vec = ns_db1['G-VECTORS'][...].data.T # Gvec database
            #
            wfc_grid = ns_db1['WFC_GRID'][...].data
            # gvec indices for each wf in Gvec database
            wfc_grid = np.rint(wfc_grid).astype(int)
            #
            # number of gvecs for each wfc
            igk = np.array(np.rint(ns_db1['WFC_NG'][...].data),dtype=int)
            #
            # COnvert gvec from yambo internal units to cart units
            G_vec = G_vec/lat_param[None,:]
            #
            # Cart to crystal units. In crystal units. Gvecs are integers
            G_vec = np.array(np.rint(G_vec@lat_vec), dtype = int)
            #
            dimensions = ns_db1['DIMENSIONS'][...].data
            self.nkpts_iBZ_total = int(np.rint(dimensions[6]))
            self.nspin     = int(np.rint(dimensions[12]))
            self.nspinor   = int(np.rint(dimensions[11]))
            self.nbnds_total = int(np.rint(dimensions[5]))

            if min(kpts_range) < 1 or max(kpts_range) > self.nkpts_iBZ:
                print("Warning : Wrong input for kpts_range, loading all kpoints")
                kpts_range = [1, self.nkpts_iBZ]
            #
            if min(bands_range) < 1 or max(bands_range) > self.nbnds_total:
                print("Warning : Wrong input for bands_range, loading all bands") 
                bands_range = [1, self.nbnds_total]
            #
            #
            self.nbnds_range = bands_range
            self.nkpts_range = kpts_range
            self.nkpts = max(kpts_range)- min(kpts_range) + 1
            # number of kpoints loaded. not the total number of kpts in ibz
            self.nbnds = max(bands_range) -min(bands_range) + 1
            # number of bnds loaded. not the total number used present in database
            #
            ns_db1.close()
        except :
            raise IOError('Cannot read ns.db1 file')
        ##
        wf = []
        for ik in tqdm(range(self.nkpts), desc="Loading Wfcs "):
            ikk = ik + min(kpts_range)
            for ispin in range(self.nspin):
                try:
                    fname = "%s_fragments_%d_1"%(filename, ispin*self.nkpts_iBZ + ikk)
                    fname = os.path.join(path,fname)
                    database = Dataset(fname,'r')
                    database_var_name = 'WF_COMPONENTS_@_SP_POL%d_K%d_BAND_GRP_1' % (ispin+1, ikk)
                    aux = database.variables[database_var_name][min(bands_range)-1:max(bands_range),...].data
                    aux = aux[...,0] + 1j*aux[...,1]
                    aux[..., igk[ikk-1]:] = 0
                    ## NM : Some of the component of wf's are not valid (i.e they are set to some dummy values) 
                    ## so we set them to be zero.
                    wf.append(aux)
                    database.close()
                except:
                    raise IOError('Could not read %s'%fname)
        
        self.ng = wf[0].shape[-1]
        self.wf = np.array(wf).reshape(self.nkpts, self.nspin, self.nbnds, self.nspinor, self.ng)

## NM : This need's to be fixed But Why do we need this?
    #def write(self,path):
    #    """
    #    Write the (new?) wavefunctions in new files
    #    """
    #    if os.path.isdir(path): shutil.rmtree(path)
    #    os.mkdir(path)

    #    #copy all the files
    #    oldpath = self.path
    #    filename = self.filename
    #    shutil.copyfile("%s/%s"%(oldpath,filename),"%s/%s"%(path,filename))
    #    for nk in range(self.nkpoints):
    #        fname = "%s_fragments_%d_1"%(filename,nk+1)
    #        shutil.copyfile("%s/%s"%(oldpath,fname),"%s/%s"%(path,fname))

    #    #edit with the new wfs
    #    wf = self.wf
    #    for nk in range(self.nkpoints):
    #        fname = "%s_fragments_%d_1"%(filename,nk+1)
    #        database = Dataset("%s/%s"%(path,fname),'r+')
    #        aux = np.array([wf[nk].real,wf[nk].imag])
    #        aux = np.moveaxis(aux,0,-1)
    #        database.variables['WF_COMPONENTS_@_SP_POL1_K%d_BAND_GRP_1'%(nk+1)][:] = aux
    #        database.close()
    #    print('New wavefunctions written in %s'%path)

    def __str__(self):
        lines = []; app = lines.append
        app(marquee(self.__class__.__name__))
        app("nkpoints: %4d - %4d" %(min(self.nkpts_range), max(self.nkpts_range)))
        app("nspin:    %4d"%self.nspin)
        app("nbands:   %4d - %4d" %(min(self.nbnds_range), max(self.nbnds_range)))
        app("ng:       %4d"%self.ng)
        return "\n".join(lines)

    def get_spin_projections(self,ik,ib):
        # Compute the expectation values of S_Z operator for electronic states
        # ik must be in [min_kpts, max_kpts]. not from [0,self.nkpts]
        # ib must be in [min_band, max_band]. not from [0,self.nbands]
        assert self.nspin == 1, "spin projections is useful only for nspin = 1"
        if self.nspinor == 1:
            return 0
        else self.nspinor == 2:
            # <psi| S_z | psi>. +1 for spin up and -1 for spin down. Maybe not be welldefined
            # in the precence of strong spin-orbit coupling.
            s_z = np.array([[1,0],[0,-1]])
            wfc_tmp = self.wf[ik-min(self.nkpts_range), 0, ib-min(self.nbnds_range)]
            s_tmp = np.einsum('ij,jg->ig',s_z, wfc_tmp,optimize=True)
            return np.dot(s_tmp.reshape(-1),wfc_tmp.reshape(-1).conj())
        else:
            print("Invalid nspinor values %d", %(self.nspinor))
            exit()
        

if __name__ == "__main__":
    ywf = YamboWFDB(path='database')
