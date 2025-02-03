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
import scipy.fft

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
        self.nbnds_range = bands_range
        self.nkpts_range = kpts_range
        #read wf 
        self.read()

    def read(self):
        path = self.path
        filename = self.filename

        bands_range = self.nbnds_range
        kpts_range = self.nkpts_range

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
            self.ngvecs = igk ## number of meaning full gvecs for each wfc
            #
            # COnvert gvec from yambo internal units to cart units
            G_vec = G_vec/lat_param[None,:]
            #
            # Cart to crystal units. In crystal units. Gvecs are integers
            G_vec = np.array(np.rint(G_vec@lat_vec), dtype = int)
            #
            ## find a default fft box grid 
            min_fft_idx = np.min(G_vec,axis=0)
            max_fft_idx = np.max(G_vec,axis=0)
            assert np.min(max_fft_idx) >= 0 and np.max(min_fft_idx) < 0, "Invalid Gvectors"
            self.fft_box = np.zeros(3,dtype=int)
            for i in range(3):
                self.fft_box[i] = max_fft_idx[i] - min_fft_idx[i] + 1
            #
            dimensions = ns_db1['DIMENSIONS'][...].data
            self.nkpts_iBZ_total = int(np.rint(dimensions[6]))
            self.nspin     = int(np.rint(dimensions[12]))
            self.nspinor   = int(np.rint(dimensions[11]))
            self.nbnds_total = int(np.rint(dimensions[5]))

            if len(kpts_range) == 0: kpts_range = [1, self.nkpts_iBZ_total]
            if len(bands_range) == 0: bands_range = [1, self.nbnds_total]

            if min(kpts_range) < 1 or max(kpts_range) > self.nkpts_iBZ_total:
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
                    fname = "%s_fragments_%d_1"%(filename, ispin*self.nkpts_iBZ_total + ikk)
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

        self.gvecs = np.zeros((self.nkpts,self.ng,3),dtype=int) # gvecs for the wfcs [ik, ng]
        for ik in tqdm(range(self.nkpts), desc="Loading miller indices "):
            ikk = ik + min(kpts_range)
            self.gvecs[ik,igk[ikk-1]:,:] = 2147483646*np.array([1,1,1])[None,:]
            ###  NM :This is a garbage number used to identify 
            ### that this is invalid gvec and should be ignored
            self.gvecs[ik,:igk[ikk-1],:] = G_vec[wfc_grid[ikk-1][:igk[ikk-1]]-1,:]

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

        assert ik >= min(self.nkpts_range) and ik <= max(self.nkpts_range), \
            "Invalid kpoint idx. Provided out side the range %d - %d" \
            %(min(self.nkpts_range), max(self.nkpts_range))

        assert ib >= min(self.nbnds_range) and ib <= max(self.nbnds_range), \
            "Invalid band idx. Provided out side the range %d - %d" \
            %(min(self.nbnds_range), max(self.nbnds_range))

        assert self.nspin == 1, "spin projections is useful only for nspin = 1"
        if self.nspinor == 1:
            return 0
        elif self.nspinor == 2:
            # <psi| S_z | psi>. +1 for spin up and -1 for spin down. Maybe not be welldefined
            # in the precence of strong spin-orbit coupling.
            s_z = np.array([[1,0],[0,-1]])
            wfc_tmp = self.wf[ik-min(self.nkpts_range), 0, ib-min(self.nbnds_range)]
            s_tmp = np.einsum('ij,jg->ig',s_z, wfc_tmp,optimize=True)
            return np.dot(s_tmp.reshape(-1),wfc_tmp.reshape(-1).conj())
        else:
            print("Invalid nspinor values %d" %(self.nspinor))
            exit()

    def get_wf(self,ik):
        ## returns :
        ## wfc and its gvecs in crystal coordinates
        ## ik is in ibz 
        ## wfc (nspin,nbnds, nspinor, ng). note ng <= self.ng (only valid gvecs are returned)
        ## 
        assert ik >= min(self.nkpts_range) and ik <= max(self.nkpts_range), \
            "Invalid kpoint idx. Provided out side the range %d - %d" \
            %(min(self.nkpts_range), max(self.nkpts_range))
        ikk = ik-min(self.nkpts_range)
        return [self.wf[ikk][...,:self.ngvecs[ik-1]], self.gvecs[ikk,:self.ngvecs[ik-1],:] ]

    def wfcG2r(self, ik, ib, grid=[]):
        ## grid is list of three number. (fft box)
        ## get the wavefunction on real space. (nspin, nspinor, Nx, Ny, Nz)
        ## The default is a box that will barely fit the sphere.
        ## Incase, the user wants to provide something, it must be greater than default box
        ### default box = self.fft_box
        assert len(grid) == 0 or len(grid) == 3, "grid must be empty list or list of 3 integers"
        assert ik >= min(self.nkpts_range) and ik <= max(self.nkpts_range), \
            "Invalid kpoint idx. Provided out side the range %d - %d" \
            %(min(self.nkpts_range), max(self.nkpts_range))

        assert ib >= min(self.nbnds_range) and ib <= max(self.nbnds_range), \
            "Invalid band idx. Provided out side the range %d - %d" \
            %(min(self.nbnds_range), max(self.nbnds_range))
        #
        if len(grid) == 0:
            grid = self.fft_box
        else :
            for i in range(3):
                assert grid[i] >= self.fft_box[i], \
                    "Invalid fft grid. FFT grid must be >= %d" %(self.fft_box[i])
        print('FFT Grid : %d %d %d' %(grid[0], grid[1], grid[2]))
        ikk = ik-min(self.nkpts_range)
        ibb = ib-min(self.nbnds_range)
        wfc_tmp = self.wf[ikk][:,ibb,:,:self.ngvecs[ik-1]]
        gvec_tmp = self.gvecs[ikk,:self.ngvecs[ik-1],:]
        #
        tmp_wfc = np.zeros((self.nspin,self.nspinor, grid[0], grid[1], grid[2]),dtype=wfc_tmp.dtype)
        #
        Nx_vals = np.where(gvec_tmp[:,0] >= 0, gvec_tmp[:,0], gvec_tmp[:,0] + grid[0])
        Ny_vals = np.where(gvec_tmp[:,1] >= 0, gvec_tmp[:,1], gvec_tmp[:,1] + grid[1])
        Nz_vals = np.where(gvec_tmp[:,2] >= 0, gvec_tmp[:,2], gvec_tmp[:,2] + grid[2])
        
        assert np.min(Nx_vals) >=0 and (np.min(Ny_vals) >=0 and np.min(Nz_vals)) >=0, "Wrong fft indices"
        tmp_wfc[:,:,Nx_vals, Ny_vals,Nx_vals] = wfc_tmp
        return scipy.fft.fftn(tmp_wfc,axes=(2,3,4))

if __name__ == "__main__":
    ywf = YamboWFDB(path='database')
