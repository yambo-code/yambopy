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
from yambopy.io.cubetools import write_cube
from yambopy.dbs.latticedb import YamboLatticeDB
from scipy.spatial import KDTree

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
            ns_db1_fname = os.path.join(path,'ns.db1')
            self.ydb = YamboLatticeDB.from_db_file(ns_db1_fname,Expand=False)
            #
            ns_db1 = Dataset(ns_db1_fname, 'r')
            #
            lat_vec = ns_db1['LATTICE_VECTORS'][...].data
            self.lat_vec = lat_vec # ai = lat_vec[:,i] 
            lat_param = ns_db1['LATTICE_PARAMETER'][...].data
            #
            ## kpoints in iBZ (crystal units)
            self.kpts_iBZ = self.ydb.iku_kpoints/lat_param[None,:]
            self.kpts_iBZ = self.kpts_iBZ @ lat_vec
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
        return self.to_real_space(wfc_tmp, gvec_tmp, grid = grid)
        
    def write2cube(self, ik, ib, grid=[]):
        wfc_r = self.wfcG2r(ik, ib, grid=grid)
        wfc_r = np.sum(np.abs(wfc_r)**2,axis=(1))
        wfc_r = wfc_r/wfc_r.max()
        for ispin in range(self.nspin):
            filename = 'wfc_k%d_bnd_%d_spin%d.cube' %(ik, ib, ispin+1)
            write_cube(filename, wfc_r[ispin], self.lat_vec, self.ydb.car_atomic_positions, \
                       self.ydb.atomic_numbers, header='Real space electronic wavefunction')

    def get_wfc_kpt(self, ik):
        # return kpoint of ik wfc in crystal units
        assert ik >= min(self.nkpts_range) and ik <= max(self.nkpts_range), \
            "Invalid kpoint idx. Provided out side the range %d - %d" \
            %(min(self.nkpts_range), max(self.nkpts_range))
        return self.kpts_iBZ[ik]

    def rotate_wfc(self, ik, time_rev, sym_mat, frac_vec=np.array([0,0,0])):
        ## Applies a symmetry operaton on wavefuntion.
        ### Given a symmetry operation x -> sym_mat * x + frac_vec 
        ### get the rotated wfc i.e U(R) * psi
        ## both sym_mat and frac_vec are in cartiasian coordinates 

        ### if space == 'g', returns wfc in reciprocal space 
        #### if space = 'r', returns wfc in real space
        
        ## ik is kpoint index in irreducble bz

        ## time_rev == true if sym_mat is a time-reversal symmetry. 
        ## In yambo the standard convertion of representing time reversal is to use -R, -tau
        ##
        ## returns [Rotated kpoints, wfc, gvecs]
        wfc_k, gvecs_k = self.get_wf(ik)
        
        su_mat = np.eye(2)
        if self.nspinor == 1: su_mat = 1
        else : su_mat = su2_mat(sym_mat,time_rev).astype(wfc_k.dtype)
        ### rotate gvecs_k
        sym_red = np.linalg.inv(self.lat_vec) @ sym_mat.T @ self.lat_vec
        sym_red = np.rint(sym_red).astype(int)
        gvec_rot = gvecs_k@sym_red
        
        if self.nspinor == 2:
            wfc_rot = np.einsum('ij,sbjg->sbig',su_mat,wfc_k,optimize=True)

        ## add phase due to fractional translation
        Rkvec = sym_red.T@self.kpts_iBZ[ik]
        tau_frac = lat@frac_vec

        kphase = np.exp( -1j * 2*np.pi * np.dot(Rkvec, tau_frac))
        gphase = kphase*np.exp( -1j * 2*np.pi * (gvec_rot@tau_frac))

        wfc_rot *= gphase[None,None,None,:]
        if time_rev: wfc_rot = wfc_rot.conj()

        return [Rkvec, wfc_rot,gvec_rot]

    
    def to_real_space(self, wfc_tmp, gvec_tmp, grid = []):
        #
        # COnvert a wfc from Gspace to real space 
        # wfc_tmp C_G
        # gvec_tmp miller indices
        # grid fft grid (if emtpy, default will used)
        #
        if len(grid) == 0: grid = self.fft_box
        cel_vol = abs(np.linalg.det(self.lat_vec))
        #
        tmp_wfc = np.zeros((self.nspin,self.nspinor, grid[0], grid[1], grid[2]),dtype=wfc_tmp.dtype)
        #
        Nx_vals = np.where(gvec_tmp[:,0] >= 0, gvec_tmp[:,0], gvec_tmp[:,0] + grid[0])
        Ny_vals = np.where(gvec_tmp[:,1] >= 0, gvec_tmp[:,1], gvec_tmp[:,1] + grid[1])
        Nz_vals = np.where(gvec_tmp[:,2] >= 0, gvec_tmp[:,2], gvec_tmp[:,2] + grid[2])

        assert np.min(Nx_vals) >=0 and (np.min(Ny_vals) >=0 and np.min(Nz_vals)) >=0, "Wrong fft indices"
        tmp_wfc[:,:,Nx_vals, Ny_vals,Nz_vals] = wfc_tmp
        return scipy.fft.ifftn(tmp_wfc,norm="forward",axes=(2,3,4))/np.sqrt(cel_vol)



def wfc_inner_product(k_bra, wfc_bra, gvec_bra, k_ket, wfc_ket, gvec_ket):
    ## k1, k2, are kpoint indices in crystal coordinates
    ## wfc1 and wfc2 are wfc coeffiencts (nspin, nbnd, nspinor, ng)
    ## gvec1 and gvec2 are miller indices (ng,3)
    ## < k2 | k1>
    #
    assert wfc_ket.shape[:3] == wfc_bra.shape[:3], "Inconsistant wfcs"
    #
    nspin, nbnd, nspinor = wfc_ket.shape[:3]
    kdiff = k_bra-k_ket
    G0 = np.rint(kdiff).astype(int)
    kdiff = kdiff-G0
    ## crystal momentum mismatch delta_{k,k'}
    if np.max(np.abs(kdiff)) > 1e-5:
        return np.zeros((nspin, nbnd, nbnd),dtype=wfc_ket.dtype)
    ket_Gtree = KDTree(gvec_ket)
    gbra_shift = gvec_bra + G0[None,:]
    dd, ii = ket_Gtree.query(gbra_shift, k=1)
    wfc_bra_tmp = np.zero(wfc_ket.shape,dtype=wfc_ket.dtype)
    bra_idx = ii[dd < 1e-6]
    wfc_bra_tmp[:,:,:,bra_idx] = wfc_bra[...,dd<1e-6].conj()
    return np.einsum('sixg,sjxg->sij',wfc_bra_tmp,wfc_ket,optimize=True)
    


def su2_mat(symm_mat,time_rev1):
    # Computes the Spinor rotation matrix
    if np.isclose(np.abs(np.trace(symm_mat)),3,atol=1e-03):
        su_mat = np.eye(2)
    else :
        sym_temp = symm_mat/np.linalg.det(symm_mat)
        
        if abs(sym_temp[2,0] +1) < 10**-5 :
            alpha = 0
            beta = -np.pi/2
            delta = np.arctan2(sym_temp[0,1],sym_temp[0,2])
        elif abs(sym_temp[2,0] - 1) < 10**-5:
            alpha = 0
            beta = np.pi/2
            delta = np.arctan2(-sym_temp[0,1],-sym_temp[0,2])
        else :
            beta = np.arcsin(sym_temp[2,0])
            alpha = np.arctan2(sym_temp[1,0]/np.cos(beta),sym_temp[0,0]/np.cos(beta))
            delta = np.arctan2(sym_temp[2,1]/np.cos(beta),sym_temp[2,2]/np.cos(beta))
        # print(180/np.pi*alpha,180/np.pi*beta,180/np.pi*delta)
        alpha = alpha/2 ; beta = -beta/2 ; delta = delta/2 ;
        spin_RX_delta=sigma_0*np.cos(delta)-1j*sigma_x*np.sin(delta)
        spin_RY_beta =sigma_0*np.cos(beta) -1j*sigma_y*np.sin(beta)
        spin_RZ_alpha=sigma_0*np.cos(alpha)-1j*sigma_z*np.sin(alpha)
        su_mat = spin_RZ_alpha@spin_RY_beta@spin_RX_delta
    if time_rev1:
        ## -ve sign because of yambo does it wrong
        su_mat = -np.array([[0, -1], [1, 0]])@su_mat
    return su_mat


if __name__ == "__main__":
    ywf = YamboWFDB(path='database')
