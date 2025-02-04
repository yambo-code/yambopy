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
import scipy

class YamboWFDB:
    """
    A class to load and manipulate wavefunctions from Yambo (ns.wf).

    Attributes:
        path (str): Path to the directory containing the wavefunction files.
        filename (str): Name of the wavefunction file (default: 'ns.wf').
        wf (numpy.ndarray): Wavefunctions stored as a 5D array [nkpoints, nspin, nbands, nspinor, ngvect].
        gvecs (numpy.ndarray): G-vectors for the wavefunctions [nkpoints, ngvect, 3].
        kpts_iBZ (numpy.ndarray): K-points in the irreducible Brillouin Zone (iBZ) in crystal coordinates.
        ngvecs (numpy.ndarray): Number of meaningful G-vectors for each wavefunction.
        fft_box (numpy.ndarray): Default FFT grid size for real-space conversion.
        nkpoints (int): Number of k-points in the iBZ.
        nspin (int): Number of spin components.
        nspinor (int): Number of spinor components.
        nbands (int): Number of bands.
        min_bnd (int): Minimum band index loaded.

    Methods:
        read(bands_range=[]): Read wavefunctions from the file.
        get_spin_projections(ik, ib, s_z=np.array([[1, 0], [0, -1]])): Compute spin projections.
        get_wf(ik): Get wavefunctions and G-vectors for a specific k-point.
        wfcG2r(ik, ib, grid=[]): Convert wavefunctions from G-space to real space.
        write2cube(ik, ib, grid=[]): Write wavefunctions to a cube file.
        get_wfc_kpt(ik): Get the k-point in crystal coordinates.
        rotate_wfc(ik, isym): Rotate wavefunctions using a symmetry operation.
        apply_symm(kvec, wfc_k, gvecs_k, time_rev, sym_mat, frac_vec=np.array([0, 0, 0])): Apply symmetry to wavefunctions.
        to_real_space(wfc_tmp, gvec_tmp, grid=[]): Convert wavefunctions to real space.
    """

    def __init__(self, path=None, save='SAVE', filename='ns.wf', bands_range=[]):
        """
        Initialize the YamboWFDB class.

        Args:
            path (str, optional): Path to the directory containing the wavefunction files. Defaults to the current directory.
            save (str, optional): Subdirectory containing the wavefunction files. Defaults to 'SAVE'.
            filename (str, optional): Name of the wavefunction file. Defaults to 'ns.wf'.
            bands_range (list, optional): Range of bands to load. Defaults to all bands.
        """
        if path is None:
            path = os.getcwd()
        self.path = os.path.join(path, save)
        self.filename = filename

        # Read wavefunctions
        self.read(bands_range=bands_range)

    def read(self, bands_range=[]):
        """
        Read wavefunctions from the file.

        Args:
            bands_range (list, optional): Range of bands to load. Defaults to all bands.
        """
        path = self.path
        filename = self.filename

        # Open the ns.db1 database to get essential data
        try:
            ns_db1_fname = os.path.join(path, 'ns.db1')
            self.ydb = YamboLatticeDB.from_db_file(ns_db1_fname, Expand=False)

            # Read G-vectors and other data from ns_db1
            ns_db1 = Dataset(ns_db1_fname, 'r')
            lat_vec = self.ydb.lat.T
            lat_param = self.ydb.alat

            # K-points in iBZ (crystal units)
            self.kpts_iBZ = self.ydb.iku_kpoints / lat_param[None, :]
            self.kpts_iBZ = self.kpts_iBZ @ lat_vec

            # G-vectors in cartesian units
            G_vec = ns_db1['G-VECTORS'][...].data.T
            wfc_grid = ns_db1['WFC_GRID'][...].data
            wfc_grid = np.rint(wfc_grid).astype(int)

            # Number of G-vectors for each wavefunction
            igk = np.array(np.rint(ns_db1['WFC_NG'][...].data), dtype=int)
            self.ngvecs = igk

            # Convert G-vectors to crystal units
            G_vec = G_vec / lat_param[None, :]
            G_vec = np.array(np.rint(G_vec @ lat_vec), dtype=int)

            # Determine the default FFT grid size
            min_fft_idx = np.min(G_vec, axis=0)
            max_fft_idx = np.max(G_vec, axis=0)
            assert np.min(max_fft_idx) >= 0 and np.max(min_fft_idx) < 0, "Invalid G-vectors"
            self.fft_box = np.zeros(3, dtype=int)
            for i in range(3):
                self.fft_box[i] = max_fft_idx[i] - min_fft_idx[i] + 1

            # Read dimensions
            dimensions = ns_db1['DIMENSIONS'][...].data
            self.nkpoints = int(np.rint(dimensions[6]))  # Number of k-points in iBZ
            self.nspin = int(np.rint(dimensions[12]))    # Number of spin components
            self.nspinor = int(np.rint(dimensions[11]))  # Number of spinor components
            self.nbands = int(np.rint(dimensions[5]))    # Number of bands

            # Validate bands_range
            if len(bands_range) == 0:
                bands_range = [0, self.nbands]
            elif min(bands_range) < 0 or max(bands_range) > self.nbands:
                print("Warning: Invalid bands_range, loading all bands.")
                bands_range = [0, self.nbands]

            self.min_bnd = min(bands_range)
            self.nbands = max(bands_range) - self.min_bnd

            ns_db1.close()
        except Exception as e:
            raise IOError(f'Cannot read ns.db1 file: {e}')

        # Load wavefunctions
        wf = []
        for ik in tqdm(range(self.nkpoints), desc="Loading Wavefunctions"):
            for ispin in range(self.nspin):
                try:
                    fname = f"{filename}_fragments_{ispin * self.nkpoints + ik + 1}_1"
                    fname = os.path.join(path, fname)
                    database = Dataset(fname, 'r')
                    database_var_name = f'WF_COMPONENTS_@_SP_POL{ispin + 1}_K{ik + 1}_BAND_GRP_1'
                    aux = database.variables[database_var_name][self.min_bnd:self.min_bnd + self.nbands, ...].data
                    aux = aux[..., 0] + 1j * aux[..., 1]
                    aux[..., igk[ik]:] = 0  # Set invalid components to zero
                    wf.append(aux)
                    database.close()
                except Exception as e:
                    raise IOError(f'Could not read {fname}: {e}')

        self.ng = wf[0].shape[-1]  # Maximum number of wavefunction components
        self.wf = np.array(wf).reshape(self.nkpoints, self.nspin, self.nbands, self.nspinor, self.ng)

        # Load G-vectors
        self.gvecs = np.zeros((self.nkpoints, self.ng, 3), dtype=int)
        for ik in tqdm(range(self.nkpoints), desc="Loading Miller Indices"):
            self.gvecs[ik, igk[ik]:, :] = 2147483646 * np.array([1, 1, 1])[None, :]  # Invalid G-vector marker
            self.gvecs[ik, :igk[ik], :] = G_vec[wfc_grid[ik][:igk[ik]] - 1, :]

    def __str__(self):
        """Return a string representation of the object."""
        lines = []
        lines.append(marquee(self.__class__.__name__))
        lines.append(f"nkpoints: {self.nkpoints:4d}")
        lines.append(f"nspin:    {self.nspin:4d}")
        lines.append(f"nbands:   [{self.min_bnd:4d}, {self.min_bnd + self.nbands:4d})")
        lines.append(f"Max wf components: {self.ng:4d}")
        return "\n".join(lines)

    def assert_k_inrange(self, ik):
        """Assert that the k-point index is valid."""
        assert 0 <= ik < self.nkpoints, "Invalid k-point index"

    def assert_bnd_range(self, ib):
        """Assert that the band index is valid."""
        assert 0 <= ib < self.nbands, "Invalid band index"

    def get_spin_projections(self, ik, ib, s_z=np.array([[1, 0], [0, -1]])):
        """
        Compute the expectation values of the spin operator for electronic states.

        Args:
            ik (int): K-point index.
            ib (int): Band index.
            s_z (numpy.ndarray, optional): Spin operator. Defaults to the z-component.

        Returns:
            numpy.ndarray: Spin projection values.
        """
        self.assert_k_inrange(ik)
        self.assert_bnd_range(ib)

        assert self.nspin == 1 and self.nspinor == 2, "Spin projections are only useful for nspin=1 and nspinor=2"
        assert s_z.shape == (2, 2), "Spin operator must be a 2x2 matrix"

        wfc_tmp = self.wf[ik, 0, ib]
        s_tmp = np.einsum('ij,jg->ig', s_z, wfc_tmp, optimize=True)
        return np.dot(s_tmp.reshape(-1), wfc_tmp.reshape(-1).conj())

    def get_wf(self, ik):
        """
        Get wavefunctions and G-vectors for a specific k-point.

        Args:
            ik (int): K-point index.

        Returns:
            list: Wavefunctions and G-vectors.
        """
        self.assert_k_inrange(ik)
        return [self.wf[ik][..., :self.ngvecs[ik]], self.gvecs[ik, :self.ngvecs[ik], :]]

    def wfcG2r(self, ik, ib, grid=[]):
        """
        Convert wavefunctions from G-space to real space.

        Args:
            ik (int): K-point index.
            ib (int): Band index.
            grid (list, optional): FFT grid size. Defaults to the default FFT box.

        Returns:
            numpy.ndarray: Wavefunctions in real space.
        """
        assert len(grid) == 0 or len(grid) == 3, "Grid must be an empty list or a list of 3 integers"
        self.assert_k_inrange(ik)
        self.assert_bnd_range(ib)

        if len(grid) == 0:
            grid = self.fft_box
        else:
            for i in range(3):
                assert grid[i] >= self.fft_box[i], f"Invalid FFT grid. Grid must be >= {self.fft_box[i]}"

        print(f'FFT Grid: {grid[0]} {grid[1]} {grid[2]}')

        wfc_tmp = self.wf[ik][:, ib, :, :self.ngvecs[ik]]
        gvec_tmp = self.gvecs[ik, :self.ngvecs[ik], :]
        return self.to_real_space(wfc_tmp, gvec_tmp, grid=grid)

    def write2cube(self, ik, ib, grid=[]):
        """
        Write wavefunctions to a cube file.

        Args:
            ik (int): K-point index.
            ib (int): Band index.
            grid (list, optional): FFT grid size. Defaults to the default FFT box.
        """
        wfc_r = self.wfcG2r(ik, ib, grid=grid)
        wfc_r = np.sum(np.abs(wfc_r) ** 2, axis=(1))
        wfc_r = wfc_r / wfc_r.max()
        for ispin in range(self.nspin):
            filename = f'wfc_k{ik + 1}_bnd_{ib + self.min_bnd + 1}_spin{ispin + 1}.cube'
            write_cube(filename, wfc_r[ispin], self.ydb.lat.T, self.ydb.car_atomic_positions,
                       self.ydb.atomic_numbers, header='Real space electronic wavefunction')

    def get_wfc_kpt(self, ik):
        """
        Get the k-point in crystal coordinates.

        Args:
            ik (int): K-point index.

        Returns:
            numpy.ndarray: K-point in crystal coordinates.
        """
        self.assert_k_inrange(ik)
        return self.kpts_iBZ[ik]

    def rotate_wfc(self, ik, isym):
        """
        Rotate wavefunctions using a symmetry operation.

        Args:
            ik (int): K-point index.
            isym (int): Symmetry operation index.

        Returns:
            list: Rotated wavefunctions and G-vectors.
        """
        wfc_k, gvecs_k = self.get_wf(ik)
        sym_mat = self.ydb.sym_car[isym]
        kvec = self.get_wfc_kpt(ik)
        time_rev = (isym >= len(self.ydb.sym_car) / (1 + int(np.rint(self.ydb.time_rev))))
        return self.apply_symm(kvec, wfc_k, gvecs_k, time_rev, sym_mat)

    def apply_symm(self, kvec, wfc_k, gvecs_k, time_rev, sym_mat, frac_vec=np.array([0, 0, 0])):
        """
        Apply symmetry to wavefunctions.

        Args:
            kvec (numpy.ndarray): K-point in crystal coordinates.
            wfc_k (numpy.ndarray): Wavefunctions for the k-point.
            gvecs_k (numpy.ndarray): G-vectors for the k-point.
            time_rev (bool): Whether to apply time reversal.
            sym_mat (numpy.ndarray): Symmetry matrix.
            frac_vec (numpy.ndarray, optional): Fractional translation vector. Defaults to [0, 0, 0].

        Returns:
            list: Rotated wavefunctions and G-vectors.
        """
        sym_red = np.linalg.inv(self.ydb.lat.T) @ sym_mat.T @ self.ydb.lat.T
        sym_red = np.rint(sym_red).astype(int)
        gvec_rot = gvecs_k @ sym_red

        if self.nspinor == 2:
            su_mat = su2_mat(sym_mat, time_rev).astype(wfc_k.dtype)
            wfc_rot = np.einsum('ij,sbjg->sbig', su_mat, wfc_k, optimize=True)

        # Add phase due to fractional translation
        Rkvec = sym_red.T @ kvec
        tau_frac = self.ydb.lat.T @ frac_vec
        kphase = np.exp(-1j * 2 * np.pi * np.dot(Rkvec, tau_frac))
        gphase = kphase * np.exp(-1j * 2 * np.pi * (gvec_rot @ tau_frac))
        wfc_rot *= gphase[None, None, None, :]

        if time_rev:
            wfc_rot = wfc_rot.conj()

        return [Rkvec, wfc_rot, gvec_rot]

    def to_real_space(self, wfc_tmp, gvec_tmp, grid=[]):
        """
        Convert wavefunctions from G-space to real space.

        Args:
            wfc_tmp (numpy.ndarray): Wavefunctions in G-space.
            gvec_tmp (numpy.ndarray): G-vectors.
            grid (list, optional): FFT grid size. Defaults to the default FFT box.

        Returns:
            numpy.ndarray: Wavefunctions in real space.
        """
        if len(grid) == 0:
            grid = self.fft_box
        cel_vol = abs(np.linalg.det(self.ydb.lat.T))

        tmp_wfc = np.zeros((self.nspin, self.nspinor, grid[0], grid[1], grid[2]), dtype=wfc_tmp.dtype)

        Nx_vals = np.where(gvec_tmp[:, 0] >= 0, gvec_tmp[:, 0], gvec_tmp[:, 0] + grid[0])
        Ny_vals = np.where(gvec_tmp[:, 1] >= 0, gvec_tmp[:, 1], gvec_tmp[:, 1] + grid[1])
        Nz_vals = np.where(gvec_tmp[:, 2] >= 0, gvec_tmp[:, 2], gvec_tmp[:, 2] + grid[2])
        #
        assert np.min(Nx_vals) >=0 and (np.min(Ny_vals) >=0 and np.min(Nz_vals)) >=0, "Wrong fft indices"
        tmp_wfc[:,:,Nx_vals, Ny_vals,Nz_vals] = wfc_tmp
        return scipy.fft.ifftn(tmp_wfc,norm="forward",axes=(2,3,4))/np.sqrt(cel_vol)



def wfc_inner_product(k_bra, wfc_bra, gvec_bra, k_ket, wfc_ket, gvec_ket, ket_Gtree=-1):
    """
    Computes the inner product between two wavefunctions in reciprocal space. <k_bra | k_ket>
    
    Parameters
    ----------
    k_bra : ndarray
        Crystal momentum of the bra wavefunction (3,) in reduced coordinates.
    wfc_bra : ndarray
        Wavefunction coefficients for the bra state with shape (nspin, nbnd, nspinor, ng).
    gvec_bra : ndarray
        Miller indices of the bra wavefunction (ng, 3) in reduced coordinates.
    k_ket : ndarray
        Crystal momentum of the ket wavefunction (3,) in reduced coordinates.
    wfc_ket : ndarray
        Wavefunction coefficients for the ket state with shape (nspin, nbnd, nspinor, ng).
    gvec_ket : ndarray
        Miller indices of the ket wavefunction (ng, 3) in reduced coordinates.
    ket_Gtree  : scipy.spatial._kdtree.KDTree (optional)
        Kdtree for gvec_ket. leave it or give -1 to internally build one
    #
    Returns
    -------
    ndarray
        Inner product matrix of shape (nspin, nbnd, nbnd). If the momenta mismatch
        is too large, returns a zero matrix.
    """
    #
    # Check consistency of wavefunction dimensions
    assert wfc_ket.shape[:3] == wfc_bra.shape[:3], "Inconsistant wfcs"
    #
    nspin, nbnd, nspinor = wfc_ket.shape[:3]
    kdiff = k_bra-k_ket
    G0 = np.rint(kdiff).astype(int)
    kdiff = kdiff-G0
    ## crystal momentum mismatch delta_{k,k'}
    if np.max(np.abs(kdiff)) > 1e-5:
        return np.zeros((nspin, nbnd, nbnd),dtype=wfc_ket.dtype)
    # Construct KDTree for nearest-neighbor search in G-vectors
    if type(ket_Gtree) != scipy.spatial._kdtree.KDTree:
        ket_Gtree = KDTree(gvec_ket)
    gbra_shift = gvec_bra + G0[None,:]
    ## get the nearest indices and their distance
    dd, ii = ket_Gtree.query(gbra_shift, k=1)
    #
    wfc_bra_tmp = np.zeros(wfc_ket.shape,dtype=wfc_ket.dtype)
    # Get only the indices that are present
    bra_idx = ii[dd < 1e-6]
    #
    wfc_bra_tmp[:,:,:,bra_idx] = wfc_bra[...,dd<1e-6].conj()
    # return the dot product
    return np.einsum('sixg,sjxg->sij',wfc_bra_tmp,wfc_ket,optimize=True)



def su2_mat(symm_mat,time_rev1):
    """
    Computes the SU(2) spinor rotation matrix for a given symmetry matrix.

    Parameters
    ----------
    symm_mat : ndarray
        3x3 rotation matrix describing the symmetry operation.
    time_rev1 : bool
        If True, applies time-reversal symmetry to the SU(2) matrix.

    Returns
    -------
    ndarray
        2x2 SU(2) matrix representing the spinor transformation.
    """
    # If the symmetry matrix corresponds to identity or inversion, return identity matrix
    if np.isclose(np.abs(np.trace(symm_mat)),3,atol=1e-03):
        su_mat = np.eye(2)
    else :
        sigma_0= np.array([1, 0 , 0, 1]).reshape(2,2)
        sigma_x= np.array([0, 1 , 1, 0]).reshape(2,2)
        sigma_y= np.array([0,-1j, 1j,0]).reshape(2,2)
        sigma_z= np.array([1 , 0, 0,-1]).reshape(2,2)
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
        #
        # Convert angles to half-angles for spinor representation
        alpha = alpha/2 ; beta = -beta/2 ; delta = delta/2 ;
        #
        # Compute SU(2) rotation matrix using Euler angle decomposition
        spin_RX_delta=sigma_0*np.cos(delta)-1j*sigma_x*np.sin(delta)
        spin_RY_beta =sigma_0*np.cos(beta) -1j*sigma_y*np.sin(beta)
        spin_RZ_alpha=sigma_0*np.cos(alpha)-1j*sigma_z*np.sin(alpha)
        su_mat = spin_RZ_alpha@spin_RY_beta@spin_RX_delta
    # Apply time-reversal symmetry if needed
    if time_rev1:
        su_mat = np.array([[0, -1], [1, 0]])@su_mat
    return su_mat

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


if __name__ == "__main__":
    ywf = YamboWFDB(path='database')



