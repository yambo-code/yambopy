#
# License-Identifier: GPL
#
# Copyright (C) 2024 The Yambo Team
#
# Authors: FP, RR
#
# This file is part of the yambopy project
#
import os
from glob import glob

import numpy as np

from yambopy import *
from yambopy.units import *
from yambopy.plot.plotting import add_fig_kwargs,BZ_Wigner_Seitz
from yambopy.lattice import replicate_red_kmesh, calculate_distances, car_red, red_car
from yambopy.tools.funcs import gaussian, lorentzian
from yambopy.tools.skw import SkwInterpolator
from yambopy import YamboExcitonDB, BrillouinZone

class ExcitonDispersion():
    """
    Class to obtain exciton information at all momenta

    - Dispersion plot (under development)
    - Plots of exciton weights in q-space

    :: Lattice is an instance of YamboLatticeDB
    :: nexcitons is the number of excitonic states - by default it is taken from the Q=1 database

    NB: so far does not support spin-polarised exciton plots (should be implemented when needed!)
    NB: only supports BSEBands option in yambo bse input, not BSEEhEny
    """

    def __init__(self,lattice,nexcitons=None,folder='.'):

        if not isinstance(lattice,YamboLatticeDB):
            raise ValueError('Invalid type for lattice argument. It must be YamboLatticeDB')
        
        files   = glob(folder+'/ndb.BS_diago_Q*')
        nqpoints  = len(files)

        # Check
        if not nqpoints==lattice.ibz_nkpoints :
            raise ValueError("Incomplete list of qpoints (%d/%d)"%(nqpoints,lattice.ibz_nkpoints)) 
    
        dbs_are_consistent, spin_is_there = self.db_check(lattice,nqpoints,folder)
        if nexcitons is None: nexcitons = self.ntransitions

        # Read
        car_qpoints         = np.zeros((nqpoints,3))
        exc_energies        = np.zeros((nqpoints,nexcitons))
        # Commented by JMC to make it work with visual-excitons GUI
        # exc_eigenvectors    = np.zeros((nqpoints,nexcitons,self.ntransitions),dtype=complex)
        # exc_tables          = np.zeros((nqpoints,self.ntransitions,5),dtype=int)
        exc_eigenvectors = np.zeros((nqpoints, nexcitons, self.nkpoints * self.nvalence * self.nconduction), dtype=complex)
        exc_tables = np.zeros((nqpoints, self.nkpoints * self.nvalence * self.nconduction, 5), dtype=int)
        for iQ in range(nqpoints):
            exc_obj = YamboExcitonDB.from_db_file(lattice,filename=folder+'/ndb.BS_diago_Q%d'%(iQ+1))
            if iQ==0: car_qpoints[iQ] = np.array([0.,0.,0.])
            else:     car_qpoints[iQ] = exc_obj.car_qpoint
            exc_energies[iQ,:]        = exc_obj.eigenvalues[:nexcitons].real
            exc_eigenvectors[iQ,:]    = exc_obj.eigenvectors[:nexcitons]
            exc_tables[iQ,:]          = exc_obj.table

        # Set up variables
        self.nqpoints     = nqpoints
        self.nexcitons    = nexcitons
        self.car_qpoints  = car_qpoints
        self.red_qpoints  = car_red(car_qpoints,lattice.rlat)
        self.lattice      = lattice
        self.exc_energies = exc_energies
        self.exc_tables   = exc_tables

        # Reshape eigenvectors if possible
        if dbs_are_consistent and not spin_is_there: self.exc_eigenvectors = self.reshape_eigenvectors(exc_eigenvectors)
        else:                                        self.exc_eigenvectors = exc_eigenvectors

        # Necessary lattice information
        self.alat = lattice.alat
        self.rlat = lattice.rlat

    def db_check(self,lattice,nqpoints,folder):
        """
        Check nexcitons and ntransitions in each database
        """
        nexcitons_each_Q    = np.zeros(nqpoints,dtype=int)
        for iQ in range(nqpoints):
            exc_obj = YamboExcitonDB.from_db_file(lattice,filename=folder+'/ndb.BS_diago_Q%d'%(iQ+1))
            nexcitons_each_Q[iQ] = exc_obj.nexcitons
            if iQ==0: tbl = exc_obj.table

        is_spin_pol = len(np.unique(tbl[:,3]))>1 or len(np.unique(tbl[:,4]))>1
        is_consistent = np.all(nexcitons_each_Q==nexcitons_each_Q[0])
        if not is_consistent: 
            print("[WARNING] BSE Hamiltonian has different dimensions for some Q.")
            print("          Taking the minimum number of transitions, be careful.")
            self.ntransitions = np.min(nexcitons_each_Q)
        else:
            if is_spin_pol: print("[WARNING] Spin-polarised excitons, only partially supported")
            self.ntransitions = nexcitons_each_Q[0]
            valence_bands    = np.unique(tbl[:,1]) - 1
            conduction_bands = np.unique(tbl[:,2]) - 1
            self.nkpoints    = np.max(tbl[:,0])
            self.nvalence    = len(valence_bands)
            self.nconduction = len(conduction_bands)

        return is_consistent, is_spin_pol

    def reshape_eigenvectors(self,eigenvectors):
        """
        eigenvectors in:  [nqpoints, nexcitons, ntransitions]
        eigenvectors out: [nqpoints, nexcitons, nkpoints, nvalence, nconduction]

        TODO: Extend to spin-polarised case
        """
        nq, nexc, nk, nv, nc = self.nqpoints, self.nexcitons, self.nkpoints, self.nvalence, self.nconduction
        reshaped_eigenvectors = np.zeros((nq,nexc,nk,nv,nc),dtype=complex)
        #print(eigenvectors[2,5,2])
        for iQ in range(nq):
            for i_exc in range(nexc): 
                reshaped_eigenvectors[iQ,i_exc,:,:,:] = eigenvectors[iQ,i_exc,:].reshape([nk,nv,nc])
        #print(reshaped_eigenvectors[2,5,0,1,0])    

        return reshaped_eigenvectors

    @add_fig_kwargs
    def plot_Aweights(self,data,plt_show=False,plt_cbar=False,**kwargs):
        """
        2D scatterplot in the q-BZ of the quantity A_{iq}(iexc,ik,ic,iv).

        Any real quantity which is a function of only the q-grid may be supplied.
        The indices iq,inu,ib1,ib2 are user-specified.

        - if plt_show plot is shown
        - if plt_cbar colorbar is shown
        - kwargs example: marker='H', s=300, cmap='viridis', etc.

        NB: So far requires a 2D system.
            Can be improved to plot BZ planes at constant k_z for 3D systems.
        """

        qpts = self.car_qpoints

        # Input check
        if len(data)!=len(qpts):
            raise ValueError('Something wrong in data dimensions (%d data vs %d qpts)'%(len(data),len(qpts)))

        # Global plot stuff
        self.fig, self.ax = plt.subplots(1, 1)
        if self.nqpoints<self.nkpoints:  c_BZ_borders='black'
        if self.nqpoints==self.nkpoints: c_BZ_borders='white'
        self.ax.add_patch(BZ_Wigner_Seitz(self.lattice,color=c_BZ_borders))


        if plt_cbar:
            if 'cmap' in kwargs.keys(): color_map = plt.get_cmap(kwargs['cmap'])
            else:                       color_map = plt.get_cmap('viridis')
        lim = 1.05*np.linalg.norm(self.rlat[0])
        self.ax.set_xlim(-lim,lim)
        self.ax.set_ylim(-lim,lim)

        # Reproduce plot also in adjacent BZs
        BZs = shifted_grids_2D(qpts,self.rlat)
        for qpts_s in BZs: plot=self.ax.scatter(qpts_s[:,0],qpts_s[:,1],c=data,**kwargs)

        if plt_cbar: self.fig.colorbar(plot)

        plt.gca().set_aspect('equal')

        if plt_show: plt.show()
        else: print("Plot ready.\nYou can customise adding savefig, title, labels, text, show, etc...")

    #####################################
    # Dispersion plot under development #
    #####################################
    def get_dispersion(self, bz):
        """
        Obtain the exciton dispersion along the path of bz (a BrillouinZone
        describing this lattice, in any units or cell): the q-points of the
        grid lying on the path.

        Similar to band plots in k-space, check YamboExcitonDB for more comments
        """
        car_qpoints = red_car(self.red_qpoints, self.lattice.rlat)
        exc_qpoints, exc_indexes, _ = bz.find_collinear_lattice_kpoints(self.lattice.lat, car_qpoints)

        # Here assuming same ordering in index expansion between k-yambopy and q-yambo...
        energies = self.exc_energies[self.lattice.kpoints_indexes]
        energies_path  = energies[exc_indexes]

        ybs_disp = YambopyBandStructure(energies_path, exc_qpoints, bz=bz)
        return ybs_disp

    def get_dispersion_interpolated(self):
        """ Interpolated with SkwInterpolator
        """

    def interpolate_dispersion(self, bz: BrillouinZone):
        lpratio = 10
        fermie = 0
        nelect = 0

        ibz_nkpoints = self.lattice.ibz_nkpoints
        ibz_kpoints = np.zeros([ibz_nkpoints, 3])
        for idx_bz, idx_ibz in enumerate(self.lattice.kpoints_indexes):
            ibz_kpoints[idx_ibz] = self.lattice.red_kpoints[idx_bz]

        ibz_energies = self.exc_energies

        na = np.newaxis

        cell = (self.lattice.lat, self.lattice.red_atomic_positions, self.lattice.atomic_numbers)

        symrel = [sym for sym, trev in zip(self.lattice.sym_rec_red, self.lattice.time_rev_list) if trev == False ]
        time_rev = False

        skw = SkwInterpolator(lpratio, ibz_kpoints, ibz_energies[na, :, :], fermie, nelect, cell, symrel, time_rev, verbose=False)

        energies = skw.interp_kpts(bz.kpoints()).eigens

        return bz.kpoints_distances(), np.transpose(energies[0])

    def plot_exciton_disp_ax(self,ax,bz,**kwargs):
        """Plot the exciton dispersion along the path of bz (a BrillouinZone) on ax."""
        ybs_disp = self.get_dispersion(bz)
        return ybs_disp.plot_ax(ax, **kwargs)

    @add_fig_kwargs
    def plot_exciton_disp(self,bz,**kwargs):
        """Plot the exciton dispersion along the path of bz (a BrillouinZone)."""
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        self.plot_exciton_disp_ax(ax,bz)
        return fig
    
    def plot_dispersion(self):
        """ Do plot
        """
        pass
    
    def __str__(self):
        lines = []; app = lines.append
        app(" Exciton Dispersion ")
        app(" Number of qpoints:                    %d"%self.nqpoints)
        app(" Number of exciton branches read:      %d"%self.nexcitons)
        app(" Total number of excitons/transitions: %d"%self.ntransitions)
        return "\n".join(lines)
