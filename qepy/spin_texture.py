# 
# License-Identifier: GPL
#
# Copyright (C) 2024 The Yambo Team
#
# Authors: JC-V
# First verion by JC-V (2024)
#
# This file is part of the yambopy project
#
#

from qepy import *
from yambopy import *
from re import findall
from scipy.interpolate import Rbf

class Spin_texture():



    def __init__(self,prefix,path='.',folder_spin='.',prefix_spin='Spin'):

        """ 
        Initialize the structure and compute the spin texture (spin-z as colormap) (complete)

        === Usage and variables ===

        >> Spin_tex = Spin_texture(prefix=db_refix,path=db_path,folder_spin=db_folder_spin,prefix_spin=db_prefix_spin)
        >> Spin_tex.plot_spin_texture(ax,limfactor=0.8,mode="raw",nband=db_nband)

        Input:
        :: prefix is the prefix of the QE directory
        :: path is the path of the QE save folder
        :: folder_spin is the path of the QE spin folder
        :: prefix_spin is the prefix of the spin files
        
        Output:
        :: Plot of the spin texture in three different ways: "raw" where the data is represented as points from QE; "interpolared" where the data is interpolated; and "arrow" where the data is represented as arrows where the direction of the arrows is defined by the S_{x} and S_{y} components.
        """

        print("=== Initializing the data ===")

        self.prefix = prefix
        self.path = path

        # Load data from QE
        self.data_xml = PwXML(prefix=self.prefix, path=self.path)
        self.nkpoints = self.data_xml.nkpoints
        self.nbands = self.data_xml.nbands

        kpoints_cart = np.array(self.data_xml.kpoints)
        kpoints_red = np.vstack([car_red([k], self.data_xml.rcell)[0] for k in kpoints_cart])
        kmesh_full, kmesh_idx = replicate_red_kmesh(kpoints_red, repx=range(-1, 2), repy=range(-1, 2))
        self.x, self.y = red_car(kmesh_full, self.data_xml.rcell)[:, :2].T

        # Load spin data
        self.spin_data = {}
        for i in range(1, 4):
            self.spin_data[i] = self.load_spin_data(folder_spin, prefix_spin, i, kmesh_idx)

    def load_spin_data(self, folder_spin, prefix_spin, index, kmesh_idx):
        """Load spin data from files."""
        filepath = f'{folder_spin}/{prefix_spin}.out.{index}'
        with open(filepath, 'r') as file:
            lines = file.readlines()

        nband = int(findall(r"[-+]?\d*\.\d+|\d+", lines[0].strip().split()[2])[0])
        nk = int(lines[0].strip().split()[-2])
        nline = nband // 10

        if nband < 10 or nband % 10 != 0:
            raise ValueError("Error: nband must be >= 10 and a multiple of 10.")
        if self.nbands != nband or self.nkpoints != nk:
            print("Warning: Dimensions are inconsistent!")

        spin = np.zeros((self.nkpoints, self.nbands))

        for ik in range(self.nkpoints):
            for ib in range(nline):
                ib1, ib2 = int(ib * 10), int((ib + 1) * 10)
                line_idx = int(ik * (nline + 1) + 2 + ib)
                spin[ik, ib1:ib2] = list(map(float, lines[line_idx].split()))
        return spin[kmesh_idx] # I think the problem is here

    def plot_spin_texture(self, ax, limfactor=0.8, mode="raw", nband=1, **kwargs):
        """Plot the spin texture."""
        lim = np.max(self.data_xml.rcell) * limfactor

        if mode == "raw":
            ax.scatter(self.x, self.y, c=self.spin_data[3][:, nband], cmap='seismic', s=20)

        elif mode == "interpolated":
            npts = kwargs.get('npts', 100)
            rbfi = Rbf(self.x, self.y, self.spin_data[3][:, nband], function='linear')
            xi = yi = np.linspace(-lim, lim, npts)
            zi = np.array([[rbfi(x, y) for x in xi] for y in yi])
            ax.imshow(zi, interpolation=kwargs.get('interp_method', 'bicubic'),
                      extent=[-lim, lim, -lim, lim], cmap='seismic')
        elif mode == "arrow":
            ax.quiver(self.x, self.y, self.spin_data[1][:, nband], self.spin_data[2][:, nband],
                      self.spin_data[3][:, nband], cmap='seismic', scale=20)

        else:
            raise ValueError("Invalid mode. Choose 'raw', 'interpolated', or 'arrow'.")

        ax.set_title(f"Spin texture - band = {nband}")
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_xlabel("k$_{x}$ (bohr$^{-1}$)")
        ax.set_ylabel("k$_{y}$ (bohr$^{-1}$)")

        print("=== Spin texture computed successfully ===")
