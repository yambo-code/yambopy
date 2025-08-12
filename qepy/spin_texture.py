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
        Reading of spin-projected values and calculation of the spin texture as a colormap, where the spin-z is represented by the intensity

        === Usage and variables ===

        === Reading spin-projected files ===
        >> Spin_data = Spin_texture(prefix=db_prefix,path=db_path,folder_spin=db_folder_spin,prefix_spin=db_prefix_spin)
        >> Spin_{index} = Spin_data.load_spin_data(folder_spin=db_folder_spin, prefix_spin=db_prefix_spin, index)
 
        === Plot texture ===
        >> Spin_data = Spin_texture(prefix=db_prefix,path=db_path,folder_spin=db_folder_spin,prefix_spin=db_prefix_spin)
        >> Spin_data.plot_spin_texture(ax,limfactor=0.8,mode="raw",nband=db_nband)

        Input:
        :: prefix is the prefix of the QE directory
        :: path is the path of the QE save folder
        :: folder_spin is the path of the QE spin folder
        :: prefix_spin is the prefix of the spin files
        :: index represents 1,2,3 = x,y,z
        
        Output:
        :: Plot of the spin texture in three different ways: 
           - "raw": The data is represented as points directly from QE 
           - "interpolated": The data is interpolated from QE and subsequently represented
           -  "arrow": The data is represented from QE as arrows, where the direction of the arrows is defined by the spin-x and spin-y components
        """

        print("=== Initializing the data ===")

        # Prefix and path of QE data
        self.prefix = prefix
        self.path = path
        self.folder_spin = folder_spin
        self.prefix_spin = prefix_spin

        # Loading data from QE
        self.data_xml = PwXML(prefix=self.prefix, path=self.path)
        self.nkpoints = self.data_xml.nkpoints
        self.nbands = self.data_xml.nbands
  
    def load_spin_data(self, index):

        print(f"=== Reading spin-{index} data from files ===")

        # Reading spin files
        filepath = f'{self.folder_spin}/{self.prefix_spin}.out.{index}'
        with open(filepath, 'r') as file:
            lines = file.readlines()

        # Defining the structure to read the data properly
        nband = int(findall(r"[-+]?\d*\.\d+|\d+", lines[0].strip().split()[2])[0])
        nk = int(lines[0].strip().split()[-2])
        nline = nband // 10

        # Warnings
        if nband < 10 or nband % 10 != 0:
            raise ValueError("Error: nband must be >= 10 and a multiple of 10.")
        if self.nbands != nband or self.nkpoints != nk:
            print("Error: Dimensions are inconsistent!")

        # Definition of spin array
        spin = np.zeros((self.nkpoints, self.nbands))

        # Reading spin values and saving them in IBZ repetition
        for ik in range(self.nkpoints):
            for ib in range(nline):
                ib1, ib2 = int(ib * 10), int((ib + 1) * 10)
                line_idx = int(ik * (nline + 1) + 2 + ib)
                spin[ik, ib1:ib2] = list(map(float, lines[line_idx].split()))

        return spin 

    def plot_spin_texture(self, ax, limfactor=0.8, mode="raw", nband=1, **kwargs):

        print("=== Plotting spin texture ===")

        # Repetition of IBZ
        kpoints_cart = np.array(self.data_xml.kpoints)
        kpoints_red = np.vstack([car_red([k], self.data_xml.rcell)[0] for k in kpoints_cart])
        kmesh_full, kmesh_idx = replicate_red_kmesh(kpoints_red, repx=range(-1, 2), repy=range(-1, 2))
        self.x, self.y = red_car(kmesh_full, self.data_xml.rcell)[:, :2].T

        # Defining plot boundaries
        lim = np.max(self.data_xml.rcell) * limfactor

        # Loading spin-z data
        self.spin_projection = {}
        self.spin_projection[3] = self.load_spin_data(3)

        # Redefining spin-z in IBZ repetition
        self.spin_z = self.spin_projection[3][kmesh_idx]

        # Types of plot
        if mode == "raw":

            ax.scatter(self.x, self.y, c=self.spin_z[:, nband], cmap='seismic', s=20)

        elif mode == "interpolated":

            npts = kwargs.get('npts', 100)
            rbfi = Rbf(self.x, self.y, self.spin_z[:, nband], function='linear')
            xi = yi = np.linspace(-lim, lim, npts)
            zi = np.array([[rbfi(x, y) for x in xi] for y in yi])
            ax.imshow(zi, interpolation=kwargs.get('interp_method', 'bicubic'),
                      extent=[-lim, lim, -lim, lim], cmap='seismic')

        elif mode == "arrow":

            # Loading spin-x and spin-y data
            for index in range(1, 3):
                self.spin_projection[index] = self.load_spin_data(index)

            # Redefining spin-x and spin.y in IBZ repetition
            self.spin_x, self.spin_y = self.spin_projection[1][kmesh_idx], self.spin_projection[2][kmesh_idx]

            ax.quiver(self.x, self.y, self.spin_x[:, nband], self.spin_y[:, nband],
                      self.spin_z[:, nband], cmap='seismic', scale=20)

        else:
            raise ValueError("Invalid mode. Choose 'raw', 'interpolated', or 'arrow'.")

        # Defining plot boundaries
        ax.set_title(f"Spin texture - band = {nband}")
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_xlabel("k$_{x}$ (bohr$^{-1}$)")
        ax.set_ylabel("k$_{y}$ (bohr$^{-1}$)")

        print("=== Spin texture plotted successfully ===")
