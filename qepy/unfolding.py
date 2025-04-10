# 
# License-Identifier: GPL
#
# Copyright (C) 2024 The Yambo Team
#
# Authors: HM, AM-S, JC-V, FP
# First verion by HM and AM-S (2017), revised and expanded by JC-V (2025)
#
# This file is part of the yambopy project
#
#
from qepy import *
import numpy as np
from sys import stdout
import h5py

class Unfolding():

    def __init__(self,prefix_pc,prefix_sc,path_pc='.',path_sc='.',spin="none",band_min=0,sc_rotated=False,compute_projections=True):
        """ 
        Initialize the structure with data from the primitive cell and the supercell and then compute the projection of the supercell band structure
        taking the primitive cell band structure as a reference

        === Usage and variables ===

        >> Unfold = Unfolding(prefix_pc=db_prefix_pc,path_pc=db_path_pc,prefix_sc=db_prefix_sc,path_sc=db_path_sc,spin="noncol",band_min = 0,sc_rotated=False,compute_projections=True)
        >> Unfold.plot_eigen_ax(ax,path=db_path_kpoints_sc,ylim = (ylim_min,ylim_max))
 
        Input:
        :: prefix_pc(sc) is the prefix of the primitive cell (supercell) QE save
        :: path_pc(sc) is the path of the primitive cell (supercell) along the Brillouin zone
        :: spin is equal to "none" or "noncol" for nspin = 1 (noSOC) and nspin 4 (SOC) calculations according to QE format
        :: band min determines the number of computed bands. If band_min = 0, it computes the projection of all the bands. If band_min != 0, it computes the projection of the nbands-band_min bands
        :: sc_rotated determines the rotation matrix if necessary. If sc_rotated = False, it takes the matrix identity. If sc_rotated = (3x3) matrix, it takes the rotation matrix defined by the user
        :: compute_projections determines whether to calculate projections or not. If compute_projections = True, it calculates the projections and save them in a .npy file.
           If compute_projections = False, it loads the previous .npy file.
        
        Output:
        :: projections.npy file to plot the EBS with an external script
        :: Plot with the EBS of the supercell taking as a reference the defined primitive cell
        """

        print("=== Initializing the data ===")

        # Prefix and path of primitive cell and supercell
        self.prefix_pc = prefix_pc  
        self.prefix_sc = prefix_sc
        self.path_pc   = path_pc 
        self.path_sc   = path_sc
       
        # Reading primitive cell and supercell database from QE
        pc_xml = PwXML(prefix=self.prefix_pc,path=self.path_pc) 
        sc_xml = PwXML(prefix=self.prefix_sc,path=self.path_sc)

        # Number of kpoints of the primitive cell and supercell
        self.nkpoints_pc = pc_xml.nkpoints 
        self.nkpoints_sc = sc_xml.nkpoints 

        # List of kpoints of the supercell 
        self.kpoints = sc_xml.kpoints  

        # Number of bands of the primitive cell and supercell
        self.nbands_pc = pc_xml.nbands  
        self.nbands_sc = sc_xml.nbands 

        # This indicates the code to compute the projections of nbands - band_min bands
        self.band_min  = band_min 

        # Condition to control the number of computed bands  
        if self.band_min > self.nbands_sc:
           raise Exception("Minimum of bands larger than total number of bands")

        # Reciprocal lattice of the primitive cell and supercell in cartesian coordiantes
        self.rcell_pc = array(pc_xml.rcell)/pc_xml.cell[0][0] 
        self.rcell_sc = array(sc_xml.rcell)/sc_xml.cell[0][0]
   
        # Eigenvalues of the primitive cell and supercell
        self.eigen_pc = array(pc_xml.eigen1) 
        self.eigen_sc = array(sc_xml.eigen1) 

        # Format to save Miller indices
        format_string = "%12.4lf %12.4lf %12.4lf" 
        n_decs = 8 

        # Distance of the k-path
        kpoints_dists = calculate_distances(self.kpoints) 
       
        # Condition to use a rotation matrix if needed 
        if sc_rotated is False:
            self.rot = np.array([ [1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0] ])

        elif isinstance(sc_rotated, np.ndarray) and sc_rotated.shape == (3,3):
            self.rot = sc_rotated

        else:
            raise Exception("Input for sc_rotated must be False or a 3x3 matrix")

        # Array to save the projections
        self.projection = zeros([self.nkpoints_sc,self.nbands_sc-self.band_min]) 

        print("=== Data initialized successfully ===")

        print("=== Calculation of the projections  ===")

        # Condition to compute the projections or to use the saved ones  
        if compute_projections == True:

           # Loop along the kpoints of the supercell  
           for ik in range(self.nkpoints_sc):
               load(ik,self.nkpoints_sc)

               # Loading the data of each kpoint of the primitive cell and supercell 
               f_pc = h5py.File('%s/%s.save/wfc%01d.hdf5' % (self.path_pc,self.prefix_pc,(ik + 1)), 'r') 
               f_sc = h5py.File('%s/%s.save/wfc%01d.hdf5' % (self.path_sc,self.prefix_sc,(ik + 1)), 'r') 

               # Dimension of the Miller indices
               self.ng_pc = int(f_pc.attrs['igwx'])
               self.ng_sc = int(f_sc.attrs['igwx']) 

               # Creating a dictionary to collect key-value pairs
               g_sc = dict()
               g_sc_int = dict()
           
               # Loading Miller indices
               mill_sc = f_sc['MillerIndices']
               mill_pc = f_pc['MillerIndices']
        
               # Loop along the number of Miller indices of the super cell
               for ig in arange(self.ng_sc):
              
                   # Definition of Miller indices in order to reconstruct the k-points (to rotate it if neccessary) 
                   # and to save them in the dictionary associated to an integer
                   h,k,l = mill_sc[ig,0], mill_sc[ig,1], mill_sc[ig,2]
                   g_sc_int[(int(h),int(k),int(l))] = ig
                   w = h*self.rcell_sc[:][0] + k*self.rcell_sc[:][1] + l*self.rcell_sc[:][2]
                   w = dot(self.rot,w)
                   w = np.around(w, decimals = n_decs) + array([0,0,0])
                   w = format_string % (abs(w[0]), abs(w[1]), abs(w[2]))
                   g_sc[w] = ig
                   
               # To check in the following loop if the k-point of the super cell is found in the primitive cell    
               g_contain = [0]*self.ng_pc

               # Loop along the number of Miller indices of the primitive cell
               for ig in arange(self.ng_pc):

                   # Definition of Miller indices in order to reconstruct the k-points
                   h,k,l = mill_pc[ig,0], mill_pc[ig,1], mill_pc[ig,2]
                   w = h*self.rcell_pc[:][0] + k*self.rcell_pc[:][1] + l*self.rcell_pc[:][2]
                   w = np.around(w, decimals = n_decs) + array([0,0,0])
                   w = format_string % (abs(w[0]), abs(w[1]), abs(w[2]))

                   # Checking if the k-point in the supercell is found in the primitive cell,
                   # if missing, the projection will be wrong.
                   try:
                       g_contain[ig] = g_sc[w]
                   except KeyError:
                       print("Missing k-point %d" % ig)
           

               # Condition to read the eigenvectors for nspin = 1 or nspin = 4 in QE (To implement spin == col, i.e., nspin = 2 in QE)
               if spin == "none" or spin == "noncol":
              
                  # Loading eigenvectors of the super cell
                  evc_sc = f_sc['evc']
                  eivecs = []

                  # Loop along the number of bands indicated by the user to be computed
                  for ib in range(self.band_min,self.nbands_sc):
                      eivec = evc_sc[ib, :]
                      
                      # Rewriting the eigenvectors to manage them properly
                      eivec_complex = [complex(eivec[i], eivec[i+1]) for i in range(0, len(eivec), 2)]

                      eivecs.append(eivec_complex)
                        
                      # Defining specifically ib == 0 since it presents one component instead of two as the rest of values
                      if ib==0:
                         x = 0.0
                         for ig in range(self.ng_sc):
                             x += eivecs[-1][ig]*eivecs[-1][ig].conjugate()
 
               # Condition to compute the projections for nspin = 1 or spin = 4 in QE (To implement spin == col, i.e., nspin = 2 in QE)
               if spin == "none" or spin == "noncol":

                  # Loop along the number of bands indicated by the user to be computed
                  for ib in range(self.nbands_sc-self.band_min): 
                      x = 0.0
                      for ig in range(self.ng_pc): 

                          # Computing the projection between the primitive cell and the super cell
                          x += eivecs[ib][g_contain[ig]]*(eivecs[ib][g_contain[ig]].conjugate())

                      # If the value is less than a threshold, the projection is set to zero (to avoid ficticious points when plotting)
                      if abs(x) < 1e-4:

                         self.projection[ik,ib] = 0.0
                     
                      else:

                         self.projection[ik,ib] = abs(x)
               
               # Saving the data to avoid recomputing the projections twice
               np.save('projections',self.projection)


        # Loading the projections when already computed
        if compute_projections == False:

           print(" * Loading projections ...  ")

           try:

              self.projection = np.load('projections.npy')

           except FileNotFoundError:

              print(f"Error: the file projections.npy does not exits")
              raise

        return print("=== Projections calculated successfully ===")

    def plot_eigen_ax(self,ax,path=[],xlim=(),ylim=()):
        """
        Provisional plot function for quick visualization of the data. Useful for small calculations. 
        For large calculations is more useful loading the .npy and plot it with another script.
        """

        print(" * Plotting EBS ...  ")

        # Getting the data of the defined path
        if path:
            if isinstance(path,Path):
                path = path.get_indexes()
            ax.set_xticks( *list(zip(*path)) )
        ax.set_ylabel('E (eV)')

        # Computing the distance of the path 
        kpoints_dists = calculate_distances(self.kpoints)

        # Defining labels
        ticks, labels = list(zip(*path))
        ax.set_xticks([kpoints_dists[t] for t in ticks])
        ax.set_xticklabels(labels)
        ax.set_ylabel('E (eV)')

        # Plotting high-symmetry vertical lines
        for t in ticks:
            ax.axvline(kpoints_dists[t],c='k',lw=2)
        ax.axhline(0,c='k',lw=1)

        # Plotting the band for the primitive cell
        for ib in range(self.nbands_pc):
           ax.plot(kpoints_dists,self.eigen_pc[:,ib],'k--',lw=0.5)

        # Plotting the projection of the super cell 
        for ib in range(self.nbands_sc-self.band_min):
           #ax.plot(kpoints_dists,self.eigen_sc[:,ib],'darkgreen',lw=0.2)
           ax.scatter(kpoints_dists,self.eigen_sc[:,ib+self.band_min],s=self.projection[:,ib]*5,color='navy',edgecolor = None)

        # Establishing x and y axis limits 
        if xlim: ax.set_xlim(xlim)
        if ylim: ax.set_ylim(ylim)

        return print("=== EBS calculated successfully  ===")

# Load bar
def load(x,n):
    bar_length = 100
    x+=1
    ratio = x/float(n)
    c = int(ratio * bar_length)
    stdout.write(" * Computing projections of the supercell onto the primitive cell ["+"="*c+" "*(bar_length-c)+"] %03.3f%%" % (ratio*100))
    if (x==n): stdout.write("\n")
    stdout.flush()
    stdout.write("\r")

