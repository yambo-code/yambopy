#
# License-Identifier: GPL
#
# Copyright (C) 2024 The Yambo Team
#
# Authors: HPC, AMS, FP, RR
#
# This file is part of the yambopy project
#
import re
import xml.etree.ElementTree as ET
from numpy import array, zeros, pi, conjugate, arange
from .lattice import Path, calculate_distances 
from yambopy.tools.string import marquee
from .auxiliary import *
from itertools import chain

RytoeV = 13.605698066
HatoeV = 2.0*RytoeV

class ProjwfcXML(object):
    """
    Class to read data from a Quantum espresso projwfc XML file.
    
    This file contains the projection of the Kohn-Sham stated in
    the atomic orbitals read from the pseudopotential
    """
    _proj_file = 'atomic_proj.xml'

    def __init__(self,prefix,output_filename='projwfc.log',path='.',qe_version='7.0'):
        """
        Initialize the structure with the path where the atomic_proj.xml is
        The zero energy is fixed at the Fermi energy
        Check conversion of Fermi energy (Ry to eV or Ha to eV)
        """
        self.qe_version   = qe_version
        self.prefix       = prefix
        self.path         = path
        self.get_qe_version(path,output_filename)
        version_number = float(re.findall(r"([0-9.]+)",self.qe_version)[0])
        if version_number<6.7: raise NotImplementedError(f"Warning: QE version {self.qe_version} not supported by this script (needs to be>=6.7).")
        if version_number>=7.0: self.order_is_l_j_mj   = True
        else:                   self.order_is_j_l_m_mj = True
        self.datafile_xml = ET.parse( "%s/%s.save/%s"%(path, prefix, self._proj_file)).getroot()

        # Read the number of BANDS
        self.nbands = int( self.datafile_xml.findall("HEADER")[0].attrib['NUMBER_OF_BANDS'] )
        #get nkpoints
        self.nkpoints = int( self.datafile_xml.findall("HEADER")[0].attrib['NUMBER_OF_K-POINTS'] )
        #get number of spin components
        self.spin_components = int( self.datafile_xml.findall("HEADER")[0].attrib['NUMBER_OF_SPIN_COMPONENTS'] )
        #get number of projections
        self.nproj = int( self.datafile_xml.findall("HEADER")[0].attrib['NUMBER_OF_ATOMIC_WFC'] )
        #get fermi
        self.fermi = float( self.datafile_xml.findall("HEADER")[0].attrib['FERMI_ENERGY'] )*RytoeV
        #self.number_electrons = int( self.datafile_xml.findall("HEADER")[0].attrib['NUMBER_OF_ELECTRONS'])
           
        #get kpoints
        self.kpoints = self.get_kpoints()

        # Read Eigenvalues

        if self.spin_components == 1: self.eigen = self.get_eigen()
        if self.spin_components == 2: self.eigen1,self.eigen2  = self.get_eigen()
        if self.spin_components == 4: self.eigen = self.get_eigen()

        # Read Atomic Orbitals Projections
          
        if self.spin_components == 1: self.proj  = self.get_proj()
        if self.spin_components == 2: self.proj1,self.proj2 = self.get_proj()
        if self.spin_components == 4: self.proj  = self.get_proj()

        #here we open the ouput file of projwfc and get the quantum numbers of the orbitals
        try:
            # Try to open the file in the specified path
            f = open(os.path.join(path, output_filename), 'r')
        except FileNotFoundError:
            try:
                # Try to open the file in the './out' subdirectory of the specified path
                f = open(os.path.join(path, './out', output_filename), 'r')
            except FileNotFoundError:
                # Raise an exception if the file is not found in both locations
                raise Exception(f"The output file of projwfc.x: {output_filename} was not found in the directory or its ./out subdirectory.")        
            
        if hasattr(self, 'order_is_l_j_mj') and self.order_is_l_j_mj:
            states = []
            #                                                                                         wfc                  l                 j                 m_j                 
            for line in re.findall('state\s+\#\s+([0-9]+):\s+atom\s+([0-9]+)\s+\(([a-zA-Z]+)\s+\),\s+wfc\s+([0-9])\s+\((?:l=([0-9.]+))? ?(?:j=([0-9.]+))? ?(?:m_j=\s+([0-9.]+))?',f.read()):
                # examples of the lines we have to read
                #  5: atom   1 (C  ), wfc  3 (l=2 m= 1)               #no spin case
                #  5: atom   1 (C  ), wfc  3 (j=1.5 l=1 m_j=-1.5)     #non collinear spin case
                istate, iatom, atype, wfc, l, j, m_j = line
                if j: j = float(j)
                if l: l = int(l)
                if m_j: m_j = float(m_j)
                states.append({'istate':int(istate), 'iatom':int(iatom), 'atype':atype, 'wfc':int(wfc), 'l':l, 'j':j, 'm_j':m_j})
            self.states = states

            f.close() 

        if hasattr(self, 'order_is_j_l_m_mj') and self.order_is_j_l_m_mj:      
            states = []
            #                                                                                        wfc                  j                 l                 m                    m_j
            for line in re.findall('state\s+\#\s+([0-9]+):\s+atom\s+([0-9]+)\s+\(([a-zA-Z]+)\s+\),\s+wfc\s+([0-9])\s+\((?:j=([0-9.]+))? ?(?:l=([0-9.]+))? ?(?:m=\s+([0-9.]+))? ?(?:m_j=([ \-0-9.]+))?',f.read()):
                # examples of the lines we have to read
                #  5: atom   1 (C  ), wfc  3 (l=2 m= 1)               #no spin case
                #  5: atom   1 (C  ), wfc  3 (j=1.5 l=1 m_j=-1.5)     #non collinear spin case
                istate, iatom, atype, wfc, j, l, m, m_j = line
                if j: j = float(j)
                if l: l = int(l)
                if m: m = int(m)
                if m_j: m_j = float(m_j)
                states.append({'istate':int(istate),'iatom':int(iatom), 'atype':atype, 'wfc':int(wfc), 'j':j, 'l':l, 'm':m, 'm_j':m_j})
            self.states = states

            f.close()
    
    def get_qe_version(self,path,output_filename):
        """
        Get the version of Quantum Espresso. v<6.7 is not supported.

        Try both data-file-schema and projwfc output
        """
        try:
            datafile_schema_xml = ET.parse( f"{path}/{self.prefix}.save/data-file-schema.xml").getroot()
            qe_version_xml = datafile_schema_xml.findall("general_info/creator")[0].attrib["VERSION"].strip()
        except FileNotFoundError: qe_version_xml = None
        try:
            with open(f"{path}/{output_filename}") as f: qe_version_output = re.findall(r"Program PROJWFC v.([0-9.]+[a-zA-Z]*)", f.read())[0]
        except FileNotFoundError: qe_version_output = None
        if qe_version_xml and qe_version_output: 
            if qe_version_xml != qe_version_output: print(f"Warning: QE version mismatch between data-file-schema.xml ({qe_version_xml}) and projwfc output ({qe_version_output}). Using the output version for compatibility.")
            self.qe_version = qe_version_output
        elif qe_version_xml:    self.qe_version = qe_version_xml
        elif qe_version_output: self.qe_version = qe_version_output
        else: print(f"Warning: QE version not found in data-file-schema.xml or projwfc output. Using default version {self.qe_version}.")

    def get_indexes(self):
        """
        Get indexes of the bands where the projection is maximal
        """
        # Selection of the bands
        proj = zeros([self.nkpoints,self.nproj],dtype=int)
        for ik in range(self.nkpoints):
            for ip in range(self.nproj):
                proj[ik,ip] = np.argmax(np.absolute(self.proj[ik,ip,:])**2)

        return proj

    def get_states_helper(self, atom_query=['all'], orbital_query=['s','p','d','f']):
        """
        Get the sates that you want based on dictionary query by providing array of atoms and orbitals, default all orbitals
        
        Returns an array with the indices of the requested states in the qe array
        """
        states =  self.states
        queried_states = []

        for state in states:
            if (atom_query == ['all']) or (state['atype'] in atom_query):
                if (state['l']==0) and 's' in orbital_query:         #s orbital
                    queried_states.append(state['istate'] - 1)
                if (state['l']==1) and 'p' in orbital_query:         #p orbital
                    queried_states.append(state['istate'] - 1)
                if (state['l']==2) and 'd' in orbital_query:         #d orbital
                    queried_states.append(state['istate'] - 1)
                if (state['l']==3) and 'f' in orbital_query:         #f orbital
                    queried_states.append(state['istate'] - 1)
                
        return queried_states


    def plot_eigen(self, ax, size=20, cmap=None, cmap2=None,color='r', color_2='b',path_kpoints=[], label_1='', label_2='',
                   selected_orbitals=[], selected_orbitals_2=[],bandmin=0,bandmax=None,alpha=1,size_projection=False,y_offset=0.0,marker='o'):
        """ 
        Plot the band structure. The size of the points is the weigth of the selected orbitals.

        Options:

            (a) Relative weight between two compositions: selected_orbitals and selected_orbitas_2
                Format >>> selected_orbitals = [0,2,4]
            (b) Colormap enters as a string

        Arguments for plot layout:
        - size:     size of markers in scatterplot        [default 20]
        - alpha:    alpha value of markers in scatterplot [default 1]
        - label_1, label_2: plot labels [default empty string]

        Under development to include also colormap and a dictionary for the
        selection of the orbitals...
        example usage to get: 
             state #   2: atom   1 (Li ), wfc  2 (l=1 m= 1)
        
        plot_eigen(ax, path_kpoints=path_kpoints, selected_orbitals=[1], color=color, size=dotsize) 

            notice python counting; state# - 1 = selected_orbital index
        """
        #from numpy import arange # redundent
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        if path_kpoints:
            if isinstance(path_kpoints,Path):
                path_kpoints = path_kpoints.get_indexes()
        if bandmax is None or bandmax > self.nbands:
            bandmax = self.nbands

        #Colormap
        if cmap:  color_map = plt.get_cmap(cmap)
        else:     color_map = plt.get_cmap('rainbow')
        if cmap2: color_map2 = plt.get_cmap(cmap2)
        else:     color_map2 = plt.get_cmap('rainbow')

        #get kpoint_dists
        kpoints_dists = calculate_distances(self.kpoints[:self.nkpoints])
 
        #make K-points labels
        ticks, labels = list(zip(*path_kpoints))
        ax.set_xticks([kpoints_dists[t] for t in ticks])
        ax.set_xticklabels(labels)
        ax.set_ylabel('E (eV)')

        #plot vertical lines
        for t in ticks: ax.axvline(kpoints_dists[t],c='k',lw=1)
        ax.axhline(0,c='k')
     
        # Plot bands for fixed size in a colormap
        if selected_orbitals_2:
           # No spin or full spinor
           if self.spin_components == 1 or self.spin_components == 4:
              w_rel = self.get_relative_weight(selected_orbitals=selected_orbitals, selected_orbitals_2=selected_orbitals_2)
              for ib in range(bandmin,bandmax):
                  eig = self.eigen[:,ib] + y_offset
                  eig_last = self.eigen[:,-1] + y_offset
                  state = self.states
                  j = state[ib]['j']
                  l = state[ib]['l']
                  if size_projection==True:
                     cax = ax.scatter(kpoints_dists,eig,s=size[:,ib],c=w_rel[:,ib],cmap=color_map,vmin=0,vmax=1,edgecolors='none',label=label_1,rasterized=True,zorder=2,marker=marker)
                  else:
                     cax = ax.scatter(kpoints_dists,eig,s=size,c=w_rel[:,ib],cmap=color_map,vmin=0,vmax=1,edgecolors='none',label=label_1,rasterized=True,zorder=2)
                     #print(f'b:{ib} J={j} L={l}')
                     #ax.textye(f'b:{ib} J={j} L={l}', ((kpoints_dists[-1]-kpoints_dists[0])/2,ib*(eig_last-eig)/(bandmax-bandmin)), textcoords='offset points', xytext=(0,10), ha='center', va='bottom',color='teal')
                     #ax.annotate(f'b:{ib} J={j} L={l}', ((kpoints_dists[-1]-kpoints_dists[0])/2,ib*(eig_last-eig)/(bandmax-bandmin)), textcoords='offset points', xytext=(0,10), ha='center', va='bottom',color='teal')

           # Spin polarized no SOC
           if self.spin_components == 2:
              w_rel1, w_rel2 = self.get_relative_weight(selected_orbitals=selected_orbitals, selected_orbitals_2=selected_orbitals_2)
              for ib in range(bandmin,bandmax):
                  eig1 = self.eigen1[:,ib] + y_offset
                  eig2 = self.eigen2[:,ib] + y_offset
                  if size_projection==True:
                     cax = ax.scatter(kpoints_dists,eig,s=size[:,ib],c=w_rel[:,ib],cmap=color_map,vmin=0,vmax=1,edgecolors='none',label=label_1,rasterized=True,zorder=2,marker=marker)
                  else:
                     cax = ax.scatter(kpoints_dists,eig1,s=size,c=w_rel1[:,ib],cmap=color_map,vmin=0,vmax=1,edgecolors='none',label=label_1,rasterized=True,zorder=2,marker=marker)
                     cax2= ax.scatter(kpoints_dists,eig2,s=size,c=w_rel2[:,ib],cmap=color_map2,vmin=0,vmax=1,edgecolors='none',label=label_1,rasterized=True,zorder=2,marker=marker)

        # Plot bands with changing size and a fixed color
        else:
            if self.spin_components == 1 or self.spin_components == 4:
               w_proj = self.get_weights(selected_orbitals=selected_orbitals)
               ib_max = np.where(w_proj==np.max(w_proj))[1][0] # needed for label
               for ib in range(bandmin,bandmax):
                   if ib==ib_max: lab = label_1
                   else:          lab = '_'+str(label_1)
                   eig = self.eigen[:,ib] + y_offset
                   cax = ax.scatter(kpoints_dists,eig,s=w_proj[:,ib]*size,c=color,edgecolors='none',alpha=alpha,label=lab,rasterized=True,zorder=2,marker=marker)

            elif self.spin_components == 2:
                 w_proj1, w_proj2 = self.get_weights(selected_orbitals=selected_orbitals)
                 ib_max1, ib_max2 = np.where(w_proj1==np.max(w_proj1))[1][0], np.where(w_proj2==np.max(w_proj2))[1][0]
                 for ib in range(bandmin,bandmax):
                     lab1, lab2 = ['_'+label_1,'_'+label_2]
                     if ib==ib_max1: lab1 = label_1
                     if ib==ib_max2: lab2 = label_2
                     eig1, eig2 = self.eigen1[:,ib], self.eigen2[:,ib]
                     cax = ax.scatter(kpoints_dists,eig1,s=w_proj1[:,ib]*size,c=color  ,edgecolors='none',alpha=alpha,label=lab1,marker=marker)
                     cax2= ax.scatter(kpoints_dists,eig2,s=w_proj2[:,ib]*size,c=color_2,edgecolors='none',alpha=alpha,label=lab2,marker=marker)

        ax.set_xlim(0, max(kpoints_dists))
        return cax
    
    def get_weights(self,selected_orbitals=[],bandmin=0,bandmax=None):
        if bandmax is None:
           bandmax = self.nbands

        if self.spin_components == 1:

           # Selection of the bands
           w_proj = zeros([self.nkpoints,self.nbands])
           for ik in range(self.nkpoints):
               for ib in range(bandmin,bandmax):
                   w_proj[ik,ib] = sum(abs(self.proj[ik,selected_orbitals,ib])**2)
           return w_proj

        if self.spin_components == 2:

           # Selection of the bands
           w_proj1 = zeros([self.nkpoints,self.nbands])
           w_proj2 = zeros([self.nkpoints,self.nbands])

           for ik in range(self.nkpoints):
               for ib in range(bandmin,bandmax):
                   w_proj1[ik,ib] = sum(abs(self.proj1[ik,selected_orbitals,ib])**2)
                   w_proj2[ik,ib] = sum(abs(self.proj2[ik,selected_orbitals,ib])**2)
           return w_proj1, w_proj2

    def get_dorbitals_projection(self,selected_orbitals=[],bandmin=0,bandmax=None):
        """
        This function return the weights for d-orbitals in the basis of a1g, e+
        and e- 
        selected_orbitals must the indices of the d-orbital in the array in the order of QE
        """
        if bandmax is None:
           bandmax = self.nbands

        #if self.spin_components == 1:

           # Selection of the bands
           #w_proj = zeros([self.nkpoints,self.nbands])
           #for ik in range(self.nkpoints):
           #    for ib in range(bandmin,bandmax):
           #        w_proj[ik,ib] = sum(abs(self.proj[ik,selected_orbitals,ib])**2)
           #return w_proj

        if self.spin_components == 2:

           # Selection of the bands
           w_proj1 = zeros([self.nkpoints,self.nbands])
           w_proj2 = zeros([self.nkpoints,self.nbands])

           for ik in range(self.nkpoints):
               for ib in range(bandmin,bandmax):

                   w_proj1[ik,ib] = sum(abs(self.proj1[ik,selected_orbitals,ib])**2)
                   w_proj2[ik,ib] = sum(abs(self.proj2[ik,selected_orbitals,ib])**2)
           return w_proj1, w_proj2

    def get_pdos(self,selected_orbitals=None,bandmin=0,bandmax=None,energy_steps=100,e_min=-10.0,e_max=5.0,Gamma=0.1):

        print(selected_orbitals)

        energy_grid = arange(e_min,e_max,(e_max-e_min)/energy_steps)
        if bandmax is None:
           bandmax = self.nbands
        if selected_orbitals is None: selected_orbitals = range(self.nproj)

        if self.spin_components == 2:
           pdos_up, pdos_dw = zeros([energy_steps,self.nproj]), zeros([energy_steps,self.nproj])
           for ik in range(self.nkpoints):
               for ib in range(bandmin,bandmax):
                   for io in selected_orbitals:
                       for ie,e in enumerate(energy_grid):
                           pdos_up[ie,io] = pdos_up[ie,io] + abs(conjugate(self.proj1[ik,io,ib])*self.proj1[ik,io,ib])*self._lorentz(self.eigen1[ik,ib],e,Gamma)

        return energy_grid, pdos_up #, dos_up       

    def _lorentz(self,x,x_0,Gamma):
        
        return (1.0/pi)*(0.5*Gamma)/((x-x_0)**2 + (0.5*Gamma)**2)
    
    def get_relative_weight(self,selected_orbitals=[],selected_orbitals_2=[],bandmin=0,bandmax=None):
        if bandmax is None:
            bandmax = self.nbands

        # No spin polarized
        if self.spin_components == 1 or self.spin_components == 4:
           # Selection of the bands
           w_rel = zeros([self.nkpoints,self.nbands])
           for ik in range(self.nkpoints):
               for ib in range(bandmin,bandmax):
                   a = sum(abs(self.proj[ik,selected_orbitals  ,ib])**2)
                   b = sum(abs(self.proj[ik,selected_orbitals_2,ib])**2)
                   w_rel[ik,ib] = a/(a+b)
           return w_rel

        # Spin polarized collinear
        if self.spin_components == 2:
           # Selection of the bands
           w_rel1 = zeros([self.nkpoints,self.nbands])
           w_rel2 = zeros([self.nkpoints,self.nbands])
           for ik in range(self.nkpoints):
               for ib in range(bandmin,bandmax):
                   a1 = sum(abs(self.proj1[ik,selected_orbitals  ,ib])**2)
                   b1 = sum(abs(self.proj1[ik,selected_orbitals_2,ib])**2)
                   w_rel1[ik,ib] = a1/(a1+b1)
                   a2 = sum(abs(self.proj2[ik,selected_orbitals  ,ib])**2)
                   b2 = sum(abs(self.proj2[ik,selected_orbitals_2,ib])**2)
                   w_rel2[ik,ib] = a2/(a2+b2)
           return w_rel1, w_rel2

    def get_kpoints(self):
        """ Return kpoints
        """
        kpoints = []
        datafile_xml = self.datafile_xml
        for word in self.datafile_xml.findall("EIGENSTATES/K-POINT"):
            kpoints.append( list( map(float, word.text.split()) )  )
        
        return kpoints

    def get_eigen(self):
        """ Return eigenvalues
        """
        datafile_xml = self.datafile_xml
        eigen = []
        eigen1 = []
        eigen2 = []

        # No spin polarized
        if self.spin_components == 1 or self.spin_components == 4:

            eigen =  [ list( map(float, word.text.split())) for word in self.datafile_xml.findall("EIGENSTATES/E") ] 

            self.eigen = np.array(eigen)*RytoeV - self.fermi
            return self.eigen

        # Spin polarized
        if self.spin_components == 2:
           
            eigen_prov =  [ list( map(float, word.text.split())) for word in self.datafile_xml.findall("EIGENSTATES/E") ] 
            eigen_aux = np.array(eigen_prov)*RytoeV 
            self.eigen1 = eigen_aux[            0:  self.nkpoints,:] - self.fermi
            self.eigen2 = eigen_aux[self.nkpoints:2*self.nkpoints,:] - self.fermi

            return self.eigen1, self.eigen2

    def write_proj(self,filename='proj'):
        """
        Write the projection array in a numpy file
        """
        np.savez(filename,proj=self.proj,weights=self.weights)
            
    def get_proj(self):
        """ Return projections
        """
        datafile_xml = self.datafile_xml

        if self.spin_components == 1 or self.spin_components == 4:
            proj  = zeros([self.nkpoints,self.nproj,self.nbands],dtype=complex)

            data_atomic_wfc = self.datafile_xml.findall("EIGENSTATES/PROJS/ATOMIC_WFC")
            for ik in range(self.nkpoints):
                for ip in range(self.nproj):
                    i_data = ik*self.nproj + ip
                    projlist = data_atomic_wfc[i_data].text.splitlines()[1:-1]
                    atom_aux = []
                    for c in projlist:
                        z = float(c.split()[0]) + 1.0j*float(c.split()[1])
                        atom_aux.append(z)
                    proj[ik,ip] = atom_aux
            
            self.proj = np.array(proj)

            return proj

        # Spin polarized
        elif self.spin_components == 2:
            print('return spin polarized projections')
           
            # Two independent spinors
            data_atomic_wfc = self.datafile_xml.findall("EIGENSTATES/PROJS/ATOMIC_WFC")
            proj1 = zeros([self.nkpoints,self.nproj,self.nbands],dtype=complex)
            proj2 = zeros([self.nkpoints,self.nproj,self.nbands],dtype=complex)
            
            for ik in range(self.nkpoints):
                for ip in range(self.nproj):
                    i_data1 = ik*self.nproj + ip
                    i_data2 = (ik+self.nkpoints)*self.nproj + ip
                    projlist1 = data_atomic_wfc[i_data1].text.splitlines()[1:-1]
                    projlist2 = data_atomic_wfc[i_data2].text.splitlines()[1:-1]
                    atom_aux1, atom_aux2 = [], []
                    for c in projlist1:
                        z = float(c.split()[0]) + 1.0j*float(c.split()[1])
                        atom_aux1.append(z)
                    for c in projlist2:
                        z = float(c.split()[0]) + 1.0j*float(c.split()[1])
                        atom_aux2.append(z)
                    proj1[ik,ip] = atom_aux1
                    proj2[ik,ip] = atom_aux2

            self.proj1 = np.array(proj1)
            self.proj2 = np.array(proj2)

            return proj1, proj2
       
    def shift_bands(self,qpcorrection,vb,cb):
        """
        Shift band structure, e.g. to account for a G0W0 run.

        The idea is to shift the bands by the qp corrections in order to use 
        the weights from a projwfc calculation.

        Note that the path used in QE must be the same as the one used in Yambo

        - vb and cb allow for tuning the bands that are to be corrected
        """

        for ib in  range(self.nbands-vb,self.bands,self.nbands+cb): 
            self.eigen[:,ib] = self.eigen[:,ib]+qpcorrection[:,ib]
 
    def add_scissor(self,n_val,scissor):
        """
        Apply a scissor operator to the bands, e.g. to account for a G0W0 run.
            - n_val: number of valence bands
            - scissor: list of three numbers [shift, cond_stretch, val_stretch]
            - bands: array of eigenvalues
        """
        def band_formats(bands):
            """ Yambo format has leading spin axis
            """
            if len(bands.shape)>2: return bands[0] #Yambo spin-polarized case not supported
            else:                  return bands

        scissored_bands = np.zeros(self.eigen.shape)

        tmp_bands    = band_formats(self.eigen)
        tmp_sc_bands = band_formats(scissored_bands)

        top_v, bottom_c = tmp_bands[:,n_val-1], tmp_bands[:,n_val]
        ind_k_dir_gap = np.argmin(bottom_c-top_v)
        ev_max, ec_min = top_v[ind_k_dir_gap], bottom_c[ind_k_dir_gap]

        for ik in range( self.nkpoints ):
            for ib in range( self.nbands ):
                if ib<n_val: tmp_sc_bands[ik][ib] = ev_max-(ev_max-tmp_bands[ik][ib])*scissor[2]
                else: tmp_sc_bands[ik][ib] = ec_min+scissor[0]+(tmp_bands[ik][ib]-ec_min)*scissor[1]
        self.eigen0 = self.eigen # save original bands
        self.eigen  = scissored_bands

    def __str__(self):
        lines = []; app = lines.append
        app(marquee(self.__class__.__name__))

        app(f"Running projwfcxml for QE version {self.qe_version}")
        app(f"nkpoints: {self.nkpoints}")
        app(f"nbands:   {self.nbands}")
        if (self.qe_version=='7.0'): #replace with order_l-j-m-mj and order_j-l-m-mj 
            for n,state in enumerate(self.states):
                app(f"n: {n} -> iatom:{state['iatom']} atype:{state['atype']} wfc:{state['wfc']} l:{state['l']} j:{state['j']} m_j:{state['m_j']}")
            return "\n".join(lines)
        else:
            for n,state in enumerate(self.states):
                app(f"n: {n} -> iatom:{state['iatom']} atype:{state['atype']} wfc:{state['wfc']} j:{state['j']} l:{state['l']} m:{state['m']} m_j:{state['m_j']}")
            return "\n".join(lines)
