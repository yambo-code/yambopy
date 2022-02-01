# Copyright (C) 2018 Alejandro Molina-Sanchez, Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
from __future__ import print_function, division
import re
import xml.etree.ElementTree as ET
from numpy import array, zeros
from .lattice import Path, calculate_distances 
from .auxiliary import *

RytoeV = 13.605698066

class ProjwfcXML(object):
    """
    Class to read data from a Quantum espresso projwfc XML file.
    
    This file contains the projection of the Kohn-Sham stated in
    the atomic orbitals read from the pseudopotential
    """
    _proj_file = 'atomic_proj.xml'

    def __init__(self,prefix,output_filename='projwfc.log',path='.',qe_version='6.1'):
        """
        Initialize the structure with the path where the atomic_proj.xml is
        """
        self.qe_version   = qe_version
        self.prefix       = prefix
        self.path         = path
        self.datafile_xml = ET.parse( "%s/%s.save/%s"%(path, prefix, self._proj_file)).getroot()
        print('Running projwfcxml for QE version %s' % qe_version)

        if self.qe_version=='6.7' or self.qe_version=='7.0':
           self.nbands = int( self.datafile_xml.findall("HEADER")[0].attrib['NUMBER_OF_BANDS'] )
           self.nkpoints = int( self.datafile_xml.findall("HEADER")[0].attrib['NUMBER_OF_K-POINTS'] )
           self.spin_components = int( self.datafile_xml.findall("HEADER")[0].attrib['NUMBER_OF_SPIN_COMPONENTS'] )
           self.nproj = int( self.datafile_xml.findall("HEADER")[0].attrib['NUMBER_OF_ATOMIC_WFC'] )
           self.fermi = float( self.datafile_xml.findall("HEADER")[0].attrib['FERMI_ENERGY'] )*RytoeV
           #self.number_electrons = int( self.datafile_xml.findall("HEADER")[0].attrib['NUMBER_OF_ELECTRONS'])
        else :
           #get nkpoints
           self.nkpoints = int(self.datafile_xml.findall("HEADER/NUMBER_OF_K-POINTS")[0].text.strip())
           # Read the number of BANDS
           self.nbands   = int(self.datafile_xml.find("HEADER/NUMBER_OF_BANDS").text)
           #get fermi
           self.fermi    = float(self.datafile_xml.find("HEADER/FERMI_ENERGY").text)*RytoeV
           #get number of projections
           self.nproj    = int(self.datafile_xml.find("HEADER/NUMBER_OF_ATOMIC_WFC").text)
           #get weights of kpoints projections
           self.weights  = list(map(float,self.datafile_xml.find("WEIGHT_OF_K-POINTS").text.split()))
           #get number of spin components
           self.spin_components = int(self.datafile_xml.find("HEADER/NUMBER_OF_SPIN_COMPONENTS").text)
           #get kpoints

        # What is the utility of all this?
        #kpoints_lines = self.datafile_xml.find("K-POINTS").text.strip().split('\n')
        #kpoints_float = [ list(map(float, kline.split())) for kline in kpoints_lines ]
        #self.kpoints  = np.array(kpoints_float)

        self.kpoints = self.get_kpoints()

        # Read Eigenvalues

        if self.spin_components == 1: self.eigen = self.get_eigen()
        if self.spin_components == 2: self.eigen1,self.eigen2  = self.get_eigen()
        if self.spin_components == 4: self.eigen = self.get_eigen()

        # Read Atomic Oribtals Projections
          
        if self.spin_components == 1: self.proj  = self.get_proj()
        if self.spin_components == 2: self.proj1,self.proj2 = self.get_proj()
        if self.spin_components == 4: self.proj  = self.get_proj()

        #here we open the ouput file of projwfc and get the quantum numbers of the orbitals
        try:
            f = open("%s/%s"%(path,output_filename),'r')
        except:
            print("The output file of projwfc.x: %s was not found"%output_filename)
            exit(1)

        states = []
        #                                                                                        wfc                  j                 l                 m                    m_j
        for line in re.findall('state\s+\#\s+([0-9]+):\s+atom\s+([0-9]+)\s+\(([a-zA-Z]+)\s+\),\s+wfc\s+([0-9])\s+\((?:j=([0-9.]+))? ?(?:l=([0-9.]+))? ?(?:m=\s+([0-9.]+))? ?(?:m_j=([ \-0-9.]+))?',f.read()):
            # examples of the lines we have to read
            #  5: atom   1 (C  ), wfc  3 (l=2 m= 1)               #no spin case
            #  5: atom   1 (C  ), wfc  3 (j=1.5 l=1 m_j=-1.5)     #non collinear spin case
            _, iatom, atype, wfc, j, l, m, m_j = line
            state = {'iatom':int(iatom), 'atype':atype, 'wfc':int(wfc)}
            if j: j = float(j)
            if l: l = int(l)
            if m: m = int(m)
            if m_j: m_j = float(m_j)
            states.append({'iatom':int(iatom), 'atype':atype, 'wfc':int(wfc), 'j':j, 'l':l, 'm':m, 'm_j':m_j})
        self.states = states

        f.close()

    def __str__(self):
        s = ""
        s += "nkpoints: %d\n"%self.nkpoints
        s += "nbands: %d\n"%self.nbands
        return s

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

    def plot_eigen(self, ax, size=20, cmap=None, cmap2=None,color='r', color_2='b',path_kpoints=[], label_1=None, label_2=None,
                   selected_orbitals=[], selected_orbitals_2=[],bandmin=0,bandmax=None,alpha=1,size_projection=False,y_offset=0.0):
        """ 
        Plot the band structure. The size of the points is the weigth of the selected orbitals.

        Options:

            (a) Relative weight between two compositions. Pass a second set of orbitals
            (b) Colormap enters as a string
            (c) spin = 1 (no spin) and 2 (collinear spin)

        Under development to include also colormap and a dictionary for the
        selection of the orbitals...
        """
        from numpy import arange
        # Careful with the path variable! I am changing this variable to path_kpoints
        # Check we are not breaking the code some where
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        if path_kpoints:
            if isinstance(path_kpoints,Path):
                path_kpoints = path_kpoints.get_indexes()
        if bandmax is None or bandmax > self.nbands:
            bandmax = self.nbands

        #Colormap
        if cmap:
          color_map = plt.get_cmap(cmap)
        else:
          color_map = plt.get_cmap('rainbow')
        if cmap2:
          color_map2 = plt.get_cmap(cmap2)
        else:
          color_map2 = plt.get_cmap('rainbow')

        # Fix here
        #get kpoint_dists
        print(self.kpoints)
        kpoints_dists = calculate_distances(self.kpoints[:self.nkpoints])
 
        #make labels
        ticks, labels = list(zip(*path_kpoints))
        ax.set_xticks([kpoints_dists[t] for t in ticks])
        ax.set_xticklabels(labels)
        ax.set_ylabel('E (eV)')

        #plot vertical lines
        for t in ticks:
            ax.axvline(kpoints_dists[t],c='k',lw=2)
        ax.axhline(0,c='k')
     
        if selected_orbitals_2:
           # No spin or full spinor
           if self.spin_components == 1 or self.spin_components == 4:
              #get weights of second set of orbitals
              w_rel = self.get_relative_weight(selected_orbitals=selected_orbitals, selected_orbitals_2=selected_orbitals_2)
              #plot bands for fixed size
              for ib in range(bandmin,bandmax):
                  eig = self.eigen[:,ib] - self.fermi + y_offset
                  if size_projection==True:
                     cax = ax.scatter(kpoints_dists,eig,s=size[:,ib],c=w_rel[:,ib],cmap=color_map,vmin=0,vmax=1,edgecolors='none',label=label_1,rasterized=True,zorder=2)
                  else:
                     cax = ax.scatter(kpoints_dists,eig,s=size,c=w_rel[:,ib],cmap=color_map,vmin=0,vmax=1,edgecolors='none',label=label_1,rasterized=True,zorder=2)
                     #plt.plot(kpoints_dists,eig,'r-')#,s=size,c=w_rel[:,ib],cmap=color_map,vmin=0,vmax=1,edgecolors='none',label=label_1,rasterized=True,zorder=2)

           # Spin polarized
           if self.spin_components == 2:
              w_rel1, w_rel2 = self.get_relative_weight(selected_orbitals=selected_orbitals, selected_orbitals_2=selected_orbitals_2)
              #plot bands for fixed size
              for ib in range(bandmin,bandmax):
                  eig1 = self.eigen1[:,ib] - self.fermi + y_offset
                  eig2 = self.eigen2[:,ib] - self.fermi + y_offset
                  if size_projection==True:
                     cax = ax.scatter(kpoints_dists,eig,s=size[:,ib],c=w_rel[:,ib],cmap=color_map,vmin=0,vmax=1,edgecolors='none',label=label_1,rasterized=True,zorder=2)
                  else:
                    #cax = ax.scatter(kpoints_dists,eig2,s=size,c=w_rel2[:,ib],cmap=color_map,vmin=0,vmax=1,edgecolors='none',label=label_1,rasterized=True,zorder=2)
                    print(len(kpoints_dists))
                    print(len(eig1))

                    cax = ax.scatter(kpoints_dists,eig1,s=size,c='r',label=label_1,rasterized=True,zorder=2)
                    cax = ax.scatter(kpoints_dists,eig2,s=size,c='b',label=label_1,rasterized=True,zorder=2)

#          if self.spin_components == 2:
#             #get weights of second set of orbitals
#             w_rel1, w_rel2 = self.get_relative_weight(selected_orbitals=selected_orbitals, selected_orbitals_2=selected_orbitals_2)
#             #plot bands for fixed size
#             for ib in range(bandmin,bandmax):
#                 eig1 = self.eigen1[:,ib] - self.fermi
#                 eig2 = self.eigen2[:,ib] - self.fermi
#                 if size_projection==True:
#                    cax = ax.scatter(kpoints_dists,eig1,s=size[:,ib],c=w_rel1[:,ib],cmap=color_map,vmin=0,vmax=1,edgecolors='none',label=label_1)
#                    cax = ax.scatter(kpoints_dists,eig2,s=size[:,ib],c=w_rel2[:,ib],cmap=color_map2,vmin=0,vmax=1,edgecolors='none',label=label_2)
#                 else:
#                    cax = ax.scatter(kpoints_dists,eig1,s=size,c=w_rel1[:,ib],cmap=color_map,vmin=0,vmax=1,edgecolors='none',label=label_1)
#                    cax = ax.scatter(kpoints_dists,eig2,s=size,c=w_rel2[:,ib],cmap=color_map2,vmin=0,vmax=1,edgecolors='none',label=label_2)

        else:
#          if self.spin_components == 1:
          #plot bands for a varying size
             print('spin non-polarized')
             w_proj = self.get_weights(selected_orbitals=selected_orbitals)
             for ib in range(bandmin,bandmax):
                 eig = self.eigen[:,ib] - self.fermi + y_offset
                 cax = ax.scatter(kpoints_dists,eig,s=w_proj[:,ib]*size,c=color,edgecolors='none',alpha=alpha,label=label_1,rasterized=True,zorder=2)
                 #cax = ax.scatter(kpoints_dists,eig,s=1.0,c=color,edgecolors='none',alpha=alpha,label=label_1,rasterized=True,zorder=2)
#
#          if self.spin_components == 2:
#          #plot bands for a varying size
#             w_proj1, w_proj2 = self.get_weights(selected_orbitals=selected_orbitals)
#             for ib in range(bandmin,bandmax):
#                 eig1 = self.eigen1[:,ib] - self.fermi
#                 eig2 = self.eigen2[:,ib] - self.fermi
#                 cax = ax.scatter(kpoints_dists,eig1,s=w_proj1[:,ib]*size,c=color  ,edgecolors='none',alpha=alpha,label=label_1)
#                 cax = ax.scatter(kpoints_dists,eig2,s=w_proj2[:,ib]*size,c=color_2,edgecolors='none',alpha=alpha,label=label_2)

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

    def get_relative_weight(self,selected_orbitals=[],selected_orbitals_2=[],bandmin=0,bandmax=None):
        if bandmax is None:
            bandmax = self.nbands

        # No spin polarized
        if self.spin_components == 1 or self.spin_components == 4:
           # Selection of the bands
           w_rel = zeros([self.nkpoints,self.nbands])
           for ik in range(self.nkpoints):
               for ib in range(bandmin,bandmax):
                   a = sum(abs(self.proj[ik,selected_orbitals,ib])**2)
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
                   a1 = sum(abs(self.proj1[ik,selected_orbitals,ib])**2)
                   b1 = sum(abs(self.proj1[ik,selected_orbitals_2,ib])**2)
                   w_rel1[ik,ib] = a1/(a1+b1)
                   a2 = sum(abs(self.proj2[ik,selected_orbitals,ib])**2)
                   b2 = sum(abs(self.proj2[ik,selected_orbitals_2,ib])**2)
                   w_rel2[ik,ib] = a2/(a2+b2)
           return w_rel1, w_rel2

    def get_kpoints(self):

        if self.qe_version == '6.1':
           kpoints = []
           k_aux =  list( self.datafile_xml.find("K-POINTS").text.split() )
           for ik in range(self.nkpoints):
               kpoints.append([float(k_aux[ik*3]),float(k_aux[ik*3+1]),float(k_aux[ik*3+2])])

        elif self.qe_version == '6.7' or self.qe_version=='7.0':
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

           if self.qe_version == '6.7' or self.qe_version=='7.0':
              eigen =  [ list( map(float, word.text.split())) for word in self.datafile_xml.findall("EIGENSTATES/E") ] 

              self.eigen = np.array(eigen)*RytoeV
              return self.eigen

           if self.qe_version == '6.1':

              for ik in range(self.nkpoints):
                  eigen.append( list(map(float, self.datafile_xml.find("EIGENVALUES/K-POINT.%d/EIG"%(ik+1)).text.split())))  # version before 6.7
              self.eigen = np.array(eigen)*RytoeV
              return self.eigen
               #exit()
               #eigen.append( list(map(float, self.datafile_xml.find("EIGENSTATES/E"%(ik+1)).text.split())))  # version 6.7
               #exit()
           #exit()

        # Spin polarized
        if self.spin_components == 2:
           
            if self.qe_version == '6.7' or self.qe_version=='7.0':
               eigen_prov =  [ list( map(float, word.text.split())) for word in self.datafile_xml.findall("EIGENSTATES/E") ] 
               eigen_aux = np.array(eigen_prov)*RytoeV 
               self.eigen1 = eigen_aux[            0:  self.nkpoints,:]
               self.eigen2 = eigen_aux[self.nkpoints:2*self.nkpoints,:]
               return self.eigen1, self.eigen2

            if self.qe_version == '6.1':
               for ik in range(self.nkpoints):
                   eigen1.append( list(map(float, self.datafile_xml.find("EIGENVALUES/K-POINT.%d/EIG.1"%(ik+1)).text.split() )))
                   eigen2.append( list(map(float, self.datafile_xml.find("EIGENVALUES/K-POINT.%d/EIG.2"%(ik+1)).text.split() )))
               self.eigen1 = np.array(eigen1)*RytoeV
               self.eigen2 = np.array(eigen2)*RytoeV

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
        proj  = zeros([self.nkpoints,self.nproj,self.nbands],dtype=complex)

        if self.spin_components == 1 or self.spin_components == 4:

           # version 6.1
           if self.qe_version == '6.1':
              for ik in range(self.nkpoints):
                  for ip in range(self.nproj):
                      projlist    = self.datafile_xml.find("PROJECTIONS/K-POINT.%d/ATMWFC.%d" % (ik+1,ip+1) ).text.splitlines()[1:-1]
                      proj[ik,ip] = [ (lambda x,y: complex(float(x),float(y)))(*c.split(',')) for c in projlist ]

              self.proj = np.array(proj)

              return proj

           # version 6.7
           elif self.qe_version == '6.7' or self.qe_version=='7.0':
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

           # Spin polarized
           if self.spin_components == 2:

              return proj

        if self.spin_components == 2:
           
            if self.qe_version == '6.1':
               proj1 = zeros([self.nkpoints,self.nproj,self.nbands],dtype=complex)
               proj2 = zeros([self.nkpoints,self.nproj,self.nbands],dtype=complex)

               for ik in range(self.nkpoints):
                   for ip in range(self.nproj):
                       projlist1    = self.datafile_xml.find("PROJECTIONS/K-POINT.%d/SPIN.1/ATMWFC.%d" % (ik+1,ip+1) ).text.splitlines()[1:-1]
                       projlist2    = self.datafile_xml.find("PROJECTIONS/K-POINT.%d/SPIN.2/ATMWFC.%d" % (ik+1,ip+1) ).text.splitlines()[1:-1]
                       proj1[ik,ip] = [ (lambda x,y: complex(float(x),float(y)))(*c.split(',')) for c in projlist1 ]
                       proj2[ik,ip] = [ (lambda x,y: complex(float(x),float(y)))(*c.split(',')) for c in projlist2 ]

               self.proj1 = np.array(proj1)
               self.proj2 = np.array(proj2)

               return proj1, proj2

            # Two independent spinors
            elif self.qe_version == '6.7' or self.qe_version=='7.0':
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
    
    #def get_overlaps(self):


    def __str__(self):
        s  = "nbands:   %d\n"%self.nbands
        s += "nkpoints: %d\n"%self.nkpoints
        for n,state in enumerate(self.states):
            s += "n: %3d -> iatom:%3d atype:%2s wfc:%d j:%s l:%s m:%s m_j:%s\n"%(n,state['iatom'],state['atype'],state['wfc'],str(state['j']),str(state['l']),str(state['m']),str(state['m_j']))
        return s
