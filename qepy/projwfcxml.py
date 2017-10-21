from __future__ import print_function, division
#
# Copyright (C) 2017  Alejandro Molina-Sanchez, Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
import xml.etree.ElementTree as ET
from   qepy.auxiliary import *
from   numpy import array, zeros
import re

RytoeV = 13.605698066

class ProjwfcXML(object):
    """
    Class to read data from a Quantum espresso projwfc XML file.
    
    This file contains the projection of the Kohn-Sham stated in
    the atomic orbitals read from the pseudopotential
    """
    _proj_file = 'atomic_proj.xml'

    def __init__(self,prefix,output_filename='projwfc.log',path='.'):
        """
        Initialize the structure with the path where the atomic_proj.xml is
        """
        self.prefix   = prefix
        self.path     = path
        self.datafile_xml = ET.parse( "%s/%s.save/%s"%(path, prefix, self._proj_file)).getroot()
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

        self.eigen = self.get_eigen()
        self.proj  = self.get_proj()

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

    def plot_eigen(self, ax, size=20, cmap=None, color='r', path=[], selected_orbitals=[], selected_orbitals_2=[]):
        """ 
        Plot the band structure. The size of the points is the weigth of the selected orbitals.

        Options:

            (a) Relative weight between two compositions. Pass a second set of orbitals
            (b) Colormap enters as a string

        Under development to include also colormap and a dictionary for the
        selection of the orbitals...
        """
        import matplotlib.pyplot as plt
        if path:
            if isinstance(path,Path):
                path = path.get_indexes()
        ticks, labels = list(zip(*path))
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels)
        ax.set_ylabel('E (eV)')

        #Colormap
        if cmap:
          color_map = plt.get_cmap(cmap)

        #plot vertical line
        for x, label in path:
            ax.axvline(x,c='k',lw=2)
        ax.axhline(0,c='k')

        #get weights
        w_proj = self.get_weights(selected_orbitals=selected_orbitals)
        
        #get weights of second set of orbitals
        if selected_orbitals_2:
          w_rel = self.get_relative_weight(selected_orbitals=selected_orbitals, selected_orbitals_2=selected_orbitals_2)
          #plot bands for fix size
          for ib in range(self.nbands):
            ax.scatter(list(range(self.nkpoints)),self.eigen[:,ib] - self.fermi,s=size,c=w_rel[:,ib],cmap=color_map,edgecolors='none')

        #plot bands for a varying size
        if not selected_orbitals_2:
          for ib in range(self.nbands):
            #ax.scatter(range(self.nkpoints),self.eigen[:,ib] - self.fermi,c='r',edgecolors='none')
            ax.scatter(list(range(self.nkpoints)),self.eigen[:,ib] - self.fermi,s=w_proj[:,ib]*size,c=color,edgecolors='none')

        ax.set_xlim(0, self.nkpoints-1)
        ax.set_ylim(auto=True)

    def get_weights(self,selected_orbitals=[]):
        # Selection of the bands
        w_proj = zeros([self.nkpoints,self.nbands])
        for ik in range(self.nkpoints):
          for ib in range(self.nbands):
            w_proj[ik,ib] = sum(abs(self.proj[ik,selected_orbitals,ib])**2)
        return w_proj

    def get_relative_weight(self,selected_orbitals=[],selected_orbitals_2=[]):
        # Selection of the bands
        w_rel = zeros([self.nkpoints,self.nbands])
        for ik in range(self.nkpoints):
          for ib in range(self.nbands):
            a = sum(abs(self.proj[ik,selected_orbitals,ib])**2)
            b = sum(abs(self.proj[ik,selected_orbitals,ib])**2)
            w_rel[ik,ib] = a/(a+b)
        return w_rel

    def get_eigen(self):
        """ Return eigenvalues
        """
        datafile_xml = self.datafile_xml
        eigen = []
        for ik in range(self.nkpoints):
          eigen.append( list(map(float, self.datafile_xml.find("EIGENVALUES/K-POINT.%d/EIG"%(ik+1)).text.split() )))
        self.eigen = np.array(eigen)*RytoeV
        return self.eigen

    def write_proj(self,filename='proj'):
        """
        Write the projection array in a numpy file
        """
        np.savez(filename,proj=self.proj,weights=self.weights)
        
    def get_proj(self):
        """ Return projections
        """
        datafile_xml = self.datafile_xml
        proj = zeros([self.nkpoints,self.nproj,self.nbands],dtype=complex)
        for ik in range(self.nkpoints):
          for ip in range(self.nproj):
             projlist    = self.datafile_xml.find("PROJECTIONS/K-POINT.%d/ATMWFC.%d" % (ik+1,ip+1) ).text.splitlines()[1:-1]
             proj[ik,ip] = [ (lambda x,y: complex(float(x),float(y)))(*c.split(',')) for c in projlist ]
        self.proj = np.array(proj)
        return proj

    def __str__(self):
        s  = "nbands:   %d\n"%self.nbands
        s += "nkpoints: %d\n"%self.nkpoints
        for n,state in enumerate(self.states):
            s += "n: %3d -> iatom:%3d atype:%2s wfc:%d j:%s l:%s m:%s m_j:%s\n"%(n,state['iatom'],state['atype'],state['wfc'],str(state['j']),str(state['l']),str(state['m']),str(state['m_j']))
        return s
