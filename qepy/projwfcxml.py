# Copyright (C) 2015  Alejandro Molina-Sanchez, Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
import xml.etree.ElementTree as ET
from   qepy.auxiliary import *
from   matplotlib     import pyplot as plt
from   numpy import array, zeros

HatoeV = 27.2107

class projwfcXML():
    """ Class to read data from a Quantum espresso XML file
    """
    _proj_file = 'atomic_proj.xml'

    def __init__(self,prefix,path='.'):
        """ Initialize the structure with the path where the atomic_proj.xml is
        """
        self.prefix   = prefix
        self.path     = path
        self.datafile_xml = ET.parse( "%s/%s.save/%s"%(path, prefix, self._proj_file)).getroot()
        #get nkpoints
        self.nkpoints = int(self.datafile_xml.findall("HEADER/NUMBER_OF_K-POINTS")[0].text.strip())
        # Read the number of BANDS
        self.nbands   = int(self.datafile_xml.find("HEADER/NUMBER_OF_BANDS").text)
        #get fermi
        self.fermi    = float(self.datafile_xml.find("HEADER/FERMI_ENERGY").text)
        #get number of projections
        self.nproj    = int(self.datafile_xml.find("HEADER/NUMBER_OF_ATOMIC_WFC").text)
 
        self.eigen = None 
        self.proj  = None 

    def __str__(self):
        s = ""
        s += "nkpoints: %d\n"%self.nkpoints
        s += "nbands: %d\n"%self.nbands
        return s

    def plot_eigen(self, path=[], selected_orbitals=[]):
        """ Plot the band structure. The size of the points is the weigth of
            the selected orbitals.
            Under development to include also colormap and a dictionary for the
            selection of the orbitals...
        """
        renormalization = 25.0
###################################################################################
# Font selection and borders
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',serif="Computer Modern Roman",size=20)
#        rcParams['axes.linewidth']    = 3
#        rcParams['xtick.major.width'] = 2
#        rcParams['ytick.major.width'] = 2
###################################################################################

        if self.eigen is None:
          self.get_eigen()
          eigen = array(self.eigen)
        if self.proj  is None:
          self.get_proj()
      
        if path:
            if isinstance(path,Path):
                path = path.get_indexes()
            plt.xticks( *zip(*path) )
        plt.ylabel('E (eV)')

        #plot vertical line
        for point in path:
            x, label = point
            plt.axvline(x)
        plt.axhline(0)

        # Selection of the bands
        w_pro          = zeros([self.nkpoints,self.nbands])
        for ik in range(self.nkpoints):
          for ib in range(self.nbands):
            w_pro[ik,ib] = sum(abs(self.proj[ik,selected_orbitals,ib])**2)

        #plot bands
        for ib in range(self.nbands):
           plt.scatter(range(self.nkpoints),eigen[:,ib]*HatoeV - self.fermi*HatoeV,s=w_pro[:,ib]*renormalization,c='r',edgecolors='none')
        plt.show()

    def get_eigen(self):
        """ Return eigenvalues
        """
        datafile_xml = self.datafile_xml
        eigen = []
        for ik in xrange(self.nkpoints):
          eigen.append( map(float, self.datafile_xml.find("EIGENVALUES/K-POINT.%d/EIG"%(ik+1)).text.split() ))
        self.eigen = eigen
        return eigen

    def get_proj(self):
        """ Return projections  
        """
        datafile_xml = self.datafile_xml
        proj = zeros([self.nkpoints,self.nproj,self.nbands],dtype=complex) 
        for ik in range(self.nkpoints):
          for ip in range(self.nproj):
             projlist    = self.datafile_xml.find("PROJECTIONS/K-POINT.%d/ATMWFC.%d" % (ik+1,ip+1) ).text.splitlines()[1:-1]
             proj[ik,ip] = [ (lambda x,y: complex(float(x),float(y)))(*c.split(',')) for c in projlist ]
        self.proj = proj
        return proj      
