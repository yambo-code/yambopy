# Copyright (C) 2015  Alejandro Molina-Sanchez, Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
import xml.etree.ElementTree as ET
from   qepy.auxiliary import *
from   numpy import array, zeros
import re

RytoeV = 13.605698066

class ProjwfcXML():
    """ Class to read data from a Quantum espresso XML file
    """
    _proj_file = 'atomic_proj.xml'

    def __init__(self,prefix,output_filename='projwfc.log',path='.'):

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

        self.eigen = self.get_eigen()
        self.proj  = self.get_proj()

        #here we open the ouput file of projwfc and get the quantum numbers of the orbitals
        try:
            f = open("%s/%s"%(path,output_filename),'r')
        except:
            print "The output file of projwfc.x: %s was not found"%output_filename
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
        """ get indexes of the bands where the projection is maximal
        """
        # Selection of the bands
        proj = zeros([self.nkpoints,self.nproj],dtype=int)
        for ik in range(self.nkpoints):
            for ip in range(self.nproj):
                proj[ik,ip] = np.argmax(np.absolute(self.proj[ik,ip,:])**2)

        return proj

    def plot_eigen(self, path=[], selected_orbitals=[]):
        """ Plot the band structure. The size of the points is the weigth of
            the selected orbitals.
            Under development to include also colormap and a dictionary for the
            selection of the orbitals...
        """
        from   matplotlib     import pyplot as plt
        renormalization = 25.0
        # Font selection and borders
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',serif="Computer Modern Roman",size=20)
        #rcParams['axes.linewidth']    = 3
        #rcParams['xtick.major.width'] = 2
        #rcParams['ytick.major.width'] = 2

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
        w_proj = zeros([self.nkpoints,self.nbands])
        for ik in range(self.nkpoints):
          for ib in range(self.nbands):
            w_proj[ik,ib] = sum(abs(self.proj[ik,selected_orbitals,ib])**2)

        #plot bands with linewithds
        #for ib in range(self.nbands):
        #    x = range(self.nkpoints)
        #    y = (self.eigen[:,ib] - self.fermi)*RytoeV
        #    points = np.array([x, y]).T.reshape(-1, 1, 2)
        #    segments = np.concatenate([points[:-1], points[1:]], axis=1)
        #    print w_proj[:,ib]
        #    lc = LineCollection(segments, linewidths=1+w_proj[:,ib]*10)
        #    plt.gca().add_collection(lc)

        #plot bands
        for ib in range(self.nbands):
           plt.scatter(range(self.nkpoints),self.eigen[:,ib]*RytoeV - self.fermi*RytoeV,s=w_proj[:,ib]*renormalization,c='r',edgecolors='none')

        plt.xlim(0, self.nkpoints-1)
        lim = 0.1
        plt.ylim(min(self.eigen.flatten()-lim)*RytoeV, max(self.eigen.flatten()+lim)*RytoeV)
        plt.show()

    def get_eigen(self):
        """ Return eigenvalues
        """
        datafile_xml = self.datafile_xml
        eigen = []
        for ik in xrange(self.nkpoints):
          eigen.append( map(float, self.datafile_xml.find("EIGENVALUES/K-POINT.%d/EIG"%(ik+1)).text.split() ))
        self.eigen = np.array(eigen)
        return eigen

    def write_proj(self,filename='proj'):
        """
        Write the projection array in a numpy file
        """
        np.savez(filename,proj=self.proj)
        
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
            s += "n -> iatom:%3d atype:%2s wfc:%d j:%s l:%s m:%s m_j:%s\n"%(state['iatom'],state['atype'],state['wfc'],str(state['j']),str(state['l']),str(state['m']),str(state['m_j']))
        return s
