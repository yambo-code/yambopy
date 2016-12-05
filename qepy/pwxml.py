# Copyright (C) 2015 Henrique Pereira Coutada Miranda, Alejandro Molina-Sanchez
# All rights reserved.
#
# This file is part of yambopy
#
import xml.etree.ElementTree as ET
from   qepy.auxiliary import *
from   numpy import array

HatoeV = 27.2107

class PwXML():
    """ Class to read data from a Quantum espresso XML file
    """
    _datafile = 'data-file.xml'
    _eig_xml  = 'eigenval.xml'

    def __init__(self,prefix,path='.'):
        """ Initlize the structure with the path where the datafile.xml is
        """
        self.prefix = prefix
        self.path   = path
        self.datafile_xml = ET.parse( "%s/%s.save/%s"%(path, prefix, self._datafile)).getroot()

        #get acell
        self.acell = [ float(x) for x in self.datafile_xml.findall("CELL/CELL_DIMENSIONS")[0].text.strip().split('\n') ]

        #get cell
        self.cell = []
        for i in xrange(1,4):
            cell_lat = self.datafile_xml.findall("CELL/DIRECT_LATTICE_VECTORS/a%d"%i)[0].text
            self.cell.append([float(x) for x in cell_lat.strip().split()])

        #get reciprocal cell
        self.rcell = []
        for i in xrange(1,4):
            rcell_lat = self.datafile_xml.findall("CELL/RECIPROCAL_LATTICE_VECTORS/b%d"%i)[0].text
            self.rcell.append([float(x) for x in rcell_lat.strip().split()])

        #get atoms
        self.natoms = int(self.datafile_xml.findall("IONS/NUMBER_OF_ATOMS")[0].text)
        self.atoms = []
        for i in xrange(1,self.natoms+1):
            atom = self.datafile_xml.findall("IONS/ATOM.%d"%i)[0].get('tau')
            self.atoms.append([float(x) for x in atom.strip().split()])

        #get nkpoints
        self.nkpoints = int(self.datafile_xml.findall("BRILLOUIN_ZONE/NUMBER_OF_K-POINTS")[0].text.strip())
        # Read the number of BANDS
        self.nbands   = int(self.datafile_xml.find("BAND_STRUCTURE_INFO/NUMBER_OF_BANDS").text)

        #get k-points
        self.kpoints = [] 
        for i in range(self.nkpoints):
          k_aux = self.datafile_xml.findall('BRILLOUIN_ZONE/K-POINT.%d'%(i+1))[0].get('XYZ')
          self.kpoints.append([float(x) for x in k_aux.strip().split()])
 
        #get fermi
        self.fermi = float(self.datafile_xml.find("BAND_STRUCTURE_INFO/FERMI_ENERGY").text)
 
        self.eigen = None 

    def get_scaled_positions(self):
        """ get the atomic positions in reduced coordinates
        """
        return car_red(self.atoms,self.cell)

    def __str__(self):
        s = ""
        s += "cell:\n"
        for c in self.cell:
            s += str(c)+"\n"
        s += "atoms:\n"
        for a in self.atoms:
            s += str(a)+'\n'
        s += "nkpoints: %d\n"%self.nkpoints
        s += "nbands: %d\n"%self.nbands
        return s

    def plot_eigen(self,path=[],xlim=(),ylim=()):
        """ plot the eigenvalues using matplotlib
        """
        import matplotlib.pyplot as plt
        
        # Font selection and borders
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',serif="Computer Modern Roman",size=25)

        if self.eigen is None:
            self.get_eigen()
      
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

        #plot bands
        eigen = array(self.eigen)
        for ib in range(self.nbands):
           plt.plot(xrange(self.nkpoints),eigen[:,ib]*HatoeV - self.fermi*HatoeV, 'r-', lw=2)

        #plot options
        if xlim:
          plt.xlim(xlim)
        if ylim:
          plt.ylim(ylim)

        plt.show()

    def write_eigen(self,fmt='gnuplot'):
        """ write eigenvalues to a text file
        """
        if self.eigen is None:
            self.get_eigen()
        if fmt=='gnuplot':
            f = open('%s.dat'%prefix,'w')
            for ib in xrange(self.nbands):
                for ik in xrange(self.nkpoints):
                    f.write("%.1lf %.4lf \n " % (ik,self.eigen[ik][ib]*HatoeV) )
                f.write("\n")
            f.close()
        else:
            print 'fmt %s not implemented'%fmt

    def get_eigen(self):
        """ Return eigenvalues in Hartree
        """
        datafile_xml = self.datafile_xml
        eigen = []
        for ik in xrange(self.nkpoints):
            for EIGENVALUES in ET.parse( "%s/%s.save/K%05d/%s" % (self.path,self.prefix,(ik + 1),self._eig_xml) ).getroot().findall("EIGENVALUES"):
                eigen.append(map(float, EIGENVALUES.text.split()))
        self.eigen  = eigen
        return eigen

