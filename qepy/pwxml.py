from __future__ import print_function
from __future__ import absolute_import
# Copyright (C) 2015 Henrique Pereira Coutada Miranda, Alejandro Molina-Sanchez
# All rights reserved.
#
# This file is part of yambopy
#
import xml.etree.ElementTree as ET
from   qepy.auxiliary import *
from   numpy import array
from   .lattice import *

HatoeV = 27.2107

class PwXML():
    """ Class to read data from a Quantum espresso XML file
    """
    _eig_xml  = 'eigenval.xml'

    def __init__(self,prefix,path='.'):
        """ Initlize the structure with the path where the datafile.xml is
        """
        self.prefix = prefix
        self.path   = path

        datafiles = {'data-file.xml':        self.read_datafile,
                     'data-file-schema.xml': self.read_datafile_schema}

        done_reading = False
        #check if the name is data-file.xml or data-file-schema.xml or whatever....
        for filename,read in datafiles.items():
            path_filename = "%s/%s.save/%s"%(path, prefix, filename)
            if os.path.isfile(path_filename):
                print("reading %s"%filename)
                done_reading = read(path_filename)
                break
        
        #trap errors
        if not done_reading:
            possible_files = " or ".join(datafiles.keys())
            raise ValueError('Failed to read %s in %s/%s.save'%(possible_files,path,prefix))

    def read_datafile(self,filename):
        """
        Read some data from the xml file in the old format of quantum espresso
        """
        self.datafile_xml = ET.parse( filename ).getroot()

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
 
        #get eigenvalues
        eigen = []
        for ik in xrange(self.nkpoints):
            for EIGENVALUES in ET.parse( "%s/%s.save/K%05d/%s" % (self.path,self.prefix,(ik + 1),self._eig_xml) ).getroot().findall("EIGENVALUES"):
                eigen.append(map(float, EIGENVALUES.text.split()))
        self.eigen  = eigen

        #get fermi
        self.fermi = float(self.datafile_xml.find("BAND_STRUCTURE_INFO/FERMI_ENERGY").text)

        return True 

    def read_datafile_schema(self,filename):
        """
        Read the data from the xml file in the new format of quantum espresso
        """
        self.datafile_xml = ET.parse( filename ).getroot()

        #get cell
        self.cell = []
        for i in xrange(1,4):
            cell_lat = self.datafile_xml.findall("input/atomic_structure/cell/a%d"%i)[0].text
            self.cell.append([float(x) for x in cell_lat.strip().split()])

        #calculate acell
        self.acell = [ np.linalg.norm(a) for a in self.cell ]

        #get reciprocal cell
        self.rcell = []
        for i in xrange(1,4):
            rcell_lat = self.datafile_xml.findall("output/basis_set/reciprocal_lattice/b%d"%i)[0].text
            self.rcell.append([float(x) for x in rcell_lat.strip().split()])

        #get atoms
        self.natoms = int(self.datafile_xml.findall("output/atomic_structure")[0].get('nat'))
        self.atoms = []
        atoms = self.datafile_xml.findall("output/atomic_structure/atomic_positions/atom")
        for i in xrange(self.natoms):
            atom = atoms[i].text
            self.atoms.append([float(x) for x in atom.strip().split()])

        #get nkpoints
        self.nkpoints = int(self.datafile_xml.findall("output/band_structure/nks")[0].text.strip())
        # Read the number of BANDS
        self.nbands = int(self.datafile_xml.findall("output/band_structure/nbnd")[0].text.strip())

        #get ks states
        kstates = self.datafile_xml.findall('output/band_structure/ks_energies')

        #get k-points
        self.kpoints = [] 
        for i in range(self.nkpoints):
            kpoint = [float(x) for x in kstates[i].findall('k_point')[0].text.strip().split()]
            self.kpoints.append( kpoint )

        #get eigenvalues
        self.eigen = []
        for k in range(self.nkpoints):
            eigen = [float(x) for x in kstates[k].findall('eigenvalues')[0].text.strip().split()]
            self.eigen.append( eigen )
        self.eigen = np.array(self.eigen)
 
        #get fermi
        self.fermi = float(self.datafile_xml.find("output/band_structure/highestOccupiedLevel").text)

        return True

    def get_scaled_positions(self):
        """ get the atomic positions in reduced coordinates
        """
        return car_red(self.atoms,self.cell)

    def __str__(self):
        s = ""
        s += "cell:\n"
        for c in self.cell:
            s += ("%12.8lf "*3)%tuple(c)+'\n'
        s += "atoms:\n"
        for a in self.atoms:
            s += ("%12.8lf "*3)%tuple(a)+'\n'
        s += "nkpoints: %d\n"%self.nkpoints
        s += "nbands:   %d\n"%self.nbands
        return s

    def plot_eigen(self,path=[],xlim=(),ylim=()):
        """ plot the eigenvalues using matplotlib
        """
        import matplotlib.pyplot as plt
        
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
        if fmt=='gnuplot':
            f = open('%s.dat'%self.prefix,'w')
            for ib in xrange(self.nbands):
                for ik in xrange(self.nkpoints):
                    f.write("%.1lf %.4lf \n " % (ik,self.eigen[ik][ib]*HatoeV) )
                f.write("\n")
            f.close()
        else:
            print('fmt %s not implemented'%fmt)
