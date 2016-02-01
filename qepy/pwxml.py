# Copyright (C) 2015 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
#
import xml.etree.ElementTree as ET
from pwpy.auxiliary import car_red

class PwXML():
    """ Class to read data from a Quantum espresso XML file
    """
    _datafile =  'data-file.xml'

    def __init__(self,prefix,path='.'):
        """ Initlize the structure with the path where the datafile.xml is
        """
        self.datafile_xml = ET.parse( "%s/%s.save/%s"%(path, prefix, self._datafile)).getroot()

        #get acell
        self.acell = [ float(x) for x in self.datafile_xml.findall("CELL/CELL_DIMENSIONS")[0].text.strip().split('\n') ]

        #get cell
        self.cell = []
        for i in xrange(1,4):
            cell_lat = self.datafile_xml.findall("CELL/DIRECT_LATTICE_VECTORS/a%d"%i)[0].text
            self.cell.append([float(x) for x in cell_lat.strip().split()])

        #get atoms
        self.natoms = int(self.datafile_xml.findall("IONS/NUMBER_OF_ATOMS")[0].text)
        self.atoms = []
        for i in xrange(1,self.natoms+1):
            atom = self.datafile_xml.findall("IONS/ATOM.%d"%i)[0].get('tau')
            self.atoms.append([float(x) for x in atom.strip().split()])

        #get nkpoints
        self.nkpoints = int(self.datafile_xml.findall("BRILLOUIN_ZONE/NUMBER_OF_K-POINTS")[0].text.strip())

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
        return s
