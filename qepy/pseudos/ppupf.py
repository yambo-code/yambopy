#
# This file is part of yambopy
#
import xml.etree.ElementTree as ET
from re import findall
import numpy as np

HatoeV = 27.2107

class PPUPF():
    # This class reads pseudopotentials in UPF format
    """ Class to read a pseudopotential file in UPF format.

	So far only NC case implemented.
    """
    
    def __init__(self,filename):

	self.filename = filename

	self.read_datafile()

    def read_datafile(self):
	"""
	Read data from UPF file
	"""
	self.ppfile_upf = ET.parse( self.filename ).getroot()
