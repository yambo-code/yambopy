#
# This file is part of yambopy
#
import xml.etree.ElementTree as ET
from re import findall
import numpy as np
from yambopy.tools.string import marquee

class PPUPF():
    """ 
    Class to read a pseudopotential file in UPF format.
   
    Reads basic info and the LOCAL and NONLOCAL grids. Both grids
    can be set to zero with the relevant function. A new pseudo file will be saved.
 
    So far only NC case implemented.
    """
    
    def __init__(self,filename):

        self.filename = filename
        self.ppfile_upf = ET.parse( filename ).getroot()

        self.read_header()

        if not self.pseudo_type=="NC": 
            raise NotImplementedError("PP is %s. Only NC pseudos implemented so far."%self.pseudo_type)

        # Mesh: PP_R and PP_RAB (not read)
        
        self.read_PP_local()

        self.read_PP_nonlocal()
        

    def read_header(self):
        """
        Read PP_header
        """
        XML_el = "PP_HEADER"

        self.element        = self.ppfile_upf.findall(XML_el)[0].get("element").strip()
        self.pseudo_type    = self.ppfile_upf.findall(XML_el)[0].get("pseudo_type").strip()
        self.number_of_proj = self.ppfile_upf.findall(XML_el)[0].get("number_of_proj").strip()

    def read_PP_local(self):
        """
        Read values of Vloc on PP mesh
        """
        XML_el = "PP_LOCAL"
   
        vloc_size = int(self.ppfile_upf.findall(XML_el)[0].get("size").strip())
        vloc_clmn = int(self.ppfile_upf.findall(XML_el)[0].get("columns").strip())
        self.vloc = self.ppfile_upf.findall(XML_el)[0].text

    def read_PP_nonlocal(self):
        """
        Read beta projectors (nonlocal part of PP)
        """
        XML_el = "PP_NONLOCAL/PP_BETA."        
        nproj = int(self.number_of_proj)

        KB_sizes, KB_clmns = [], []
        self.betas = []
        for iproj in range(nproj):
            KB_sizes.append( int(self.ppfile_upf.findall(XML_el+"%d"%(iproj+1))[0].get("size").strip()) )
            KB_clmns.append( int(self.ppfile_upf.findall(XML_el+"%d"%(iproj+1))[0].get("columns").strip()) )
            self.betas.append(  self.ppfile_upf.findall(XML_el+"%d"%(iproj+1))[0].text )

    def copy_pseudo_file(self,tag='copy'):
        """
        Copy pseudo file
        """
        import shutil

        if self.filename[-3:]=='upf': ext = '.upf'
        if self.filename[-3:]=='UPF': ext = '.UPF'
        new_pp_name = self.filename[:-4]+'_'+tag+ext

        shutil.copyfile(self.filename, new_pp_name)

        self.new_pp_name = new_pp_name
        new_ppupf = ET.parse(new_pp_name)

        return new_ppupf

    def modify_data_string(self,data_string):
        """
        Identify numbers in the data string and set them to zero
        """ 
        import re
        from copy import deepcopy
    
        new_string = deepcopy(data_string)

        new_string = re.sub('\d', '0', new_string)
        
        return new_string 

    def set_vloc_to_zero(self):
        """
        We remove the local part of the pseudo by setting the projectors
        to zero.

        Saved in 'filename_novloc.UPF'
        """
        new_upf = self.copy_pseudo_file(tag='novloc')

        XML_el = "PP_LOCAL"
        
        new_vloc = self.modify_data_string(self.vloc)
        new_upf.findall(XML_el)[0].text = new_vloc

        new_upf.write(self.new_pp_name)

    def set_betas_to_zero(self):
        """
        We remove the nonlocal part of the pseudo by setting the projectors
        to zero.

        Saved in 'filename_noproj.UPF'
        """
        new_upf = self.copy_pseudo_file(tag='noproj')

        XML_el = "PP_NONLOCAL/PP_BETA."
        nproj = int(self.number_of_proj)

        for iproj in range(nproj): 
            new_beta = self.modify_data_string(self.betas[iproj])
            new_upf.findall(XML_el+"%d"%(iproj+1))[0].text = new_beta

        new_upf.write(self.new_pp_name)

    def __str__(self):
        lines = []; app = lines.append
        app(marquee(self.__class__.__name__))
        app("Pseudopotential: %s"%self.filename)
        app("Element             : %s"%self.element)
        app("Pseudo type         : %s"%self.pseudo_type)
        app("Number of projectors: %s"%self.number_of_proj)

        return "\n".join(lines)

