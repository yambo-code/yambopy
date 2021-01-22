#
# This file is part of yambopy
#
import xml.etree.ElementTree as ET
import numpy as np

"""
List of functions to read xml data from .xml files
"""

def get_xml_attrib(xml_tree,tag,attrib,repeated=False):
    """ 
    Extract xml attribute
    
    If repeated = True, first occurrence of attribute is stored.
    Else, its last occurrence is stored.
    """
    attribute =  None
    for element in xml_tree.iter(tag=tag):
        attribute = element.attrib[attrib]
        if repeated: break

    if attribute is None: raise ValueError('xml attribute %s not found'%attrib)
    else: return attribute
    
def get_xml_data(xml_tree,tag,as_type=float,repeated=False):
    """ 
    Extract xml data as specified type
    
    If repeated = True, first occurrence of attribute is stored.
    Else, its last occurrence is stored.
    """    
    data = None
    for element in xml_tree.iter(tag=tag):
        data_as_str = element.text.split()
        if len(data_as_str)>1: data = [as_type(d_str) for d_str in data_as_str]
        else: data = as_type(data_as_str[0])
        if repeated: break
        
    if data is None: raise ValueError('xml data within tag %s not found'%tag)
    else: return data

def get_xml_nk_bands(xml_tree):
    """
    Function to specifically get kpoint (cartesian) coordinates and corresponding eigenvalues (in Hartree)
    """    
    k_points_car  = []
    k_eigenvalues = []
    k_occupations = []
    for ks_energies in xml_tree.iter(tag='ks_energies'):
        k_points_car.append(  get_xml_data(ks_energies,'k_point',as_type=float)     )
        k_eigenvalues.append( get_xml_data(ks_energies,'eigenvalues',as_type=float) )
        k_occupations.append( get_xml_data(ks_energies,'occupations',as_type=float) )

    k_points_car  = np.array(k_points_car)
    k_eigenvalues = np.array(k_eigenvalues)
    k_occupations = np.array(k_occupations)

    return k_points_car, k_eigenvalues, k_occupations
