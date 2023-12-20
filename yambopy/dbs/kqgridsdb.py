# Copyright (c) 2022
# All rights reserved.
#
# This file is part of the yambopy project
#
import os
import numpy as np
from netCDF4 import Dataset
from yambopy.tools.string import marquee

class YamboBZgridsDB(object):
    """
    This class reads the database ndb.kindx. 
    It is mainly used to maintain consistency in the indices of the k/q-point
    expansions from IBZ to full BZ between yambo and yambopy.
    
    Actual kpoint coordinates are found as in YamboLatticeDB attributes.
    Actual qpoint coordinates are found as YamboEm1sDB, YamboElphonDB, YamboExcitonDB, YamboExcitonPhononDB attributes.
    
    For reference check src/bz_ops/bz_samp_indexes.F in the yambo source.
    NB: The -1 are needed only to match python and fortran indices in python.
    
    ! ikbz=(ik,is) --<--:--<-- okbz=(ok,os) = (IK-Q)
    !                   :
    !                  /:\ iqbz=(iq,is)
    !                   :
    !
    ! iq_is = ik_is-ok_os-Go
    !
    ! qindx_X(iq,ikbz,1)-1=okbz
    ! qindx_X(iq,ikbz,2)-1=iGo
    !
    ! qindx_B(okbz,ikbz,1)-1=iqbz
    ! qindx_B(okbz,ikbz,2)-1=iGo
    !
    ! qindx_S(ik,iqbz,1)-1=okbz
    ! qindx_S(ik,iqbz,2)-1=iGo
    !
    ! qindx_C(ikbz,iqbz,1)-1=okbz
    ! qindx_C(ikbz,iqbz,2)-1=iGo    

    NB: the indices included in the 
    """
    def __init__(self,filename):
            
        if not os.path.isfile(filename):
            raise FileNotFoundError("error opening %s in YamboKptsDB"%filename)

        with Dataset(filename) as database:

            if 'CH_GRIDS' in database.variables: grid_str = 'CH_GRIDS'
            elif 'GRIDS_CH' in database.variables: grid_str = 'GRIDS_CH'
            self.grid_types = database.variables[grid_str][0].tobytes().decode('utf-8').strip()
            
            if 'X' in self.grid_types: self.qindx_X =  database.variables['Qindx'][:].T
            if 'B' in self.grid_types: self.qindx_B =  database.variables['Bindx'][:].T
            if 'S' in self.grid_types: self.qindx_S =  database.variables['Sindx'][:].T
            if 'C' in self.grid_types: self.qindx_C =  database.variables['Cindx'][:].T

    def get_string(self,mark="="):
        lines = []; app = lines.append
        app( marquee(self.__class__.__name__,mark=mark) )
        app( "Tables:                   %s"%self.grid_types)
        if 'S' in self.grid_types:
            app( "number of kpoints (IBZ):  %s"%self.qindx_S.shape[0])
        if 'X' in self.grid_types:
            app( "number of kpoints (BZ):   %s"%self.qindx_X.shape[1])
            app( "number of qpoints (IBZ):  %s"%self.qindx_X.shape[0])    
        if 'S' in self.grid_types:
            app( "number of qpoints (BZ):   %s"%self.qindx_S.shape[1])
        return '\n'.join(lines)
    
    def __str__(self):
        return self.get_string()
