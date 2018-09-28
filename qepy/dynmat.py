# Copyright (C) 2018 Henrique Pereira Coutada Miranda, Alejandro Molina Sanchez
# All rights reserved.
#
# This file is part of yambopy
#
from __future__ import print_function, division
import os
import re
from math import sqrt
import numpy as np
from qepy.auxiliary import *

eVtocm1 = 8065.54429
cm1toeV = 1.0/eVtocm1
ha2ev  = 27.211396132
ev2ha  = 1.0/ha2ev
Thz2cm1 = 33.35641
cm12Thz = 1.0/33.35641

__all__ = [
"DynmatIn",
]

class DynmatIn(object):
    """
    Generate an manipulate quantum espresso input files for matdyn.x
    """
    def __init__(self):
        self.variable = dict()
        self.qpoints     = [] 

    @classmethod
    def from_prefix(cls,prefix,asr='simple'):
        dm = DynmatIn()
        dm['asr'] = "'%s'"%asr
        dm['fildyn'] = "'%s.dyn1'"%prefix
        dm['filout'] = "'%s.modes'"%prefix
        return dm

    def write(self,filename):
        with open(filename,'w') as f:
            f.write(str(self))

    def __str__(self):
        s = '&input\n'
        s += self.stringify_group('',self.variable) #print variable
        if len(self.qpoints) > 0:
          s+=("%d\n"%len(self.qpoints))
          for q in self.qpoints:
            s+=("%12.8lf %12.8lf %12.8lf")%tuple(q[:3])+"\n"
        return s

    def __setitem__(self,key,value):
        self.variable[key] = value

    def __getitem__(self,key):
        return self.variable[key]

    def stringify_group(self, keyword, group):
        if group != {}:
            string='\n'
            for keyword in group:
                string += "%20s = %s\n" % (keyword, group[keyword])
            string += "/\n"
            return string
        else:
            return ''

